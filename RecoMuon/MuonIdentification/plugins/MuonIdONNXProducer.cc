// Event-wide inference from the raw features used by DeepMuonReco's Phase-2 model.
#include <array>
#include <atomic>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"

namespace {
  std::unique_ptr<cms::Ort::ONNXRuntime> makeRuntime(const std::string& modelPath, const std::string& backendName) {
    cms::Ort::Backend backend;
    if (backendName == "cpu")
      backend = cms::Ort::Backend::cpu;
    else if (backendName == "cuda")
      backend = cms::Ort::Backend::cuda;
    else
      throw cms::Exception("Configuration")
          << "MuonIdONNXProducer backend must be 'cpu' or 'cuda', got '" << backendName << "'";
    if (modelPath.empty())
      throw cms::Exception("Configuration") << "MuonIdONNXProducer requires a nonempty modelPath";
    try {
      auto options = cms::Ort::ONNXRuntime::defaultSessionOptions(backend);
      options.SetIntraOpNumThreads(1);
      options.SetInterOpNumThreads(1);
      options.SetExecutionMode(ORT_SEQUENTIAL);
      return std::make_unique<cms::Ort::ONNXRuntime>(modelPath, &options);
    } catch (const ::Ort::Exception& error) {
      throw cms::Exception("MuonIdONNX") << "Failed to initialize ONNX backend '" << backendName << "' with model '"
                                         << modelPath << "': " << error.what();
    }
  }

  template <typename Segment, typename Geometry>
  bool appendSegment(const Segment& segment, const Geometry& geometry, std::vector<float>& features) {
    const auto* detector = geometry.idToDet(segment.geographicalId());
    if (detector == nullptr)
      throw cms::Exception("MuonIdONNX") << "No geometry for segment detector " << segment.geographicalId().rawId();
    const auto position = detector->surface().toGlobal(segment.localPosition());
    const auto direction = detector->surface().toGlobal(segment.localDirection());
    const std::array<float, 6> values = {{float(position.x()),
                                          float(position.y()),
                                          float(position.z()),
                                          float(direction.x()),
                                          float(direction.y()),
                                          float(direction.z())}};
    for (float value : values) {
      if (!std::isfinite(value))
        return false;
    }
    features.insert(features.end(), values.begin(), values.end());
    return true;
  }
}  // namespace

class MuonIdONNXProducer : public edm::global::EDProducer<> {
public:
  explicit MuonIdONNXProducer(const edm::ParameterSet& config)
      : tracksToken_(consumes<reco::TrackCollection>(config.getParameter<edm::InputTag>("tracks"))),
        dtToken_(consumes<DTRecSegment4DCollection>(config.getParameter<edm::InputTag>("dtSegments"))),
        cscToken_(consumes<CSCSegmentCollection>(config.getParameter<edm::InputTag>("cscSegments"))),
        gemToken_(consumes<GEMSegmentCollection>(config.getParameter<edm::InputTag>("gemSegments"))),
        dtGeometryToken_(esConsumes()),
        cscGeometryToken_(esConsumes()),
        gemGeometryToken_(esConsumes()),
        minPt_(config.getParameter<double>("minPt")),
        minP_(config.getParameter<double>("minP")),
        maxAbsEta_(config.getParameter<double>("maxAbsEta")),
        runtime_(
            makeRuntime(config.getParameter<std::string>("modelPath"), config.getParameter<std::string>("backend"))) {
    produces<edm::ValueMap<float>>();
  }

  void produce(edm::StreamID, edm::Event& event, const edm::EventSetup& setup) const override {
    ++events_;
    const auto tracks = event.getHandle(tracksToken_);
    // -1 means outside the model's training selection and must never cause rejection.
    std::vector<float> scores(tracks->size(), -1.f);
    std::vector<unsigned int> indices;
    cms::Ort::FloatArrays inputs(8);
    bool finiteFeatures = true;
    for (unsigned int i = 0; i < tracks->size(); ++i) {
      const auto& track = tracks->at(i);
      // Match the ntuplizer's track_is_good_track mask, including its strict boundaries.
      if (!(track.pt() > minPt_ && track.p() > minP_ && std::abs(track.eta()) < maxAbsEta_))
        continue;
      const std::array<float, 7> values = {{float(track.px()),
                                            float(track.py()),
                                            float(track.pz()),
                                            float(track.vx()),
                                            float(track.vy()),
                                            float(track.vz()),
                                            float(track.charge())}};
      for (float value : values)
        finiteFeatures &= std::isfinite(value);
      inputs[0].insert(inputs[0].end(), values.begin(), values.end());
      indices.push_back(i);
    }

    if (!indices.empty()) {
      const auto dt = event.getHandle(dtToken_);
      const auto csc = event.getHandle(cscToken_);
      const auto gem = event.getHandle(gemToken_);
      if (!dt.isValid() || !csc.isValid() || !gem.isValid()) {
        ++missingCollectionEvents_;
        if (!warnedMissing_.exchange(true))
          edm::LogWarning("DeepMuonReco")
              << "MuonIdONNX: missing segment collection(s) treated as empty: DT=" << !dt.isValid()
              << " CSC=" << !csc.isValid() << " GEM=" << !gem.isValid();
      }
      if (dt.isValid() && dt->size() != 0) {
        const auto& geometry = setup.getData(dtGeometryToken_);
        for (const auto& segment : *dt)
          finiteFeatures &= appendSegment(segment, geometry, inputs[2]);
      }
      if (csc.isValid() && csc->size() != 0) {
        const auto& geometry = setup.getData(cscGeometryToken_);
        for (const auto& segment : *csc)
          finiteFeatures &= appendSegment(segment, geometry, inputs[4]);
      }
      if (gem.isValid() && gem->size() != 0) {
        const auto& geometry = setup.getData(gemGeometryToken_);
        for (const auto& segment : *gem) {
          // The ntuplizer stores only the ME0 subset in its gem_seg_* branches.
          if (GEMDetId(segment.geographicalId()).isME0())
            finiteFeatures &= appendSegment(segment, geometry, inputs[6]);
        }
      }

      if (finiteFeatures) {
        static const std::vector<std::string> inputNames = {"tracker_track",
                                                            "tracker_track_data_mask",
                                                            "dt_segment",
                                                            "dt_segment_data_mask",
                                                            "csc_segment",
                                                            "csc_segment_data_mask",
                                                            "gem_segment",
                                                            "gem_segment_data_mask"};
        std::vector<std::vector<int64_t>> shapes;
        for (unsigned int object = 0; object < 4; ++object) {
          const int64_t dim = object == 0 ? 7 : 6;
          auto& features = inputs[2 * object];
          auto& mask = inputs[2 * object + 1];
          mask.assign(features.size() / dim, 1.f);
          // ONNX attention accepts a masked placeholder for a genuinely empty collection.
          if (mask.empty()) {
            features.assign(dim, 0.f);
            mask.assign(1, 0.f);
          }
          const auto count = static_cast<int64_t>(mask.size());
          shapes.push_back({1, count, dim});
          shapes.push_back({1, count});
        }
        const auto outputs = runtime_->run(inputNames, inputs, shapes, {"logits"});
        ++inferenceEvents_;
        if (outputs.size() != 1 || outputs[0].size() != indices.size())
          throw cms::Exception("MuonIdONNX") << "Model logits do not match the selected track count";
        evaluatedTracks_ += indices.size();
        for (unsigned int i = 0; i < indices.size(); ++i) {
          const float logit = outputs[0][i];
          if (!std::isfinite(logit)) {
            ++nonfiniteScores_;
            warnNonfinite();
            continue;
          }
          const float expValue = std::exp(-std::abs(logit));
          scores[indices[i]] = logit >= 0.f ? 1.f / (1.f + expValue) : expValue / (1.f + expValue);
        }
      } else {
        ++nonfiniteInputEvents_;
        warnNonfinite();
      }
    }

    auto output = std::make_unique<edm::ValueMap<float>>();
    edm::ValueMap<float>::Filler filler(*output);
    filler.insert(tracks, scores.begin(), scores.end());
    filler.fill();
    event.put(std::move(output));
  }

  void endJob() override {
    edm::LogVerbatim("DeepMuonReco") << "MuonIdONNX inference: events=" << events_.load()
                                     << " inferenceEvents=" << inferenceEvents_.load()
                                     << " evaluatedTracks=" << evaluatedTracks_.load()
                                     << " nonfiniteScores=" << nonfiniteScores_.load()
                                     << " nonfiniteInputEvents=" << nonfiniteInputEvents_.load()
                                     << " missingCollectionEvents=" << missingCollectionEvents_.load();
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
    desc.add<edm::InputTag>("dtSegments", edm::InputTag("dt4DSegments"));
    desc.add<edm::InputTag>("cscSegments", edm::InputTag("cscSegments"));
    desc.add<edm::InputTag>("gemSegments", edm::InputTag("gemSegments"));
    desc.add<std::string>("modelPath", "");
    desc.add<std::string>("backend", "cpu");
    desc.add<double>("minPt", 0.5);
    desc.add<double>("minP", 2.5);
    desc.add<double>("maxAbsEta", 3.0);
    descriptions.add("muonIdONNX", desc);
  }

private:
  void warnNonfinite() const {
    if (!warnedNonfinite_.exchange(true))
      edm::LogWarning("DeepMuonReco") << "MuonIdONNX: nonfinite input or output; affected tracks bypass the filter";
  }

  const edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  const edm::EDGetTokenT<DTRecSegment4DCollection> dtToken_;
  const edm::EDGetTokenT<CSCSegmentCollection> cscToken_;
  const edm::EDGetTokenT<GEMSegmentCollection> gemToken_;
  const edm::ESGetToken<DTGeometry, MuonGeometryRecord> dtGeometryToken_;
  const edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeometryToken_;
  const edm::ESGetToken<GEMGeometry, MuonGeometryRecord> gemGeometryToken_;
  const double minPt_, minP_, maxAbsEta_;
  const std::unique_ptr<cms::Ort::ONNXRuntime> runtime_;
  mutable std::atomic<unsigned long long> events_{0}, inferenceEvents_{0}, evaluatedTracks_{0}, nonfiniteScores_{0},
      nonfiniteInputEvents_{0}, missingCollectionEvents_{0};
  mutable std::atomic<bool> warnedMissing_{false}, warnedNonfinite_{false};
};

DEFINE_FWK_MODULE(MuonIdONNXProducer);
