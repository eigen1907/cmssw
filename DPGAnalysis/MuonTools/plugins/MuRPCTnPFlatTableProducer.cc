/** \class MuRPCTnPFlatTableProducer MuRPCTnPFlatTableProducer.cc DPGAnalysis/MuonTools/plugins/MuRPCTnPFlatTableProducer.cc
 *
 * Helper class : the RPC Tag and probe efficiency  filler
 *
 * adapted from
 * https://gitlab.cern.ch/jhgoh/RPCDPGAnalysis/-/blob/855eb87b/SegmentAndTrackOnRPC/plugins/RPCTrackerMuonProbeProducer.cc
 * https://gitlab.cern.ch/jhgoh/RPCDPGAnalysis/-/blob/15977fcf/SegmentAndTrackOnRPC/plugins/MuonHitFromTrackerMuonAnalyzer.cc
 *
 * \author Seungjin Yang (Kyung Hee U.)
 */
#include "DPGAnalysis/MuonTools/interface/MuBaseFlatTableProducer.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"

#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"

#include "DataFormats/Math/interface/deltaR.h"


namespace rpctnp {
  enum class MuonIdType {
    kTight = 0,
    kSoft,
    kTracker,
    kRPC,
    kNone,
  };

  struct Result {
    reco::MuonRef tag;
    reco::MuonRef probe;
    std::vector<std::pair<reco::MuonChamberMatch, const RPCRecHit*> > measurements;
  };
}


class MuRPCTnPFlatTableProducer : public MuBaseFlatTableProducer {
public:
  MuRPCTnPFlatTableProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);

protected:
  void fillTable(edm::Event&) final;
  void getFromES(const edm::Run&, const edm::EventSetup&) final;
  void getFromES(const edm::EventSetup&) final;

private:
  std::set<std::string> findFilterModules(const edm::Event&,
                                    const edm::Handle<edm::TriggerResults>&);

  std::vector<math::XYZTLorentzVector> findTriggerObjectsMomenta(
      const edm::Handle<trigger::TriggerEvent>&,
      const std::set<std::string>&);

  reco::MuonRef findTagMuon(const edm::Handle<reco::MuonCollection>&,
                            const reco::Vertex&,
                            const std::vector<math::XYZTLorentzVector>&);

  reco::MuonRef findProbeMuon(const edm::Handle<reco::MuonCollection>&,
                              const reco::Vertex&,
                              const reco::MuonRef&);

  bool checkIfMuonIsOutsideRoll(const reco::MuonChamberMatch&);

  std::pair<const RPCRecHit*, double> findClosestHit(
      const reco::MuonChamberMatch&,
      const edm::Handle<RPCRecHitCollection>&);

  rpctnp::Result performTagAndProbe(const edm::Event&);

  // utility member functions
  rpctnp::MuonIdType convertStrToMuonIdType(std::string);

  bool checkIfMuonPassId(const reco::Muon&,
                   const reco::Vertex&,
                   const rpctnp::MuonIdType&);

  double computeMuonPFRelIso(const reco::MuonPFIsolation&,
                             const double);

  bool doTriggerMatching(const reco::Muon&,
                         const std::vector<math::XYZTLorentzVector>&);

  bool isFiducial(const reco::MuonChamberMatch&);

  // member data
  const nano_mu::EDTokenHandle<reco::VertexCollection> m_primary_vertex_collection_token;
  const nano_mu::EDTokenHandle<edm::TriggerResults> m_trigger_results_token;
  const nano_mu::EDTokenHandle<trigger::TriggerEvent> m_trigger_event_token;
  const nano_mu::EDTokenHandle<reco::MuonCollection> m_muon_collection_token;
  const nano_mu::EDTokenHandle<RPCRecHitCollection> m_rpc_rec_hit_collection_token;
  //// tag muon
  const double m_tag_muon_min_pt;
  const double m_tag_muon_max_abs_eta;
  const rpctnp::MuonIdType m_tag_muon_id_type;
  const double m_tag_muon_max_rel_iso;
  const std::vector<std::string> m_tag_muon_trigger_matching_triggers;
  const double m_tag_muon_trigger_matching_max_delta_r;
  const double m_tag_muon_trigger_matching_max_delta_pt_rel;
  //// probe muon
  const double m_probe_muon_min_pt;
  const double m_probe_muon_max_abs_eta;
  const rpctnp::MuonIdType m_probe_muon_id_type;
  //// dimuon
  const bool m_dimuon_check_opposite_sign;
  const double m_dimuon_min_delta_r;
  const double m_dimuon_min_mass;
  const double m_dimuon_max_mass;

  nano_mu::ESTokenHandle<RPCGeometry, MuonGeometryRecord, edm::Transition::BeginRun> m_rpc_geometry;
  HLTConfigProvider m_hlt_config_provider;

  enum class Column {
    // tag
    kTagPt = 0,
    kTagEta,
    kTagPhi,
    // probe
    kProbePt,
    kProbeEta,
    kProbePhi,
    kProbeTime,
    kProbeDXDZ,
    kProbeDYDZ,
    // dimuon
    kDimuonPt,
    kDimuonMass,
    // rpc
    kRegion,
    kRing,
    kStation,
    kSector,
    kLayer,
    kSubsector,
    kRoll,
    // flags
    kIsFiducial,
    kIsMatched,
    // hit
    kClusterSize,
    kBunchX,
    //
    kResidualX,
    kResidualY,
    kPullX,
    kPullY,
    kPullXV2,
    kPullYV2,
  };


  const std::vector<std::tuple<Column, std::string, std::string> > FLOAT_COLUMNS = {
    // tag muon
    {Column::kTagPt, "tag_pt", "Tag Muon pT [GeV]"},
    {Column::kTagEta, "tag_eta", "Tag Muon eta"},
    {Column::kTagPhi, "tag_phi", "Tag Muon phi [rad]"},
    // probe
    {Column::kProbePt, "probe_pt", "Probe Muon pT [GeV]"},
    {Column::kProbeEta, "probe_eta", "Probe Muon eta"},
    {Column::kProbePhi, "probe_phi", "Probe Muon eta"},
    {Column::kProbeTime, "probe_time", "Probe Muon time [s]"},
    {Column::kProbeDXDZ, "probe_dxdz", "Probe Muon dX/dZ"},
    {Column::kProbeDYDZ, "probe_dydz", "Probe Muon dY/dZ"},
    // dimuon
    {Column::kDimuonPt, "dimuon_pt", "Tag+Probe pT [GeV]"},
    {Column::kDimuonMass, "dimuon_mass", "Tag+Probe mass [GeV]"},
    // hit
    {Column::kResidualX, "residual_x", "Residual X [cm]"},
    {Column::kResidualY, "residual_y", "Residual Y [cm]"},
    {Column::kPullX, "pull_x", "Pull X"},
    {Column::kPullY, "pull_y", "Pull Y"},
    {Column::kPullXV2, "pull_x_v2", "Pull X V2"},
    {Column::kPullYV2, "pull_y_v2", "Pull Y V2"},
  };

  const std::vector<std::tuple<Column, std::string, std::string> > INT8_COLUMNS = {
    {Column::kRegion, "region", "region"},
    {Column::kRing, "ring", "ring"},
    {Column::kStation, "station", "station"},
    {Column::kSector, "sector", "sector"},
    {Column::kLayer, "layer", "layer"},
    {Column::kSubsector, "subsector", "subsector"},
    {Column::kRoll, "roll", "roll"},
  };

  const std::vector<std::tuple<Column, std::string, std::string> > INT_COLUMNS = {
    //
    {Column::kClusterSize, "cls", "Cluster Size"},
    {Column::kBunchX, "bx", "Bunch Crossing"},
  };

  const std::vector<std::tuple<Column, std::string, std::string> > BOOL_COLUMNS = {
    {Column::kIsFiducial, "is_fiducial", "is a probe muon in the fiducial area"},
    {Column::kIsMatched, "is_matched", "is matched"},
  };
};


MuRPCTnPFlatTableProducer::MuRPCTnPFlatTableProducer(const edm::ParameterSet& config)
    : MuBaseFlatTableProducer(config),
      m_primary_vertex_collection_token{config, consumesCollector(), "primaryVertexCollectionTag"},
      m_trigger_results_token{config, consumesCollector(), "triggerResultsTag"},
      m_trigger_event_token{config, consumesCollector(), "triggerEventTag"},
      m_muon_collection_token{config, consumesCollector(), "muonCollectionTag"},
      m_rpc_rec_hit_collection_token{config, consumesCollector(), "rpcRecHitCollectionTag"},
      // tag muon
      m_tag_muon_min_pt{config.getParameter<double>("tagMuonMinPt")},
      m_tag_muon_max_abs_eta{config.getParameter<double>("tagMuonMaxAbsEta")},
      m_tag_muon_id_type{convertStrToMuonIdType(config.getParameter<std::string>("tagMuonIdType"))},
      m_tag_muon_max_rel_iso{config.getParameter<double>("tagMuonMaxRelIso")},
      m_tag_muon_trigger_matching_triggers{config.getParameter<std::vector<std::string> >("tagMuonTriggerMatchingPaths")},
      m_tag_muon_trigger_matching_max_delta_r{config.getParameter<double>("tagMuonTriggerMatchingMaxDeltaR")},
      m_tag_muon_trigger_matching_max_delta_pt_rel{config.getParameter<double>("tagMuonTriggerMatchingMaxDeltaPtRel")},
      // probe muon
      m_probe_muon_min_pt{config.getParameter<double>("probeMuonMinPt")},
      m_probe_muon_max_abs_eta{config.getParameter<double>("probeMuonMaxAbsEta")},
      m_probe_muon_id_type{convertStrToMuonIdType(config.getParameter<std::string>("probeMuonIdType"))},
      // dimuon
      m_dimuon_check_opposite_sign{config.getParameter<bool>("dimuonCheckOppositeSign")},
      m_dimuon_min_delta_r{config.getParameter<double>("dimuonMinDeltaR")},
      m_dimuon_min_mass{config.getParameter<double>("dimuonMinMass")},
      m_dimuon_max_mass{config.getParameter<double>("dimuonMaxMass")},
      // non-constant
      m_rpc_geometry{consumesCollector()} {
  produces<nanoaod::FlatTable>();
}

void MuRPCTnPFlatTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  // MuBaseFlatTableProducer
  desc.add<std::string>("name", "rpcTnP");
  // MuRPCTnPFlatTableProducer
  //// input tags
  desc.add<edm::InputTag>("primaryVertexCollectionTag", edm::InputTag{"offlinePrimaryVertices"});
  desc.add<edm::InputTag>("triggerResultsTag", edm::InputTag{"TriggerResults::HLT"});
  desc.add<edm::InputTag>("triggerEventTag", edm::InputTag{"hltTriggerSummaryAOD::HLT"});
  desc.add<edm::InputTag>("muonCollectionTag", edm::InputTag{"muons"});
  desc.add<edm::InputTag>("rpcRecHitCollectionTag", edm::InputTag{"rpcRecHits"});
  //// tag muon
  desc.add<double>("tagMuonMinPt", 30.0);
  desc.add<double>("tagMuonMaxAbsEta", 2.4);
  desc.add<std::string>("tagMuonIdType", "tight");
  desc.add<double>("tagMuonMaxRelIso", 0.25);
  desc.add<std::vector<std::string> >("tagMuonTriggerMatchingPaths",
                                      {"HLT_IsoMu24", "HLT_IsoMu27", "HLT_Mu50", "HLT_Mu55"});
  desc.add<double>("tagMuonTriggerMatchingMaxDeltaR", 0.1);
  desc.add<double>("tagMuonTriggerMatchingMaxDeltaPtRel", 0.5);
  //// probe muon
  desc.add<double>("probeMuonMinPt", 10.0);
  desc.add<double>("probeMuonMaxAbsEta", 2.1);
  desc.add<std::string>("probeMuonIdType", "tracker");
  //// muon pair
  desc.add<bool>("dimuonCheckOppositeSign", true);
  desc.add<double>("dimuonMinDeltaR", 0.1);
  desc.add<double>("dimuonMinMass", 70.0);
  desc.add<double>("dimuonMaxMass", 110.0);

  desc.setAllowAnything();
  descriptions.addWithDefaultLabel(desc);
}

void MuRPCTnPFlatTableProducer::getFromES(const edm::Run& run, const edm::EventSetup& environment) {
  m_rpc_geometry.getFromES(environment);

  bool changed = true;
  m_hlt_config_provider.init(run, environment, "HLT", changed);
}

void MuRPCTnPFlatTableProducer::getFromES(const edm::EventSetup& environment) {
}

void MuRPCTnPFlatTableProducer::fillTable(edm::Event& ev) {
  /////////////////////////////////////////////////////////////////////////////
  // STEP 1: TAG AND PROBE
  /////////////////////////////////////////////////////////////////////////////
  const rpctnp::Result result = performTagAndProbe(ev);
  const size_t size = result.measurements.size();

  /////////////////////////////////////////////////////////////////////////////
  // STEP 2: CREATE COLUMNS
  /////////////////////////////////////////////////////////////////////////////
  std::map<Column, std::vector<float> > float_columns;
  std::map<Column, std::vector<int8_t> > int8_columns;
  std::map<Column, std::vector<int> > int_columns;
  std::map<Column, std::vector<bool> > bool_columns;

  for (const auto& [key, name, descr] : FLOAT_COLUMNS) {
    float_columns.insert({key, {}});
  }
  for (const auto& [key, name, descr] : INT8_COLUMNS) {
    int8_columns.insert({key, {}});
  }
  for (const auto& [key, name, descr] : INT_COLUMNS) {
    int_columns.insert({key, {}});
  }
  for (const auto& [key, name, descr] : BOOL_COLUMNS) {
    bool_columns.insert({key, {}});
  }

  for (auto& [key, value] : float_columns) {
    value.reserve(size);
  }
  for (auto& [key, value] : int8_columns) {
    value.reserve(size);
  }
  for (auto& [key, value] : int_columns) {
    value.reserve(size);
  }
  for (auto& [key, value] : bool_columns) {
    value.reserve(size);
  }

  if (result.probe.isNonnull()) {
    float_columns.at(Column::kTagPt).assign(size, result.tag->pt());
    float_columns.at(Column::kTagEta).assign(size, result.tag->eta());
    float_columns.at(Column::kTagPhi).assign(size, result.tag->phi());

    float_columns.at(Column::kProbePt).assign(size, result.probe->pt());
    float_columns.at(Column::kProbeEta).assign(size, result.probe->eta());
    float_columns.at(Column::kProbePhi).assign(size, result.probe->phi());
    float_columns.at(Column::kProbeTime).assign(size, result.probe->time().timeAtIpInOut);

    const math::XYZTLorentzVector dimuon = result.tag->p4() + result.probe->p4();
    float_columns.at(Column::kDimuonPt).assign(size, dimuon.pt());
    float_columns.at(Column::kDimuonMass).assign(size, dimuon.mass());
  }

  for (const auto& [muon_chamber_match, hit] : result.measurements) {
    const RPCDetId det_id{muon_chamber_match.id};
    const RPCGeomServ geom_serv{det_id};

    float_columns.at(Column::kProbeDXDZ).push_back(muon_chamber_match.dXdZ);
    float_columns.at(Column::kProbeDYDZ).push_back(muon_chamber_match.dYdZ);

    // RPC Detector Information
    int8_columns.at(Column::kRegion).push_back(det_id.region());
    int8_columns.at(Column::kRing).push_back(det_id.ring());
    int8_columns.at(Column::kStation).push_back(det_id.station());
    int8_columns.at(Column::kSector).push_back(det_id.sector());
    int8_columns.at(Column::kLayer).push_back(det_id.layer());
    int8_columns.at(Column::kSubsector).push_back(det_id.subsector());
    int8_columns.at(Column::kRoll).push_back(det_id.roll());

    // Hit Information
    bool_columns.at(Column::kIsFiducial).push_back(isFiducial(muon_chamber_match));
    bool_columns.at(Column::kIsMatched).push_back(hit != nullptr);

    if (hit) {
      const LocalPoint probe_pos{muon_chamber_match.x, muon_chamber_match.y, 0.f};
      const LocalPoint hit_pos = hit->localPosition();
      const LocalError hit_pos_err = hit->localPositionError();
      //
      const float err_x = std::sqrt(hit_pos_err.xx() + std::pow(muon_chamber_match.xErr, 2));
      const float err_y = std::sqrt(hit_pos_err.yy() + std::pow(muon_chamber_match.yErr, 2));

      const float residual_x = hit_pos.x() - probe_pos.x();
      const float residual_y = hit_pos.y() - probe_pos.y();
      const float pull_x = residual_x / std::sqrt(hit_pos_err.xx());
      const float pull_y = residual_y / std::sqrt(hit_pos_err.yy());
      const float pull_x_v2 = residual_x / err_x;
      const float pull_y_v2 = residual_y / err_y;

      int_columns.at(Column::kClusterSize).push_back(hit->clusterSize());
      int_columns.at(Column::kBunchX).push_back(hit->BunchX());
      float_columns.at(Column::kResidualX).push_back(residual_x);
      float_columns.at(Column::kResidualY).push_back(residual_y);
      float_columns.at(Column::kPullX).push_back(pull_x);
      float_columns.at(Column::kPullY).push_back(pull_y);
      float_columns.at(Column::kPullXV2).push_back(pull_x_v2);
      float_columns.at(Column::kPullYV2).push_back(pull_y_v2);

    } else {
      int_columns.at(Column::kClusterSize).push_back(DEFAULT_DOUBLE_VAL);
      int_columns.at(Column::kBunchX).push_back(DEFAULT_DOUBLE_VAL);
      float_columns.at(Column::kResidualX).push_back(DEFAULT_DOUBLE_VAL);
      float_columns.at(Column::kResidualY).push_back(DEFAULT_DOUBLE_VAL);
      float_columns.at(Column::kPullX).push_back(DEFAULT_DOUBLE_VAL);
      float_columns.at(Column::kPullY).push_back(DEFAULT_DOUBLE_VAL);
      float_columns.at(Column::kPullXV2).push_back(DEFAULT_DOUBLE_VAL);
      float_columns.at(Column::kPullYV2).push_back(DEFAULT_DOUBLE_VAL);
    }
  } // measurements

  /////////////////////////////////////////////////////////////////////////////
  // STEP4: Fill Table
  /////////////////////////////////////////////////////////////////////////////
  auto table = std::make_unique<nanoaod::FlatTable>(size, m_name, false, false);
  table->setDoc("RPC Tag & Probe segment efficiency  information");

  for (const auto& [key, name, descr] : FLOAT_COLUMNS) {
    const auto& vec = float_columns.at(key);
    addColumn(table, name, vec, descr);
  }
  for (const auto& [key, name, descr] : INT8_COLUMNS) {
    const auto& vec = int8_columns.at(key);
    addColumn(table, name, vec, descr);
  }
  for (const auto& [key, name, descr] : INT_COLUMNS) {
    const auto& vec = int_columns.at(key);
    addColumn(table, name, vec, descr);
  }
  for (const auto& [key, name, descr] : BOOL_COLUMNS) {
    const auto& vec = bool_columns.at(key);
    addColumn(table, name, vec, descr);
  }

  ev.put(std::move(table));
}


// Returns a set of favourite and fired HLT filter module names.
// See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideHLTAnalysis
// TODO support wildcards like "HLT_IsoMu2*"
std::set<std::string> MuRPCTnPFlatTableProducer::findFilterModules(const edm::Event& event,
                                                             const edm::Handle<edm::TriggerResults>& trigger_results_handle) {

  // FIXME rename
  std::set<std::string> module_set;

  const std::vector<std::string>& trigger_vec = event.triggerNames(*trigger_results_handle).triggerNames();

  for (size_t idx = 0; idx < trigger_vec.size(); ++idx) {
    if (not trigger_results_handle->accept(idx) ) {
      continue;
    }
    const std::string trigger = m_hlt_config_provider.removeVersion(trigger_vec.at(idx));
    const std::vector<std::string>& module_vec = m_hlt_config_provider.saveTagsModules(idx);

    for (size_t favourite_idx = 0; favourite_idx < m_tag_muon_trigger_matching_triggers.size(); ++favourite_idx) {
      const std::string favourite_trigger = m_tag_muon_trigger_matching_triggers.at(favourite_idx);
      if (trigger == favourite_trigger) {
        module_set.insert(module_vec.begin(), module_vec.end());
      }
    } // target trigger
  } // trigger
  return module_set;
}

// Returns trigger objects' four momenta
std::vector<math::XYZTLorentzVector>
MuRPCTnPFlatTableProducer::findTriggerObjectsMomenta(
    const edm::Handle<trigger::TriggerEvent>& trigger_event_handle,
    const std::set<std::string>& filter_module_set) {

  std::vector<math::XYZTLorentzVector> momenta;

  const trigger::TriggerObjectCollection& trigger_object_collection = trigger_event_handle->getObjects();

  for (trigger::size_type module_idx = 0; module_idx < trigger_event_handle->sizeFilters(); ++module_idx) {
    const std::string& filter_label = trigger_event_handle->filterLabel(module_idx);

    if (filter_module_set.count(filter_label) == 0) {
      // TODO LogDebug
      continue;
    }

    for (const trigger::size_type key : trigger_event_handle->filterKeys(module_idx)) {
      momenta.push_back(trigger_object_collection.at(key).particle().p4());
    }
  }

  return momenta;
}

// Returns a tag muon
reco::MuonRef MuRPCTnPFlatTableProducer::findTagMuon(
    const edm::Handle<reco::MuonCollection>& muon_collection_handle,
    const reco::Vertex& primary_vertex,
    const std::vector<math::XYZTLorentzVector>& trigger_object_momentum_vec) {

  reco::MuonRef tag_muon;
  for (size_t muon_idx = 0; muon_idx < muon_collection_handle->size(); ++muon_idx) {
    const auto& mu = muon_collection_handle->at(muon_idx);
    // Basic kinematic cuts
    if (mu.pt() < m_tag_muon_min_pt) continue;
    if (std::abs(mu.eta()) > m_tag_muon_max_abs_eta) continue;
    // Muon ID
    if (mu.track().isNull()) continue;
    if (not checkIfMuonPassId(mu, primary_vertex, m_tag_muon_id_type)) continue;
    // Tight muon isolation
    if (computeMuonPFRelIso(mu.pfIsolationR03(), mu.pt()) > m_tag_muon_max_rel_iso)  continue;

    if (not doTriggerMatching(mu, trigger_object_momentum_vec)) {
      // FIXME LogDebug
      // why break?
      break;
    }

    if (tag_muon.isNull() or (mu.pt() > tag_muon->pt())) {
      tag_muon = reco::MuonRef(muon_collection_handle, muon_idx);
    }
  }
  return tag_muon;
}


// Returns a probe muon
reco::MuonRef MuRPCTnPFlatTableProducer::findProbeMuon(
    const edm::Handle<reco::MuonCollection>& muon_collection_handle,
    const reco::Vertex& primary_vertex,
    const reco::MuonRef& tag_muon) {

  reco::MuonRef probe_muon;

  for (size_t muon_idx = 0; muon_idx < muon_collection_handle->size(); ++muon_idx) {
    if (muon_idx == tag_muon.index()) continue;

    const auto& mu = muon_collection_handle->at(muon_idx);

    if (mu.pt() < m_probe_muon_min_pt) continue;
    if (std::abs(mu.eta()) > m_probe_muon_max_abs_eta) continue;
    if (not checkIfMuonPassId(mu, primary_vertex, m_probe_muon_id_type)) continue;
    if (mu.track()->originalAlgo() == reco::TrackBase::muonSeededStepOutIn) continue;

    // dimuon
    if (m_dimuon_check_opposite_sign and (mu.charge() == tag_muon->charge())) continue;
    if (deltaR(mu, tag_muon->p4()) < m_dimuon_min_delta_r) continue;
    const double dimuon_mass = (tag_muon->p4() + mu.p4()).mass();
    if (dimuon_mass < m_dimuon_min_mass) continue;
    if (dimuon_mass > m_dimuon_max_mass) continue;

    if (probe_muon.isNull() or (mu.pt() > probe_muon->pt())) {
      probe_muon = reco::MuonRef(muon_collection_handle, muon_idx);
    }
  }

  return probe_muon;
}


// Returns a boolean indicating whether the probe muon is inside the area of
// RPCRoll. This method assumes that reco::MuonChamberMatch::detector() == MuonSubdetId::RPC.
bool MuRPCTnPFlatTableProducer::checkIfMuonIsOutsideRoll(const reco::MuonChamberMatch& muon_chamber_match) {
  const RPCDetId det_id{muon_chamber_match.id};
  const RPCRoll* roll = m_rpc_geometry->roll(det_id);

  if (roll == nullptr) {
    // FIXME
    edm::LogError("") << "RPCRoll is nullptr: " << det_id;
    return false;
  }

  const LocalPoint muon_pos{muon_chamber_match.x, muon_chamber_match.y, 0.};
  return not roll->surface().bounds().inside(muon_pos);
}


// Returns the closest RPCRecHit to the probe muon passing through RPCRoll and
// the distance between them.
// The distance is defined as the residual x in the local x coordinate.
// TODO delta phi?
std::pair<const RPCRecHit*, double>
MuRPCTnPFlatTableProducer::findClosestHit(const reco::MuonChamberMatch& muon_chamber_match,
                                   const edm::Handle<RPCRecHitCollection>& rpc_rec_hit_collection_handle) {
  const RPCDetId rpc_id{muon_chamber_match.id};
  const RPCRoll* roll = m_rpc_geometry->roll(rpc_id);
  if (roll == nullptr) {
    // FIXME
    edm::LogError("") << "RPCRoll is nullptr: " << rpc_id;
    return std::make_pair(nullptr, -1.0);
  }

  const LocalPoint muon_pos{muon_chamber_match.x, muon_chamber_match.y, 0.};
  if (not roll->surface().bounds().inside(muon_pos)) {
    return std::make_pair(nullptr, -2.0);
  }

  const double muon_x = roll->centreOfStrip(roll->strip(muon_pos)).x();

  const RPCRecHit* matched_hit = nullptr;
  double min_distance = 1e6; // large number

  const RPCRecHitCollection::range& hit_range = rpc_rec_hit_collection_handle->get(rpc_id);
  for (edm::OwnVector<RPCRecHit>::const_iterator hit = hit_range.first; hit != hit_range.second; ++hit) {
    const double hit_x = hit->localPosition().x();
    const double distance = std::abs(hit_x - muon_x);
    if (distance < min_distance) {
      matched_hit = &(*hit);
      min_distance = distance;
    }
  } // hit loop
  return std::make_pair(matched_hit, min_distance);
}


// Returns a tag-and-probe result
rpctnp::Result MuRPCTnPFlatTableProducer::performTagAndProbe(const edm::Event& event) {
  rpctnp::Result result;

  const edm::Handle<reco::VertexCollection> primary_vertex_collection_handle = m_primary_vertex_collection_token.conditionalGet(event);
  const edm::Handle<edm::TriggerResults> trigger_results_handle = m_trigger_results_token.conditionalGet(event);
  const edm::Handle<trigger::TriggerEvent> trigger_event_handle = m_trigger_event_token.conditionalGet(event);
  const edm::Handle<reco::MuonCollection> muon_collection_handle = m_muon_collection_token.conditionalGet(event);
  const edm::Handle<RPCRecHitCollection> rpc_rec_hit_collection_handle = m_rpc_rec_hit_collection_token.conditionalGet(event);

  if (not primary_vertex_collection_handle.isValid()) return result;
  if (not trigger_results_handle.isValid()) return result;
  if (not trigger_event_handle.isValid()) return result;
  if (not muon_collection_handle.isValid()) return result;
  if (not rpc_rec_hit_collection_handle.isValid()) return result;

  /////////////////////////////////////////////////////////////////////////////
  //
  /////////////////////////////////////////////////////////////////////////////
  if (primary_vertex_collection_handle->empty()) {
    LogDebug("") << "no primary vertex";
    return result;
  }
  const reco::Vertex& primary_vertex = primary_vertex_collection_handle->at(0);

  const auto filter_module_set = findFilterModules(event, trigger_results_handle);
  if (filter_module_set.empty()) {
    LogDebug("") << "empty";
    return result;
  }

  const auto& trigger_object_momentum_vec = findTriggerObjectsMomenta(trigger_event_handle, filter_module_set);
  if (trigger_object_momentum_vec.empty()) {
    LogDebug("") << "no trigger object";
    return result;
  }

  result.tag = findTagMuon(muon_collection_handle, primary_vertex,
                                trigger_object_momentum_vec);
  if (result.tag.isNull()) {
    LogDebug("") << "no tag muon";
    return result;
  }

  result.probe = findProbeMuon(muon_collection_handle,
                                                 primary_vertex,
                                                 result.tag);
  if (result.probe.isNull()) {
    LogDebug("") << "no probe muon";
    return result;
  }

  for (const reco::MuonChamberMatch& muon_chamber_match : result.probe->matches()) {
    if (muon_chamber_match.detector() != MuonSubdetId::RPC) {
      continue;
    }

    if (checkIfMuonIsOutsideRoll(muon_chamber_match)) {
      LogDebug("") << "muon is outside the roll";
      continue;
    }

    const auto [hit, distance] = findClosestHit(muon_chamber_match, rpc_rec_hit_collection_handle);
    result.measurements.push_back(std::make_pair(muon_chamber_match, hit));
  }

  return result;
}

// Returns a rpctnp::MuonIdType enum corresponding to name in the case-insensitive way.
rpctnp::MuonIdType MuRPCTnPFlatTableProducer::convertStrToMuonIdType(std::string name) {
  std::transform(name.begin(), name.end(), name.begin(), ::tolower);

  rpctnp::MuonIdType id_type;
  if      (name == "tight")   id_type = rpctnp::MuonIdType::kTight;
  else if (name == "soft")    id_type = rpctnp::MuonIdType::kSoft;
  else if (name == "tracker") id_type = rpctnp::MuonIdType::kTracker;
  else if (name == "rpc")     id_type = rpctnp::MuonIdType::kRPC;
  else if (name == "none")    id_type = rpctnp::MuonIdType::kNone;
  else {
    edm::LogError("") << "got an unexpected muon id type: "
                      << name
                      << "The muon id is set to rpctnp::MuonIdType::kNone";
    id_type = rpctnp::MuonIdType::kNone;
  }

  return id_type;
}


// Returns the relative isolation with the beta correction with Particle Flow
// objects
double MuRPCTnPFlatTableProducer::computeMuonPFRelIso(const reco::MuonPFIsolation& iso,
                                                      const double muon_pt) {
  const double charged_pt_sum = iso.sumChargedHadronPt;
  const double neutral_pt_sum = iso.sumNeutralHadronEt + iso.sumPhotonEt - 0.5 * iso.sumPUPt;
  const double pt_sum = charged_pt_sum + std::max(0.0, neutral_pt_sum);
  return pt_sum / muon_pt;
}


// Returns a boolean indicating whether the tag muon is the triggering muon.
bool MuRPCTnPFlatTableProducer::doTriggerMatching(
    const reco::Muon& mu,
    const std::vector<math::XYZTLorentzVector>& trigger_object_momentum_vec) {

  bool found = false;
  for (const auto& trigger_object_p4 : trigger_object_momentum_vec) {
    const double delta_r = deltaR(mu, trigger_object_p4);
    const double delta_pt_rel = std::abs(mu.pt() - trigger_object_p4.pt()) / mu.pt();

    if ((delta_r < m_tag_muon_trigger_matching_max_delta_r) and (delta_pt_rel < m_tag_muon_trigger_matching_max_delta_pt_rel)) {
      found = true;
      break;
    }
  }
  return found;
}


// Returns a boolean indicating whether a muon passes the identification
// requirement.
bool MuRPCTnPFlatTableProducer::checkIfMuonPassId(const reco::Muon& mu,
                                                  const reco::Vertex& primary_vertex,
                                                  const rpctnp::MuonIdType& muon_id_type) {
  bool passed = false;

  switch (muon_id_type) {
    case rpctnp::MuonIdType::kTight: {
      passed = muon::isTightMuon(mu, primary_vertex);
      break;
    }
    case rpctnp::MuonIdType::kSoft: {
      passed = muon::isSoftMuon(mu, primary_vertex);
      break;
    }
    case rpctnp::MuonIdType::kTracker: {
      passed = mu.isTrackerMuon() and muon::isGoodMuon(mu, muon::TMOneStationLoose);
      break;
    }
    case rpctnp::MuonIdType::kRPC: {
      passed = mu.isRPCMuon() and muon::isGoodMuon(mu, muon::RPCMuLoose);
      break;
    }
    case rpctnp::MuonIdType::kNone: {
      passed = true;
      break;
    }
    default: {
      edm::LogError("") << "got an unexpected MuonIdType of "
                        << static_cast<int>(muon_id_type)
                        << "and falls through to default";
      break;
    }
  }

  return passed;
}

// Returns a boolean indicating whether a probe muon is in the fiducial area of
// RPCRoll.
bool MuRPCTnPFlatTableProducer::isFiducial(const reco::MuonChamberMatch& match) {
  const LocalPoint pos{match.x, match.y, 0.};

  const RPCDetId det_id{match.id};
  const RPCRoll* roll = m_rpc_geometry->roll(det_id);
  if (roll == nullptr) {
    edm::LogError("") << "skip this MuonChamberMatch because the RPCRoll is nullptr: " << det_id;
    return false;
  }
  const Bounds& bounds = roll->surface().bounds();

  const bool fiducial_y = (std::abs(pos.y()) <= bounds.length() / 2 - 8);

  bool fiducial_x = false;
  if (det_id.region() == 0) {
    fiducial_x = (std::abs(pos.x()) <= bounds.width() / 2 - 8);

  } else {
    const double wt = bounds.width();
    const double w0 = bounds.widthAtHalfLength();
    const double length = bounds.length();

    const double slope = (wt - w0) / length;
    const double w2_at_y = slope * pos.y() + w0 / 2;
    fiducial_x = (std::abs(pos.x()) <= w2_at_y - 8);
  }

  return fiducial_x and fiducial_y;
}


#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(MuRPCTnPFlatTableProducer);
