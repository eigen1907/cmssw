#include "RecoLocalMuon/RPCRecHit/plugins/RPCRecHitPhase2Producer.h"

#include "DataFormats/MuonDetId/interface/RPCDetId.h"

#include "Geometry/RPCGeometry/interface/RPCRoll.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

RPCRecHitPhase2Producer::RPCRecHitPhase2Producer(const edm::ParameterSet& config)
    : rpcDigiPhase2Token_(consumes<RPCDigiPhase2Collection>(config.getParameter<edm::InputTag>("rpcDigiPhase2Label"))),
      irpcDigiToken_(consumes<IRPCDigiCollection>(config.getParameter<edm::InputTag>("irpcDigiLabel"))),
      maskedStripsToken_(esConsumes<RPCMaskedStrips, RPCMaskedStripsRcd, edm::Transition::BeginRun>()),
      deadStripsToken_(esConsumes<RPCDeadStrips, RPCDeadStripsRcd, edm::Transition::BeginRun>()),
      rpcGeomToken_(esConsumes<RPCGeometry, MuonGeometryRecord>()),
      useIRPC_(config.getParameter<bool>("useIRPC")),
      rpcMaskedStripsObj_(std::make_unique<RPCMaskedStrips>()),
      rpcDeadStripsObj_(std::make_unique<RPCDeadStrips>()) {
  produces<RPCRecHitPhase2Collection>();
}

void RPCRecHitPhase2Producer::beginRun(const edm::Run&, const edm::EventSetup& setup) {
  rpcMaskedStripsObj_->MaskVec = setup.getData(maskedStripsToken_).MaskVec;
  rpcDeadStripsObj_->DeadVec = setup.getData(deadStripsToken_).DeadVec;
}

RollMask RPCRecHitPhase2Producer::buildMaskForRoll(const RPCDetId& rpcId) const {
  RollMask mask;
  const int rawId = rpcId.rawId();

  for (const auto& item : rpcMaskedStripsObj_->MaskVec) {
    if (item.rawId == rawId) {
      mask.set(item.strip - 1);
    }
  }
  for (const auto& item : rpcDeadStripsObj_->DeadVec) {
    if (item.rawId == rawId) {
      mask.set(item.strip - 1);
    }
  }

  return mask;
}

void RPCRecHitPhase2Producer::buildFromRPCPhase2(const RPCRoll& roll,
                                                 const RPCDetId& rpcId,
                                                 const RPCDigiPhase2Collection::Range& range,
                                                 RPCRecHitPhase2Collection& output) const {
  const RollMask mask = buildMaskForRoll(rpcId);
  const RPCClusterPhase2Container clusters = rpcClusterizer_.doAction(range, mask);

  edm::OwnVector<RPCRecHitPhase2> hits;
  for (const auto& cluster : clusters) {
    hits.push_back(algo_.build(roll, rpcId, cluster));
  }

  if (!hits.empty()) {
    output.put(rpcId, hits.begin(), hits.end());
  }
}

void RPCRecHitPhase2Producer::buildFromIRPC(const RPCRoll& roll,
                                            const RPCDetId& rpcId,
                                            const IRPCDigiCollection::Range& range,
                                            RPCRecHitPhase2Collection& output) const {
  const RollMask mask = buildMaskForRoll(rpcId);
  const IRPCClusterContainer clusters = irpcClusterizer_.doAction(range, mask);

  edm::OwnVector<RPCRecHitPhase2> hits;
  for (const auto& cluster : clusters) {
    hits.push_back(algo_.build(roll, rpcId, cluster));
  }

  if (!hits.empty()) {
    output.put(rpcId, hits.begin(), hits.end());
  }
}

void RPCRecHitPhase2Producer::produce(edm::Event& event, const edm::EventSetup& setup) {
  const auto& rpcGeom = setup.getData(rpcGeomToken_);

  edm::Handle<RPCDigiPhase2Collection> rpcPhase2Digis;
  event.getByToken(rpcDigiPhase2Token_, rpcPhase2Digis);

  edm::Handle<IRPCDigiCollection> irpcDigis;
  if (useIRPC_) {
    event.getByToken(irpcDigiToken_, irpcDigis);
  }

  auto recHitCollection = std::make_unique<RPCRecHitPhase2Collection>();

  if (rpcPhase2Digis.isValid()) {
    for (auto rpcdgIt = rpcPhase2Digis->begin(); rpcdgIt != rpcPhase2Digis->end(); ++rpcdgIt) {
      const RPCDetId& rpcId = (*rpcdgIt).first;
      const RPCRoll* roll = rpcGeom.roll(rpcId);
      if (roll == nullptr) {
        edm::LogError("BadDigiInput") << "Failed to find RPCRoll for ID " << rpcId;
        continue;
      }

      const RPCDigiPhase2Collection::Range& range = (*rpcdgIt).second;
      buildFromRPCPhase2(*roll, rpcId, range, *recHitCollection);
    }
  } else {
    edm::LogWarning("MissingInput") << "RPCRecHitPhase2Producer: missing RPCDigiPhase2Collection";
  }

  if (useIRPC_) {
    if (irpcDigis.isValid()) {
      for (auto irpcIt = irpcDigis->begin(); irpcIt != irpcDigis->end(); ++irpcIt) {
        const RPCDetId& rpcId = (*irpcIt).first;
        const RPCRoll* roll = rpcGeom.roll(rpcId);
        if (roll == nullptr) {
          edm::LogError("BadDigiInput") << "Failed to find RPCRoll for ID " << rpcId;
          continue;
        }

        const IRPCDigiCollection::Range& range = (*irpcIt).second;
        buildFromIRPC(*roll, rpcId, range, *recHitCollection);
      }
    } else {
      edm::LogWarning("MissingInput") << "RPCRecHitPhase2Producer: missing IRPCDigiCollection";
    }
  }

  event.put(std::move(recHitCollection));
}

void RPCRecHitPhase2Producer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("rpcDigiPhase2Label", edm::InputTag("simMuonRPCDigisPhase2"));
  desc.add<edm::InputTag>("irpcDigiLabel", edm::InputTag("simMuonIRPCDigis"));
  desc.add<bool>("useIRPC", true);
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(RPCRecHitPhase2Producer);