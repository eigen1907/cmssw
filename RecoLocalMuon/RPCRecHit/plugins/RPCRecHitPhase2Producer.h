#ifndef RecoLocalMuon_RPCRecHitPhase2Producer_h
#define RecoLocalMuon_RPCRecHitPhase2Producer_h

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/RPCDigi/interface/RPCDigiPhase2Collection.h"
#include "DataFormats/RPCDigi/interface/IRPCDigiCollection.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitPhase2Collection.h"

#include "CondFormats/RPCObjects/interface/RPCMaskedStrips.h"
#include "CondFormats/RPCObjects/interface/RPCDeadStrips.h"
#include "CondFormats/DataRecord/interface/RPCMaskedStripsRcd.h"
#include "CondFormats/DataRecord/interface/RPCDeadStripsRcd.h"

#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "RPCClusterizerPhase2.h"
#include "RPCRecHitPhase2Algo.h"
#include "RPCRollMask.h"

class RPCRecHitPhase2Producer : public edm::stream::EDProducer<> {
public:
  explicit RPCRecHitPhase2Producer(const edm::ParameterSet& config);
  ~RPCRecHitPhase2Producer() override = default;

  void beginRun(const edm::Run&, const edm::EventSetup&) override;
  void produce(edm::Event& event, const edm::EventSetup& setup) override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  enum class MaskSource { File, EventSetup };

  RollMask buildMaskForRoll(const RPCDetId& rpcId) const;

  void buildFromRPCPhase2(const RPCRoll& roll,
                          const RPCDetId& rpcId,
                          const RPCDigiPhase2Collection::Range& range,
                          RPCRecHitPhase2Collection& output) const;

  void buildFromIRPC(const RPCRoll& roll,
                     const RPCDetId& rpcId,
                     const IRPCDigiCollection::Range& range,
                     RPCRecHitPhase2Collection& output) const;

private:
  edm::EDGetTokenT<RPCDigiPhase2Collection> rpcDigiPhase2Token_;
  edm::EDGetTokenT<IRPCDigiCollection> irpcDigiToken_;

  edm::ESGetToken<RPCMaskedStrips, RPCMaskedStripsRcd> maskedStripsToken_;
  edm::ESGetToken<RPCDeadStrips, RPCDeadStripsRcd> deadStripsToken_;
  edm::ESGetToken<RPCGeometry, MuonGeometryRecord> rpcGeomToken_;

  bool useIRPC_;
  MaskSource maskSource_;
  MaskSource deadSource_;

  std::unique_ptr<RPCMaskedStrips> rpcMaskedStripsObj_;
  std::unique_ptr<RPCDeadStrips> rpcDeadStripsObj_;

  RPCClusterizerPhase2 clusterizer_;
  RPCRecHitPhase2Algo algo_;
};

#endif