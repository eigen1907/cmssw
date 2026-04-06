#ifndef RecoLocalMuon_RPCRecHit_RPCRecHitPhase2Producer_h
#define RecoLocalMuon_RPCRecHit_RPCRecHitPhase2Producer_h

/*
 *  Producer for RPCRecHitPhase2 from RPCDigiPhase2 and IRPCDigi.
 *
 *  RPCDigiPhase2 -> RPCClusterPhase2 -> RPCRecHitPhase2
 *  IRPCDigi -> IRPCCluster -> RPCRecHitPhase2
 *
 *  \author J. Shin -- Kyung Hee University
 */

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/RPCDigi/interface/RPCDigiPhase2Collection.h"
#include "DataFormats/RPCDigi/interface/IRPCDigiCollection.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitPhase2Collection.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"

#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "RPCClusterizerPhase2.h"
#include "IRPCClusterizer.h"
#include "RPCRecHitPhase2Algo.h"

class RPCRoll;

class RPCRecHitPhase2Producer : public edm::stream::EDProducer<> {
public:
  explicit RPCRecHitPhase2Producer(const edm::ParameterSet& config);
  ~RPCRecHitPhase2Producer() override = default;

  void produce(edm::Event& event, const edm::EventSetup& setup) override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void buildFromRPCPhase2(const RPCRoll& roll,
                          const RPCDetId& rpcId,
                          const RPCDigiPhase2Collection::Range& range,
                          RPCRecHitPhase2Collection& output) const;

  void buildFromIRPC(const RPCRoll& roll,
                     const RPCDetId& rpcId,
                     const IRPCDigiCollection::Range& range,
                     RPCRecHitPhase2Collection& output) const;

  static bool isIRPCDetId(const RPCDetId& rpcId);

private:
  edm::EDGetTokenT<RPCDigiPhase2Collection> rpcDigiPhase2Token_;
  edm::EDGetTokenT<IRPCDigiCollection> irpcDigiToken_;
  edm::ESGetToken<RPCGeometry, MuonGeometryRecord> rpcGeomToken_;

  bool useIRPC_;

  RPCClusterizerPhase2 rpcClusterizer_;
  IRPCClusterizer irpcClusterizer_;
  RPCRecHitPhase2Algo algo_;
};

#endif