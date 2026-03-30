#ifndef RecoLocalMuon_RPCClusterizerPhase2_h
#define RecoLocalMuon_RPCClusterizerPhase2_h

#include "RecoLocalMuon/RPCRecHit/plugins/RPCRollMask.h"
#include "RecoLocalMuon/RPCRecHit/plugins/RPCClusterPhase2.h"
#include "RecoLocalMuon/RPCRecHit/plugins/RPCClusterPhase2Container.h"

#include "DataFormats/RPCDigi/interface/RPCDigiPhase2Collection.h"

class RPCClusterizerPhase2 {
public:
  RPCClusterizerPhase2() = default;
  ~RPCClusterizerPhase2() = default;

  RPCClusterPhase2Container doAction(const RPCDigiPhase2Collection::Range& digiRange, const RollMask& mask) const;

private:
  RPCClusterPhase2Container makeInitialClusters(const RPCDigiPhase2Collection::Range& digiRange) const;
  RPCClusterPhase2Container reclusterWithMask(const RPCClusterPhase2Container& input, const RollMask& mask) const;
  bool canMergeAcrossMaskedGap(const RPCClusterPhase2& left, const RPCClusterPhase2& right, const RollMask& mask) const;
};

#endif