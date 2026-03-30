#ifndef RecoLocalMuon_IRPCClusterizer_h
#define RecoLocalMuon_IRPCClusterizer_h

#include "RecoLocalMuon/RPCRecHit/plugins/RPCRollMask.h"
#include "RecoLocalMuon/RPCRecHit/plugins/IRPCCluster.h"
#include "RecoLocalMuon/RPCRecHit/plugins/IRPCClusterContainer.h"

#include "DataFormats/RPCDigi/interface/IRPCDigiCollection.h"

class IRPCClusterizer {
public:
  IRPCClusterizer() = default;
  ~IRPCClusterizer() = default;

  IRPCClusterContainer doAction(const IRPCDigiCollection::Range& digiRange, const RollMask& mask) const;

private:
  IRPCClusterContainer makeInitialClusters(const IRPCDigiCollection::Range& digiRange) const;
  IRPCClusterContainer reclusterWithMask(const IRPCClusterContainer& input, const RollMask& mask) const;
  bool canMergeAcrossMaskedGap(const IRPCCluster& left, const IRPCCluster& right, const RollMask& mask) const;
};

#endif