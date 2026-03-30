#include "RecoLocalMuon/RPCRecHit/plugins/RPCClusterizerPhase2.h"

#include "DataFormats/RPCDigi/interface/RPCDigiPhase2Time.h"

#include <iterator>

RPCClusterPhase2Container RPCClusterizerPhase2::doAction(const RPCDigiPhase2Collection::Range& digiRange,
                                                         const RollMask& mask) const {
  const RPCClusterPhase2Container initial = makeInitialClusters(digiRange);
  return reclusterWithMask(initial, mask);
}

RPCClusterPhase2Container RPCClusterizerPhase2::makeInitialClusters(
    const RPCDigiPhase2Collection::Range& digiRange) const {
  RPCClusterPhase2Container initialClusters;
  RPCClusterPhase2Container finalClusters;

  if (std::distance(digiRange.first, digiRange.second) == 0) {
    return finalClusters;
  }

  for (auto digi = digiRange.first; digi != digiRange.second; ++digi) {
    RPCClusterPhase2 cl(digi->strip(), digi->strip(), digi->bx());

    // Skeleton-only hook:
    // touch the Phase2 timing helper so the input contract is already wired,
    // but intentionally do not use it yet in the actual recHit logic.
    RPCDigiPhase2Time phase2Time(*digi);
    const float dummyTime = phase2Time.time();
    (void)dummyTime;

    initialClusters.insert(cl);
  }

  if (initialClusters.empty()) {
    return finalClusters;
  }

  RPCClusterPhase2 current = *initialClusters.begin();
  for (auto it = std::next(initialClusters.begin()); it != initialClusters.end(); ++it) {
    if (current.isAdjacent(*it)) {
      current.merge(*it);
    } else {
      finalClusters.insert(current);
      current = *it;
    }
  }
  finalClusters.insert(current);

  return finalClusters;
}

RPCClusterPhase2Container RPCClusterizerPhase2::reclusterWithMask(const RPCClusterPhase2Container& input,
                                                                  const RollMask& mask) const {
  RPCClusterPhase2Container output;
  if (input.empty()) {
    return output;
  }

  auto left = input.begin();
  auto right = std::next(left);
  RPCClusterPhase2 current = *left;

  while (right != input.end()) {
    if (canMergeAcrossMaskedGap(current, *right, mask)) {
      current.merge(*right);
    } else {
      output.insert(current);
      current = *right;
    }
    ++right;
  }

  output.insert(current);
  return output;
}

bool RPCClusterizerPhase2::canMergeAcrossMaskedGap(const RPCClusterPhase2& left,
                                                   const RPCClusterPhase2& right,
                                                   const RollMask& mask) const {
  if (left.bx() != right.bx()) {
    return false;
  }
  if (right.firstStrip() - left.lastStrip() != 2) {
    return false;
  }

  const int maskedStrip = left.lastStrip() + 1;
  if (maskedStrip <= 0 || maskedStrip > static_cast<int>(mask.size())) {
    return false;
  }

  return mask.test(maskedStrip - 1);
}