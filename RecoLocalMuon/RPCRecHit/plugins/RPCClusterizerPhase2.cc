#include "RPCClusterizerPhase2.h"

RPCClusterPhase2Container RPCClusterizerPhase2::doAction(const RPCDigiPhase2Collection::Range& digiRange,
                                                         const RollMask& mask) const {
  RPCClusterPhase2Container initial = makeInitialClusters(digiRange);
  return reclusterWithMask(initial, mask);
}

RPCClusterPhase2Container RPCClusterizerPhase2::makeInitialClusters(const RPCDigiPhase2Collection::Range& digiRange) const {
  RPCClusterPhase2Container initialClusters;
  RPCClusterPhase2Container finalClusters;

  if (std::distance(digiRange.first, digiRange.second) == 0) {
    return finalClusters;
  }

  for (auto digi = digiRange.first; digi != digiRange.second; ++digi) {
    RPCClusterPhase2 cl(digi->strip(), digi->strip(), digi->bx());

    // FIXME: define how RPCDigiPhase2 sub-BX timing should be propagated to cluster time.
    // FIXME: once RPCDigiPhase2 API is frozen, convert sbx/fine timing into the cluster time observable here.
    // Example sketch:
    //   cl.addTime(convertPhase2Timing(*digi));

    // FIXME: existing Phase2 RPC digi may not provide local-y information.
    // Keep local-y empty for now.

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
  if (input.empty()) return output;

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
  if (left.bx() != right.bx()) return false;
  if (right.firstStrip() - left.lastStrip() != 2) return false;

  const int maskedStrip = left.lastStrip() + 1;

  // NOTE: RollMask is 0-based while strip numbering is 1-based in RPC.
  // This is the same convention used by the legacy producer when filling the mask bitset.
  return mask.test(maskedStrip - 1);
}