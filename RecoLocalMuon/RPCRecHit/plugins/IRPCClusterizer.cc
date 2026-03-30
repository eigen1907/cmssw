#include "RecoLocalMuon/RPCRecHit/plugins/IRPCClusterizer.h"

#include "DataFormats/RPCDigi/interface/IRPCDigiTime.h"

#include <iterator>

IRPCClusterContainer IRPCClusterizer::doAction(const IRPCDigiCollection::Range& digiRange, const RollMask& mask) const {
  const IRPCClusterContainer initial = makeInitialClusters(digiRange);
  return reclusterWithMask(initial, mask);
}

IRPCClusterContainer IRPCClusterizer::makeInitialClusters(const IRPCDigiCollection::Range& digiRange) const {
  IRPCClusterContainer initialClusters;
  IRPCClusterContainer finalClusters;

  if (std::distance(digiRange.first, digiRange.second) == 0) {
    return finalClusters;
  }

  for (auto digi = digiRange.first; digi != digiRange.second; ++digi) {
    IRPCCluster cl(digi->strip(), digi->strip(), digi->bx());

    // Dummy-but-wired implementation:
    // use the helper so bx/sbx/fine inputs are already consumed through the
    // IRPC timing model, but keep the final logic intentionally simple.
    IRPCDigiTime irpcTime(*digi);
    const float t = irpcTime.time();
    const float tLR = irpcTime.timeLR();
    const float tHR = irpcTime.timeHR();
    const float y = irpcTime.coordinateY();

    cl.addTime(t);
    cl.addY(y);

    // Keep these touched so the full input contract is already exercised.
    (void)tLR;
    (void)tHR;

    initialClusters.insert(cl);
  }

  if (initialClusters.empty()) {
    return finalClusters;
  }

  IRPCCluster current = *initialClusters.begin();
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

IRPCClusterContainer IRPCClusterizer::reclusterWithMask(const IRPCClusterContainer& input, const RollMask& mask) const {
  IRPCClusterContainer output;
  if (input.empty()) {
    return output;
  }

  auto left = input.begin();
  auto right = std::next(left);
  IRPCCluster current = *left;

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

bool IRPCClusterizer::canMergeAcrossMaskedGap(const IRPCCluster& left,
                                              const IRPCCluster& right,
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