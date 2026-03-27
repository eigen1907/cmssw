#include "RecoLocalMuon/RPCRecHit/plugins/RPCRecHitPhase2Algo.h"
#include "Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h"

bool RPCRecHitPhase2Algo::compute(const RPCRoll& roll,
                                  const RPCClusterPhase2& cluster,
                                  LocalPoint& point,
                                  LocalError& error,
                                  float& time,
                                  float& timeErr) const {
  const float fstrip = roll.centreOfStrip(cluster.firstStrip()).x();
  const float lstrip = roll.centreOfStrip(cluster.lastStrip()).x();
  const float centreOfCluster = 0.5f * (fstrip + lstrip);

  const double y = cluster.hasY() ? cluster.y() : 0.0;
  point = LocalPoint(centreOfCluster, y, 0.0);

  if (!cluster.hasY()) {
    error = LocalError(roll.localError((cluster.firstStrip() + cluster.lastStrip()) / 2.0));
  } else {
    float ex2 = roll.localError((cluster.firstStrip() + cluster.lastStrip()) / 2.0).xx();

    const float stripLen = roll.specificTopology().stripLength();
    const float maxDy = stripLen / 2.f - std::abs(cluster.y());

    if (roll.id().region() != 0) {
      const auto& topo = dynamic_cast<const TrapezoidalStripTopology&>(roll.topology());
      const double angle = topo.stripAngle((cluster.firstStrip() + cluster.lastStrip()) / 2.0);
      const double x = centreOfCluster - y * std::tan(angle);

      point = LocalPoint(x, y, 0.0);

      const double scale = topo.localPitch(point) / topo.pitch();
      ex2 *= scale * scale;
    }

    error = LocalError(ex2, 0.f, maxDy * maxDy / 3.f);
  }

  if (cluster.hasTime()) {
    time = cluster.time();
    timeErr = cluster.timeRMS();
  } else {
    // FIXME: once Phase2 timing convention is fixed, revisit whether "0, -1" is still the best default.
    time = 0.f;
    timeErr = -1.f;
  }

  return true;
}

RPCRecHitPhase2 RPCRecHitPhase2Algo::build(const RPCRoll& roll,
                                           const RPCDetId& rpcId,
                                           const RPCClusterPhase2& cluster) const {
  LocalPoint point;
  LocalError error;
  float time = 0.f;
  float timeErr = -1.f;

  const bool ok = compute(roll, cluster, point, error, time, timeErr);

  // FIXME: decide whether to throw, skip, or return a default-constructed hit on compute failure.
  (void)ok;

  RPCRecHitPhase2 hit(rpcId,
                      cluster.bx(),
                      cluster.firstStrip(),
                      cluster.clusterSize(),
                      point,
                      error);

  hit.setTimeAndError(time, timeErr);
  return hit;
}