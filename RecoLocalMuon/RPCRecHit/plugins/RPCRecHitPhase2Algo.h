#ifndef RecoLocalMuon_RPCRecHitPhase2Algo_h
#define RecoLocalMuon_RPCRecHitPhase2Algo_h

#include "Geometry/RPCGeometry/interface/RPCRoll.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometrySurface/interface/LocalError.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitPhase2.h"

#include "RPCClusterPhase2.h"

class RPCRecHitPhase2Algo {
public:
  RPCRecHitPhase2Algo() = default;
  ~RPCRecHitPhase2Algo() = default;

  bool compute(const RPCRoll& roll,
               const RPCClusterPhase2& cluster,
               LocalPoint& point,
               LocalError& error,
               float& time,
               float& timeErr) const;

  RPCRecHitPhase2 build(const RPCRoll& roll,
                        const RPCDetId& rpcId,
                        const RPCClusterPhase2& cluster) const;
};

#endif