/**
  * 
  * See header file for a description of this class.
  *
  * \author Jongwon Shin - Kyung Hee University
  */

#include "DataFormats/RPCRecHit/interface/RPCRecHitPhase2.h"

RPCRecHitPhase2::RPCRecHitPhase2(const RPCDetId& rpcId, int bx)
    : RecHit2DLocalPos(rpcId),
      theRPCId(rpcId),
      theBx(bx),
      theFirstStrip(99),
      theClusterSize(99),
      theLocalPosition(),
      theLocalError(),
      theTime(0.f),
      theTimeError(-1.f) {}

RPCRecHitPhase2::RPCRecHitPhase2()
    : RecHit2DLocalPos(),
      theRPCId(),
      theBx(99),
      theFirstStrip(99),
      theClusterSize(99),
      theLocalPosition(),
      theLocalError(),
      theTime(0.f),
      theTimeError(-1.f) {}

RPCRecHitPhase2::RPCRecHitPhase2(const RPCDetId& rpcId, int bx, const LocalPoint& pos)
    : RecHit2DLocalPos(rpcId),
      theRPCId(rpcId),
      theBx(bx),
      theFirstStrip(99),
      theClusterSize(99),
      theLocalPosition(pos),
      theTime(0.f),
      theTimeError(-1.f) {
  // FIXME: Is this value correct?      
  float stripResolution = 3.0f;
  // FIXME: Is it really needed to set the local error here?
  theLocalError = LocalError(stripResolution * stripResolution, 0.f, 0.f);
}

RPCRecHitPhase2::RPCRecHitPhase2(const RPCDetId& rpcId, int bx, const LocalPoint& pos, const LocalError& err)
    : RecHit2DLocalPos(rpcId),
      theRPCId(rpcId),
      theBx(bx),
      theFirstStrip(99),
      theClusterSize(99),
      theLocalPosition(pos),
      theLocalError(err),
      theTime(0.f),
      theTimeError(-1.f) {}

RPCRecHitPhase2::RPCRecHitPhase2(
    const RPCDetId& rpcId, int bx, int firstStrip, int clustSize, const LocalPoint& pos, const LocalError& err)
    : RecHit2DLocalPos(rpcId),
      theRPCId(rpcId),
      theBx(bx),
      theFirstStrip(firstStrip),
      theClusterSize(clustSize),
      theLocalPosition(pos),
      theLocalError(err),
      theTime(0.f),
      theTimeError(-1.f) {}

RPCRecHitPhase2::~RPCRecHitPhase2() {}

RPCRecHitPhase2* RPCRecHitPhase2::clone() const { return new RPCRecHitPhase2(*this); }

std::vector<const TrackingRecHit*> RPCRecHitPhase2::recHits() const {
  std::vector<const TrackingRecHit*> nullvector;
  return nullvector;
}

std::vector<TrackingRecHit*> RPCRecHitPhase2::recHits() {
  std::vector<TrackingRecHit*> nullvector;
  return nullvector;
}

bool RPCRecHitPhase2::operator==(const RPCRecHitPhase2& hit) const {
  return this->geographicalId() == hit.geographicalId();
}

std::ostream& operator<<(std::ostream& os, const RPCRecHitPhase2& hit) {
  os << "pos: " << hit.localPosition().x();
  os << " +/- " << std::sqrt(hit.localPositionError().xx());
  return os;
}