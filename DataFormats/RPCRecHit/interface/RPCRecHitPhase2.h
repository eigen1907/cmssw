#ifndef DataFormats_RPCRecHitPhase2_H
#define DataFormats_RPCRecHitPhase2_H

/** \class RPCRecHitPhase2
  *
  * RecHit for Resistive Plate Chamber, after Phase2 upgrade
  *
  * \author Jongwon Shin - Kyung Hee University
*/

#include "DataFormats/TrackingRecHit/interface/RecHit2DLocalPos.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"

class RPCRecHitPhase2 : public RecHit2DLocalPos {
public:
  RPCRecHitPhase2(const RPCDetId& rpcId, int bx);
  
  RPCRecHitPhase2();
  
  RPCRecHitPhase2(const RPCDetId& rpcId, int bx, const LocalPoint& pos);
  
  RPCRecHitPhase2(const RPCDetId& rpcId, int bx, const LocalPoint& pos, const LocalError& err);
  
  RPCRecHitPhase2(const RPCDetId& rpcId, int bx, int firstStrip, int clustSize, const LocalPoint& pos, const LocalError& err);
  
  ~RPCRecHitPhase2() override;

  LocalPoint localPosition() const override { return theLocalPosition; }
  
  LocalError localPositionError() const override { return theLocalError; }

  RPCRecHitPhase2* clone() const override;

  std::vector<const TrackingRecHit*> recHits() const override;
  
  std::vector<TrackingRecHit*> recHits() override;

  void setPosition(LocalPoint pos) { theLocalPosition = pos; }
  
  void setError(LocalError err) { theLocalError = err; }
  
  void setPositionAndError(LocalPoint pos, LocalError err) {
    theLocalPosition = pos;
    theLocalError = err;
  }

  void setTimeAndError(float time, float err) {
    theTime = time;
    theTimeError = err;
  }

  RPCDetId rpcId() const { return theRPCId; }
  
  int BunchX() const { return theBx; }
  
  int firstClusterStrip() const { return theFirstStrip; }
  
  int clusterSize() const { return theClusterSize; }

  float time() const { return theTime; }
  
  float timeError() const { return theTimeError; }

  bool operator==(const RPCRecHitPhase2& hit) const;

private:
  RPCDetId theRPCId;
  int theBx;
  int theFirstStrip;
  int theClusterSize;

  LocalPoint theLocalPosition;
  LocalError theLocalError;

  float theTime, theTimeError;
};
#endif

std::ostream& operator<<(std::ostream& os, const RPCRecHitPhase2& hit);
