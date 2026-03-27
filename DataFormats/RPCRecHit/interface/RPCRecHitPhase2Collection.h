#ifndef DataFormats_RPCRecHitPhase2Collection_H
#define DataFormats_RPCRecHitPhase2Collection_H

/** \class RPCRecHitCollection
  * 
  * Collection of RPCRecHit for storage in the event
  *
  * \author Jongwon Shin - Kyung Hee University
  */

#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitPhase2.h"
#include "DataFormats/Common/interface/RangeMap.h"
#include "DataFormats/Common/interface/ClonePolicy.h"
#include "DataFormats/Common/interface/OwnVector.h"
#include <functional>

typedef edm::RangeMap<RPCDetId,
                      edm::OwnVector<RPCRecHitPhase2, edm::ClonePolicy<RPCRecHitPhase2>>,
                      edm::ClonePolicy<RPCRecHitPhase2>>
    RPCRecHitPhase2Collection;

#endif