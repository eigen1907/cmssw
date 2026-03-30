#include "RecoLocalMuon/RPCRecHit/plugins/IRPCCluster.h"

#include <algorithm>
#include <cmath>

IRPCCluster::IRPCCluster()
    : fstrip_(0), lstrip_(0), bx_(0), sumTime_(0.f), sumTime2_(0.f), nTime_(0), sumY_(0.f), sumY2_(0.f), nY_(0) {}

IRPCCluster::IRPCCluster(int firstStrip, int lastStrip, int bx)
    : fstrip_(firstStrip),
      lstrip_(lastStrip),
      bx_(bx),
      sumTime_(0.f),
      sumTime2_(0.f),
      nTime_(0),
      sumY_(0.f),
      sumY2_(0.f),
      nY_(0) {}

IRPCCluster::~IRPCCluster() = default;

int IRPCCluster::firstStrip() const { return fstrip_; }
int IRPCCluster::lastStrip() const { return lstrip_; }
int IRPCCluster::clusterSize() const { return lstrip_ - fstrip_ + 1; }
int IRPCCluster::bx() const { return bx_; }

bool IRPCCluster::hasTime() const { return nTime_ > 0; }
float IRPCCluster::time() const { return hasTime() ? sumTime_ / nTime_ : 0.f; }

float IRPCCluster::timeRMS() const {
  if (!hasTime()) {
    return -1.f;
  }
  const float variance = std::max(0.f, sumTime2_ / nTime_ - time() * time());
  return std::sqrt(variance);
}

bool IRPCCluster::hasY() const { return nY_ > 0; }
float IRPCCluster::y() const { return hasY() ? sumY_ / nY_ : 0.f; }

float IRPCCluster::yRMS() const {
  if (!hasY()) {
    return -1.f;
  }
  const float variance = std::max(0.f, sumY2_ / nY_ - y() * y());
  return std::sqrt(variance);
}

void IRPCCluster::addTime(float time) {
  ++nTime_;
  sumTime_ += time;
  sumTime2_ += time * time;
}

void IRPCCluster::addY(float y) {
  ++nY_;
  sumY_ += y;
  sumY2_ += y * y;
}

void IRPCCluster::merge(const IRPCCluster& other) {
  const int first = std::min(firstStrip(), other.firstStrip());
  const int last  = std::max(lastStrip(), other.lastStrip());

  fstrip_ = static_cast<uint16_t>(first);
  lstrip_ = static_cast<uint16_t>(last);

  nTime_ += other.nTime_;
  sumTime_ += other.sumTime_;
  sumTime2_ += other.sumTime2_;

  nY_ += other.nY_;
  sumY_ += other.sumY_;
  sumY2_ += other.sumY2_;
}

bool IRPCCluster::operator<(const IRPCCluster& other) const {
  if (bx_ != other.bx()) {
    return bx_ < other.bx();
  }
  return firstStrip() < other.firstStrip();
}

bool IRPCCluster::operator==(const IRPCCluster& other) const {
  return bx_ == other.bx() && firstStrip() == other.firstStrip() && clusterSize() == other.clusterSize();
}

bool IRPCCluster::isAdjacent(const IRPCCluster& other) const {
  return (other.lastStrip() + 1 == this->firstStrip()) && (other.bx() == this->bx());
}