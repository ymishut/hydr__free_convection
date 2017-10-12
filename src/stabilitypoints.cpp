#include "stabilitypoints.h"

#include "maindatastorage.h"

#include <cmath>

namespace conv_flow {
StabilityPoints::StabilityPoints(MainDataStorage *msd)
  : mds_(msd) {}

void StabilityPoints::addPoint(stPoint sp) {
  points_.insert(sp);
}

void StabilityPoints::pushToStorage() const {
  std::lock_guard <std::mutex> (mds_->getMutex());
  mds_->points.insert(points_.begin(), points_.end());
}
}  // namespace conv_flow
