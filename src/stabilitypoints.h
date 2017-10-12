#ifndef SRC_STABILITYPOINTS_H_
#define SRC_STABILITYPOINTS_H_

#include <cstddef>

#include <set>

namespace conv_flow {
class MainDataStorage;
struct stPoint;

//================================
// StabilityPoints
//================================

class StabilityPoints {
 public:
  StabilityPoints(MainDataStorage *msd);
  void addPoint(stPoint sp);
  void pushToStorage() const;

 private:
  MainDataStorage *mds_;
  std::set<stPoint> points_;
};
}  // namespace comv_flow
#endif  // SRC_STABILITYPOINTS_H_
