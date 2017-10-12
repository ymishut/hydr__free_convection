#ifndef SRC_INPUTDATA_BY_FILE_H_
#define SRC_INPUTDATA_BY_FILE_H_

#include "maindatastorage.h"

#include <boost/noncopyable.hpp>

#include <array>
#include <string>
#include <utility>
#include <vector>

namespace conv_flow {
struct galerkInput;

//================================
// InputData
//================================

class InputData : private boost::noncopyable {
 public:
  explicit InputData(const std::string &filename);
  std::string               getCalcMethods()      const;
  std::pair<FlowID, double> getFlow()             const;
  std::array<double, 3>     getReBoundaries()     const;
  std::array<double, 3>     getWavnumBoundaries() const;
  std::array<double, 3>     getCalcParams();
  galerkInput               getGalerkParams();

 public:
  class InputExcept;

 private:
  std::array<double, 3> castArray(size_t i, std::string *str) const;
  std::array<double, 3> castBoundary(size_t i) const;

 private:
  std::vector<std::string> inputParameters_;
};

//================================
// InputDataExcept
//================================

class InputData::InputExcept : public std::exception{
  std::string message_;

 public:
  InputExcept(std::string &&message) noexcept
    : message_(message) {}

  InputExcept(const char *message)
    : message_(message) {}

  virtual const char *what() const noexcept {
    return message_.c_str();
  }

  ~InputExcept() noexcept {}
};
}  // namespace conv_flow
#endif  // SRC_INPUTDATA_BY_FILE_H_
