#ifndef SRC_CALCULATION_H_
#define SRC_CALCULATION_H_

#include "galerk_method.h"
#include "maindatastorage.h"
#include "nachts_method.h"

#include <array>
#include <complex>
#include <exception>
#include <memory>
#include <string>
#include <vector>

namespace conv_flow {
class FlowParameters;

//================================
// Calculation
//================================

class Calculation {
  typedef std::complex<double> compd;

 public:
  void calculating();
  void calculating(const double Re, const double wan,
                   const std::array<compd, 4> &stcond);
  void calculating(const double Re, const double wan);
  void setGalerkinInput(std::unique_ptr<galerkInput> newgin);
  void setNacgtsInput(std::unique_ptr<nachtsInput> newnin);

  Calculation(const Calculation &clc);
  Calculation &operator=(const Calculation &clc);
  Calculation(Calculation &&r)              = default;
  Calculation &operator=(Calculation &&r)   = default;

 public:
  class SelectOptions;

 private:
  std::shared_ptr<FlowParameters> flop_;
  std::unique_ptr<galerkInput>    gin_;
  std::unique_ptr<nachtsInput>    nin_;

 private:
  class CalculationExcept;

 private:
  Calculation(std::shared_ptr<FlowParameters> flop,
              std::unique_ptr<galerkInput> gin,
              std::unique_ptr<nachtsInput> nin);
  // few thread for various Re numbers
  void concurrNachts(const std::array<compd, 4> &stcond);
  // few thread for various stcond
  void concurrNachts(nachtsInput *nip, 
                     std::vector<std::array<compd, 4>> &stcond);
};

//================================
// SelectOptions - builder
//================================

class Calculation::SelectOptions {
 public:
  SelectOptions(FlowID fid, const double Pr, size_t Jm,
                const double zerosEdge, double dh);

  SelectOptions &useGalerkin(const size_t polynomPower,
                             const size_t matOrder,
                             const size_t startOrder);
  SelectOptions &useGalerkin(galerkInput&& gi);
  SelectOptions &useNachtsheim(
      const double Re_min, double Re_max, const double dRe,
      const double w_min, const double w_max, const double dw);
  Calculation build();

 private:
  std::shared_ptr<FlowParameters> flop;
  std::unique_ptr<galerkInput>    gin;
  std::unique_ptr<nachtsInput>    nin;
};

//================================
// CalculationExcept
//================================

class Calculation::CalculationExcept : public std::exception {
  std::string message_;

 public:
  CalculationExcept(std::string &&message) noexcept
    : message_(message) {}

  CalculationExcept(const char *message)
    : message_(message) {}

  virtual const char *what() const noexcept {
    return message_.c_str();
  }

  ~CalculationExcept() noexcept {}
};
}  // namespace conv_flow
#endif  // SRC_CALCULATION_H_

