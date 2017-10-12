#ifndef SRC_NACHTS_METHOD_H_
#define SRC_NACHTS_METHOD_H_

#include <boost/noncopyable.hpp>

#include <array>
#include <complex>
#include <memory>
#include <vector>

namespace conv_flow {
class MainDataStorage;
class StabilityPoints;

//================================
// flowparameters struct
//================================
// fisi is struct of a "perturbation" parameters
// [0] - f
// [1] - f'
// ...etc
// [5] - s
// ... [7] s''
struct flowparameters : private boost::noncopyable {
  std::array<std::vector<std::complex<double>>, 8> fisi;
 
  explicit flowparameters(size_t n);
};

//================================
// nachtsInput struct
//================================

struct nachtsInput {
  const double Re, Re_max, dRe,
               wavnum_min, wavnum_max, dwavnum;

  nachtsInput(double Re, double REM, double dRe,
              double w_min, double w_max, double dw);
};

//================================
// nachtsMethod
//================================

class nachtsMethod {
  typedef std::complex<double> compd;

 public:
  typedef std::array<compd, 6> state_type_nach;
  typedef std::function<compd(const double, const compd *,
                              const compd, const size_t)> adamFunc_nach;

 public:
  // method calculate stability points of a flow and write them
  //  in stabilitypoints object ("sp")
  // one basicflow object can creat lot stabilitypoints objects
  // (i think it is flexible way
  // to work with many little Re-wavenumber fields).
  // input: Re_wavenumbounds={ReMax,ReMin,dRe,wavnumMax,wavnumMin,dwavenum},
  //      cbase,f3base,s1base - Galerkin method results
  explicit nachtsMethod(
      const MainDataStorage *mds, const bool uniqueResult = false);
  bool NachtsMethodCalculate(const nachtsInput &inp,
      const std::array<compd, 4> &baseValues, StabilityPoints &sp);
  const std::vector<compd> &getFlowparameters(size_t i) const;
  compd getPhasevelocity() const;

 private:
  const MainDataStorage *mds_;
  std::array<std::unique_ptr<flowparameters>, 4> parametersABC_;
  std::array<compd, 4> coef_;
  double wavnum_;
  compd  phvel_;
  const bool uniq_;

 private:
  struct pushresult_observer;
  struct nachtsheinSystem;

 private:
  void AdamsMethods(flowparameters &fp, const double dx,
                    size_t IMax, flowparameters *fc);
};

//================================
// pushresult_observer struct
//================================

struct nachtsMethod::pushresult_observer {
  size_t nOrder;
  nachtsMethod *nm;
  pushresult_observer(size_t nOrd, nachtsMethod *nm);
  void operator()(const nachtsMethod::state_type_nach &x, double t);
};

//================================
// nachtsheinSystem struct
//================================

struct nachtsMethod::nachtsheinSystem {
  bool phVelocDerivated;
  nachtsMethod *nm;
  nachtsheinSystem(const bool pvd, nachtsMethod *nm);
  void operator()(const nachtsMethod::state_type_nach &x,
                  nachtsMethod::state_type_nach &dxdt, double t);
};
}  // namespace conv_flow
#endif  // SRC_NACHTS_METHOD_H_
