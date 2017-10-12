#ifndef SRC_BASICFLOW_H_
#define SRC_BASICFLOW_H_

#include "maindatastorage.h"

#include <boost/noncopyable.hpp>

#include <array>
#include <complex>
#include <memory>
#include <vector>

namespace conv_flow {
struct nachtsInput;
struct galerkInput;

//================================
// FlowParameters
//================================

class FlowParameters : private boost::noncopyable {
 protected:
  typedef std::complex<double> compd;
  typedef std::function<double (double, const double*,
              double, size_t)> adamFunc_bf;

 public:
  void useNachts(const nachtsInput &inp,
                 const std::array<compd, 4> &baseValues);
  void useGaler(const double Re, const double wavnum, galerkInput igp,
                std::vector<std::array<compd, 4> > &baseValues);
  const std::vector<double> &getFlowParameters(size_t i) const;
  const MainDataStorage *getMDS() const;
  virtual ~FlowParameters() {}

 protected:
  double *F, *F1, *F2, *F3, *H, *H1, *H2;
  std::unique_ptr<MainDataStorage> mdsp_;

 protected:
  FlowParameters(FlowID fid, const double Pr, size_t Jm,
                 const double zerosEdge, double dh);
  virtual void calculate() = 0;
  void checkResult();
};

//================================
// PolhausenFlow
//================================
// free convection boundary layer flow
class PolhausenFlow final : public FlowParameters {
 public:
  PolhausenFlow(const double Pr, size_t Jm, const double zerosEdge, double dh);

 private:
  void calculate() override;
  void RungeKutta();
  void AdamsIntegrate(adamFunc_bf f, size_t j);

 private:
  struct pushresult_observer;
  struct polhausSystem;
};

//================================
// SiegelFlow
//================================
// free convection with v = 0
class SiegelFlow final : public FlowParameters {
 public:
  SiegelFlow(const double Pr, size_t Jm, const double zerosEdge, double dh);

 private:
  void calculate() override;
};

//================================
// pushresult_observer struct
//================================

struct PolhausenFlow::pushresult_observer {
  const std::unique_ptr<conv_flow::MainDataStorage> &mds;
  explicit pushresult_observer(const std::unique_ptr<MainDataStorage> &mdsp_);

  void operator()(const std::array<double, 5> &x, double t);
};

//================================
// polhausSystem struct
//================================

struct PolhausenFlow::polhausSystem {
  const std::unique_ptr<conv_flow::MainDataStorage> &mds;
  explicit polhausSystem(
      const std::unique_ptr<conv_flow::MainDataStorage> &mdsp_);

  void operator()(const std::array<double, 5> &x,
                  std::array<double, 5> &dxdt, double);
};
}  // namespace conv_flow
#endif  // SRC_BASICFLOW_H_
