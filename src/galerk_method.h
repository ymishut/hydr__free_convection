#ifndef SRC_GALERK_METHOD_H_
#define SRC_GALERK_METHOD_H_

#include "laguerre_polynomials.h"

#include <armadillo>

#include <array>
#include <complex>
#include <string>
#include <vector>

namespace conv_flow {
class MainDataStorage;

//================================
// galerkInput
//================================

struct galerkInput {
  galerkInput(const size_t pp,
              const size_t mo,
              const size_t so);

  const size_t polynomsPower,
               countBasicFunc,
               startOrder;
};

//================================
// galerkMethod
//================================

class galerkMethod {
  typedef std::complex <double> compd;

 public:
  galerkMethod(galerkInput gi, MainDataStorage *mainDataStorage);
  void GalerkMethodCalculate(const double Re, const double wavnum,
      std::vector<std::array<compd, 4>> &baseValues);
  static bool setTriangular(arma::cx_mat &A, const size_t n);

#ifdef UNIT_TEST
  int chooseEigenvalue_test(std::vector<compd> &vec_ev) const;
  std::vector<size_t> setResult_test(std::vector<compd> &veceigvs) const;
#endif

 public:
  class galerkMethodExcept;

 private:
  const size_t polynomsPower_,
               countBasicFunc_;
        size_t startOrder_;
  MainDataStorage *mds_;

 private:
  void getBasicValues(arma::cx_mat &A, arma::cx_mat &B,
      std::array<compd, 4> &bv, const double expc);
  void fillMatrices(arma::cx_mat &A, arma::cx_mat &B,
      const std::array<double, 4>  &coef,
      const double wavnum, const double expc);
  // choose optimal from one calculation
  int chooseEigenvalue(std::vector<compd> &eigenValues);
  void setBasisFunc(laguerre_polynoms::Laguerre_polynomial_basis &fiLPB,
                    laguerre_polynoms::Laguerre_polynomial_basis &siLPB,
                    const double expc,
                    const size_t polOrder,
                    const size_t xPower);
  //  choose optimal from all previous
  //  set final strcond
  std::vector<size_t> setResult(std::vector<compd> &veceigvs);
};

//================================
// galerkMethodExcept
//================================

class galerkMethod::galerkMethodExcept : public std::exception {
  std::string message_;

 public:
  galerkMethodExcept(std::string &&message) noexcept
    : message_(message) {}

  galerkMethodExcept(const char *message)
    : message_(message) {}

  virtual const char *what() const noexcept {
    return message_.c_str();
  }

  ~galerkMethodExcept() noexcept {}
};
}  // namespace conv_flow
#endif  // SRC_GALERK_METHOD_H_
