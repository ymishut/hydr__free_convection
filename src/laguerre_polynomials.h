#ifndef SRC_LAGUERRE_POLYNOMIALS_H_
#define SRC_LAGUERRE_POLYNOMIALS_H_

#include <boost/noncopyable.hpp>

#include <cstddef>

#include <complex>
#include <exception>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace laguerre_polynoms {
//================================
// Laguerre_polynomial struct
//================================

struct Laguerre_polynomial {
  double a, c;
  int n, p, px;

  Laguerre_polynomial(double a, double c, int n, int p, int px);
  Laguerre_polynomial();
  void swap(Laguerre_polynomial &LP) noexcept;
  std::vector<double> getClosedForm() const;
};

std::ostream &operator<<(std::ostream &outstr,
                         const Laguerre_polynomial &lp);

bool operator==(const Laguerre_polynomial &lpl,
                const Laguerre_polynomial &lpr);

//================================
// Laguerre_polynomial_basis
//================================

class Laguerre_polynomial_basis : private boost::noncopyable {
 public:
  Laguerre_polynomial_basis();
  void addToBasis(double a, double c, int n, int p, int px);
  void addToBasis(Laguerre_polynomial &&LP);
  Laguerre_polynomial_basis *nextDerivate() const;
  // unpack basis to closed form (SUM coef[i] * x^i).
  double getValueCL(const double x) const;
  std::complex<double> getValue(const std::complex<double> coef,
                                const double x) const;
  void derivateLaguerrePolyns(const size_t derivateOrder);
  static double integrateLaguerrePolynsBasis(
                    Laguerre_polynomial_basis &LPBf,
                    Laguerre_polynomial_basis &LPBph,
                    const size_t derivOrder);
  static double integrateLaguerrePolynsBasis(
                    Laguerre_polynomial_basis &LPBf,
                    Laguerre_polynomial_basis &LPBph,
                    const size_t derivOrder, const double *U,
                    const double dh, const size_t JMax);

#ifdef UNIT_TEST
  std::multimap<int, Laguerre_polynomial> getPolynoms() const;
  std::vector<double> outputClosedForm() const;
#endif  // UNIT_TEST

 private:
  size_t depth_;
  mutable std::vector<double>                closedForm_;
  std::unique_ptr<Laguerre_polynomial_basis> nextDeriv_;
  std::multimap<int, Laguerre_polynomial>    polynoms_;

 private:
  explicit Laguerre_polynomial_basis(const size_t dep);
  bool addToBasisCheck(const double c);
  void derivateLaguerrePolyns(const size_t derivateOrder,
                              const size_t derivateOrderFinish);
  void getClosedForm() const;
  void derivateLaguerrePolynsSub(Laguerre_polynomial &LP);
  static double integrateLaguerrePair(Laguerre_polynomial &f,
                                      Laguerre_polynomial &ph);
  static double L_appr(const double x, const int n, const int a);
  static double integrateLaguerrePolynomials(
      const int i, const int k, const int a);
  static void equatingLaguerreXandPX(Laguerre_polynomial &forig,
                                     Laguerre_polynomial *f);

 private:
  // inteface_exceptions
  class LPBexception_empty;
  // calculation_exceptions
  class LPBexception_calculation;
};

//================================
// LPBexception_empty
//================================

class Laguerre_polynomial_basis::LPBexception_empty :
    public std::exception {
 public:
  LPBexception_empty() {}

  virtual const char *what() const noexcept {
    return " Attempt to use object Laguerre_polynomial_basis with empty"
           " basis. Call addToBasis(Laguerre_polynomial)";
  }

  ~LPBexception_empty() noexcept {}
};

//================================
// LPBexception_calculation
//================================

class Laguerre_polynomial_basis::LPBexception_calculation :
    public std::exception {
  std::string message_;

 public:
  LPBexception_calculation(std::string &&message) noexcept
    : message_(message) {}

  LPBexception_calculation(const char *message)
    : message_(message) {}

  virtual const char *what() const noexcept {
    return message_.c_str();
  }

  ~LPBexception_calculation() noexcept {}
};
}  // namespace laguerre_polynoms
#endif  // SRC_LAGUERRE_POLYNOMIALS_H_
