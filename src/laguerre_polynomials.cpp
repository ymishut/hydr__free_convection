#include "laguerre_polynomials.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>

extern double ACCURACY;
std::vector<double> factresults {1.0, 1.0};

double factorial(const size_t i) {
  if (i < factresults.size())
    return factresults[i];
  for (size_t j = factresults.size() - 1; j < i ; ++j) {
      factresults.push_back(factresults[j] * (1+j));
    }
  return factresults[i];
}

double binomialcoef(size_t n, size_t k) {
  if (k > n)
    return 1;
  return factorial(n) / factorial(n - k) / factorial(k);
}

namespace laguerre_polynoms {
//================================
// Laguerre_polynomial::getClosedForm
//================================

std::vector<double>
Laguerre_polynomial::getClosedForm() const {
  if (n < 0)
    return std::vector<double>();
  std::vector<double> result;
  for (int i = 0; i < px; ++i)
    result.push_back(0.0);
  for (int i = 0; i <= n; ++i)
    result.push_back(a*std::pow(-1, i)*binomialcoef(n+p, n-i)/
        factorial(i));
  return result;
}

//================================
// Laguerre_polynomial_basis::getClosedForm
//================================

void Laguerre_polynomial_basis::getClosedForm() const {
  if (!closedForm_.empty())
    return;
  std::vector <double> result;
  for_each(polynoms_.begin(), polynoms_.end(),
      [&result](const std::pair<const int,
                Laguerre_polynomial> &pr) {
          std::vector<double> vc = pr.second.getClosedForm();
          size_t end_it = std::min(result.size(), vc.size());
          std::transform(result.begin(), result.begin() + end_it,
              vc.begin(), result.begin(), std::plus<double>());
          if (vc.size() > result.size())
            result.insert(result.end(),
                vc.begin()+result.size(), vc.end());
        });
  closedForm_.swap(result);
}

//================================
// Laguerre_polynomial_basis::getValueCL
//================================

double Laguerre_polynomial_basis::getValueCL(
    const double x) const {
  if (polynoms_.empty())
    throw LPBexception_empty();
  if (closedForm_.empty())
    this->getClosedForm();
  if (closedForm_.empty())
    throw LPBexception_empty();
  if (x == 0.0)
    return closedForm_[0];
  double result = closedForm_[0];
  int i = 1;
  for_each(closedForm_.begin()+1, closedForm_.end(),
      [&result, &i, x](double coef) {
          result += coef*std::pow(x, i++);
        });
  return result*std::exp(polynoms_.begin()->second.c * x);
}

//================================
// Laguerre_polynomial_basis::getValue
//================================

std::complex<double> Laguerre_polynomial_basis::getValue(
    const std::complex<double> coef, const double x) const {
  std::complex<double> result(0.0, 0.0);
  for_each(polynoms_.begin(), polynoms_.end(),
      [&result, &coef, x](const std::pair<const int,
                          Laguerre_polynomial> &pr) {
          result += (coef*pr.second.a *L_appr(x, pr.second.n, pr.second.p) *
              std::exp(pr.second.c*x)*std::pow(x, pr.second.px));
        });
  return result;
}

//================================
// Laguerre_polynomial_basis ctors
//================================

Laguerre_polynomial_basis::Laguerre_polynomial_basis(const size_t dep)
  : depth_(dep), nextDeriv_(nullptr) {}

Laguerre_polynomial_basis::Laguerre_polynomial_basis()
  : depth_(0), nextDeriv_(nullptr) {}

//================================
// Laguerre_polynomial_basis::addToBasisCheck()
//================================

bool Laguerre_polynomial_basis::addToBasisCheck(const double c) {
  double lpb_c = polynoms_.begin()->second.c;
  return (((c + ACCURACY) > lpb_c) == ((lpb_c + ACCURACY) > c));
}

//================================
// Laguerre_polynomial_basis::addToBasis(...)
//================================

void Laguerre_polynomial_basis::addToBasis(
    double a, double c, int n, int p, int px) {
  if (!polynoms_.empty()) {
      if (!addToBasisCheck(c))
        throw LPBexception_calculation(
            " All basis member must have same exp arg ");

        // it's maybe not mistake, but it's better not to do this
        assert(nextDeriv_ == nullptr);
    }
  polynoms_.emplace(n+p+px, Laguerre_polynomial(a, c, n, p, px));
}

void Laguerre_polynomial_basis::addToBasis(Laguerre_polynomial &&LP) {
  if (!polynoms_.empty()) {
      if (!addToBasisCheck(LP.c))
        throw LPBexception_calculation(
            " All basis member must have same exp arg ");
        // it's maybe not mistake, but it's better not to do this
        assert(nextDeriv_ == nullptr);
    }
  polynoms_.emplace(LP.n + LP.p + LP.px, LP);
}

//================================
// Laguerre_polynomial_basis::nextDerivate
//================================

Laguerre_polynomial_basis *
Laguerre_polynomial_basis::nextDerivate() const {
  return (nextDeriv_ == nullptr) ? nullptr : nextDeriv_.get();
}

//================================
// Laguerre_polynomial_basis::derivateLaguerrePolynsSub
//================================

void Laguerre_polynomial_basis::derivateLaguerrePolynsSub(
    Laguerre_polynomial &LP) {
  std::array<Laguerre_polynomial, 3> arLP =
  {{
    Laguerre_polynomial( -LP.a,     LP.c, LP.n-1, LP.p+1, LP.px),
    Laguerre_polynomial(LP.a*LP.c,  LP.c, LP.n,   LP.p,   LP.px),
    Laguerre_polynomial(LP.a*LP.px, LP.c, LP.n,   LP.p,   LP.px-1)
  }};
  if (LP.px == 0)
    arLP[2].n = -1;
  for (auto &x : arLP) {
      if (x.n < 0)
        continue;
      int id  = x.n+x.p+x.px;
      auto lb = nextDeriv_->polynoms_.lower_bound(id);
      auto ub = nextDeriv_->polynoms_.upper_bound(id);
      bool included = false;
      for ( ; lb != ub; ++lb) {
          // overloaded operator==( L_p &a, &b) do not compare a.a with b.a
          if (lb->second == x) {
              lb->second.a += x.a;
              included = true;
              break;
            }
        }
      if (!included)
        nextDeriv_->polynoms_.emplace(id, std::move(x));
    }
}

//================================
// Laguerre_polynomial_basis::derivateLaguerrePolyns
//================================
// public
void Laguerre_polynomial_basis::derivateLaguerrePolyns(
    const size_t derivateOrder) {
  if (derivateOrder == 0)
    return;
  if (polynoms_.empty())
    throw LPBexception_empty();
  if (derivateOrder > 30) {
      throw LPBexception_calculation("Too greate order of derivation");
    }
  derivateLaguerrePolyns(1, derivateOrder);
}

// private
void Laguerre_polynomial_basis::derivateLaguerrePolyns(
    const size_t derivateOrder, const size_t derivateOrderFinish) {
  nextDeriv_  = std::unique_ptr<Laguerre_polynomial_basis>(
      new Laguerre_polynomial_basis(derivateOrder));
  auto LPB_it = polynoms_.begin();
  for ( ; LPB_it != polynoms_.end(); ++LPB_it)
    derivateLaguerrePolynsSub(LPB_it->second);
  if (derivateOrder < derivateOrderFinish)
    nextDeriv_->derivateLaguerrePolyns(derivateOrder + 1, derivateOrderFinish);
}

//================================
// Laguerre_polynomial_basis::integrateLaguerrePolynsBasis
//================================

double Laguerre_polynomial_basis::integrateLaguerrePolynsBasis(
    Laguerre_polynomial_basis &LPBf, Laguerre_polynomial_basis &LPBph,
    const size_t derivOrder) {
  if (derivOrder > 30)
    throw LPBexception_calculation("Too greate order of derivation");
  Laguerre_polynomial_basis *tempLPBf = &LPBf;
  double result = 0;
  while (true) {
      if (derivOrder == tempLPBf->depth_)
        break;
      if (tempLPBf->nextDeriv_ == nullptr)
        tempLPBf->derivateLaguerrePolyns(tempLPBf->depth_ + 1, derivOrder);
      tempLPBf = tempLPBf->nextDeriv_.get();
    }
  auto LPBf_it  = tempLPBf->polynoms_.begin();
  for ( ; LPBf_it != tempLPBf->polynoms_.end(); ++LPBf_it) {
      auto LPBph_it = LPBph.polynoms_.begin();
      for ( ; LPBph_it != LPBph.polynoms_.end(); ++LPBph_it)
        result += integrateLaguerrePair(LPBf_it->second, LPBph_it->second);
    }
  return result;
}


double Laguerre_polynomial_basis::integrateLaguerrePolynsBasis(
    Laguerre_polynomial_basis &LPBf, Laguerre_polynomial_basis &LPBph,
    const size_t derivOrder, const double *U, double dh, size_t JMax) {
  if ((U == nullptr) || (U == NULL))
    throw LPBexception_calculation("integrateLaguerrePolynsBasis get nullptr");
  if (derivOrder > 30)
    throw LPBexception_calculation("Too greate order of derivation");
  size_t dhU      = 4,
         dhUm2    = dhU << 1;

  double result   = 0.0,
         dx       = dhU*dh,
         dxm2     = dx*2,
         x;
         JMax     -= dhUm2;
  Laguerre_polynomial_basis *tempLPBf = &LPBf;
  while (true) {
      if (derivOrder == tempLPBf->depth_)
        break;

      // if next derivate was not calculated - calculate
      if (tempLPBf->nextDeriv_ == nullptr) {
          tempLPBf->derivateLaguerrePolyns(tempLPBf->depth_ + 1, derivOrder);
        }
      tempLPBf = tempLPBf->nextDeriv_.get();
    }
  for (auto i = LPBph.polynoms_.begin(); i != LPBph.polynoms_.end(); ++i) {
      Laguerre_polynomial *ph = &(i->second);
      if ((ph->n < 0) || (std::abs(ph->a) < ACCURACY))
        continue;
      for (auto j = tempLPBf->polynoms_.begin();
                j != tempLPBf->polynoms_.end(); ++j) {
          Laguerre_polynomial *f = &(j->second);
          if ((f->n < 0) || (std::abs(f->a) < ACCURACY))
            continue;
          x = 0.0;
          for (size_t k = 0; k < JMax; k += dhUm2) {
              result += dx * 0.66666667 *
                  (U[k]* L_appr(x, f->n, f->p)* L_appr(x, ph->n, ph->p)*
                  f->a*ph->a * std::pow(x, f->px+ph->px)* std::exp(-x) +
                  U[k+dhU]* L_appr(x+dx, f->n, f->p)* L_appr(x+dx, ph->n, ph->p)*
                  f->a*ph->a *  std::pow(x+dx, f->px+ph->px)* std::exp(-x-dx) +
                  U[k+dhUm2]*L_appr(x+dxm2, f->n, f->p)* L_appr(x+dxm2, ph->n, ph->p)*
                  f->a*ph->a * std::pow(x+dxm2, f->px+ph->px)* std::exp(-x-dxm2));
              x += dxm2;
            }
        }
    }
  return result;
}


//================================
// Laguerre_polynomial_basis::equatingLaguerreXandPX
//================================

void Laguerre_polynomial_basis::equatingLaguerreXandPX(
    Laguerre_polynomial &forig, Laguerre_polynomial *f) {
  if (forig.p == forig.px) {
      f[forig.n] = forig;
      for (int i = 0; i < forig.n; ++i) {
          f[i] = forig; f[i].a = 0.0;
        }
      return;
    } else if (forig.p > forig.px) {
      for (int i = 0; i < (forig.n + 1); ++i) {
          f[i].a = forig.a *
              binomialcoef(forig.p - forig.px + forig.n - i - 1, forig.n-i);
          f[i].n = i; f[i].px = forig.px;
          f[i].p = forig.px; f[i].c = forig.c;
        }
      return;
    } else {
      // count of zeros elements
      int por = forig.px - forig.p;
      for (int i = forig.n; i >= 0; --i) {
          if ((forig.n - i) <= por)
            f[i].a = forig.a *
                binomialcoef(por, forig.n - i) * std::pow(-1, forig.n - i);
          else
            f[i].a = 0.0;
          f[i].n = i; f[i].px = forig.px;
          f[i].p = forig.px; f[i].c = forig.c;
        }
    }
}

//================================
// Laguerre_polynomial_basis::integrateLaguerrePolynomials
//================================

double Laguerre_polynomial_basis::integrateLaguerrePolynomials(
    const int i, const int k, const int a) {
  double result = 0.0;
  if (i < 0)
    return result;
  if (i == k) {
      result = 1.0;
      for (int n = 1; n <= a; ++n)
        result *= (i + n);
    }
  return result;
}

//================================
// Laguerre_polynomial_basis::L_appr
//================================

double Laguerre_polynomial_basis::L_appr(
    const double x, const int n, const int a) {
  double result = 0.0;
  for (int i = 0; i <= n; ++i)
    result += std::pow(-1, i) * std::pow(x, i) * binomialcoef(n + a, n-i) / 
        factorial(i);
  return result;
}

//================================
// Laguerre_polynomial_basis::integrateLaguerrePair
//================================

double Laguerre_polynomial_basis::integrateLaguerrePair(
    Laguerre_polynomial &f, Laguerre_polynomial &ph) {
  double result = 0.0;
  if (f.n < 0)
    return result;
  std::vector<Laguerre_polynomial> Ln(f.n+1);
  std::vector<Laguerre_polynomial> Lk(ph.n+1);
  Laguerre_polynomial f1(ph.a, ph.c, ph.n, ph.p, f.px);
  equatingLaguerreXandPX(f,  &Ln[0]);
  equatingLaguerreXandPX(f1, &Lk[0]);
  int h = std::min(f.n, ph.n);
  for (int i = 0; i <= h; ++i)
    result += Ln[i].a*Lk[i].a*integrateLaguerrePolynomials(i, i, Ln[i].px);
  return result;
}

//================================
// operator==
//================================

bool operator==(const Laguerre_polynomial &lpl,
                const Laguerre_polynomial &lpr) {
  bool isEq = (std::tie(lpl.n, lpl.p, lpl.px) ==
               std::tie(lpr.n, lpr.p, lpr.px));
  if (isEq)
    isEq = (((lpl.c+ACCURACY) > lpr.c) == ((lpr.c+ACCURACY) > lpl.c));
  return isEq;
}

//================================
// operator<<
//================================

std::ostream &operator<<(std::ostream &outstr,
                         const Laguerre_polynomial &lp) {
  outstr << lp.a << "Lag(" << lp.n << ", "<< lp.p << ") * x^" << lp.px
         << " * exp(x *(" << lp.c << "))\n";
  return outstr;
}

//================================
// laguerre_polynoms ctors
//================================

Laguerre_polynomial::Laguerre_polynomial(
    double a, double c, int n, int p, int px)
  : a(a), c(c), n(n), p(p), px(px) {}

Laguerre_polynomial::Laguerre_polynomial()
  : n(-1) {}

//================================
// laguerre_polynoms swap
//================================

void Laguerre_polynomial::swap(Laguerre_polynomial &LP) noexcept {
  Laguerre_polynomial temp(std::move(LP));
  LP = std::move(*this);
  *this = temp;
}

#ifdef UNIT_TEST

//================================
// Laguerre_polynomial_basis::equatingLaguerreXandPX
//================================


std::multimap<int, Laguerre_polynomial>
    Laguerre_polynomial_basis::getPolynoms() const {
  return polynoms_;
}

//================================
// Laguerre_polynomial_basis::outputClosedForm
//================================

std::vector<double>
Laguerre_polynomial_basis::outputClosedForm() const {
  if (closedForm_.empty())
    getClosedForm();
  return closedForm_;
}

#endif  // UNIT_TEST
}  // namespace laguerre_polynoms
