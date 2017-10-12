#include "galerk_method.h"

#include "maindatastorage.h"
#include "some_math.h"
#include "wide_output.h"

#include <boost/bind.hpp>

#include <algorithm>

#define imaginaryOne std::complex<double>(0.0, 1.0)

extern bool WIDE_OUTPUT;

namespace conv_flow {
//================================
// galerkMethod ctor
//================================

galerkMethod::galerkMethod(galerkInput gi, MainDataStorage *mds)
  : polynomsPower_(gi.polynomsPower),
    countBasicFunc_(gi.countBasicFunc),
    startOrder_(gi.startOrder), mds_(mds) {
  if (mds == nullptr)
    throw galerkMethodExcept("Galerkin Method");
}

//================================
// galerkMethod::chooseEigenvalue
//================================

int galerkMethod::chooseEigenvalue(std::vector<compd> &eigenValues) {
  std::vector<compd> vec(eigenValues);

  vec.erase(std::remove_if(vec.begin(), vec.end(),
      [](const compd& c) {
          return ((c.real() < 0.0));
        }));
  vec.erase(std::remove_if(vec.begin(), vec.end(),
      [](const compd& c) {
          return ((std::abs(c.imag()/c.real()) > 1.0));
        }));

  auto el = std::max_element(vec.begin(), vec.end(),
      [](const compd& l, const compd& r) {
          return l.imag() < r.imag();
        });
  int ind = (el == vec.end()) ? -1 : int{el - vec.begin()};
  if (ind == -1)
    return ind;
  return int{std::find_if(eigenValues.begin()+ind, eigenValues.end(),
                 boost::bind(&cmplx_eq<double, double>, *el, _1, 0.0002))
                 - eigenValues.begin()};
}

//================================
// galerkMetho::setBasisFunc
//================================
// method sets Laguerre basis functions - combinate of Laguerre structs
void galerkMethod::setBasisFunc(
    laguerre_polynoms::Laguerre_polynomial_basis &fiLPB,
    laguerre_polynoms::Laguerre_polynomial_basis &siLPB,
    const double expc, const size_t polOrder, const size_t xPower) {
  if (polynomsPower_ == 2) {
      fiLPB.addToBasis(1.0, expc, polOrder, polynomsPower_, xPower);
      if (xPower != 0) {
          siLPB.addToBasis(1.0, expc, polOrder, static_cast<int>(polynomsPower_)-1,
              static_cast<int>(xPower)-1);
          return;
        }
      siLPB.addToBasis(1.0, expc, polOrder, static_cast<int>(polynomsPower_)-1, 0);
      return;
    } else if (polynomsPower_ == 0) {
      fiLPB.addToBasis(1.0, expc, polOrder, polynomsPower_, 0);
      fiLPB.addToBasis(-static_cast<int>(polOrder), expc, 1, polynomsPower_, 0);
      fiLPB.addToBasis(static_cast<int>(polOrder)-1, expc, 0, polynomsPower_, 0);
      siLPB.addToBasis(1.0, expc, static_cast<int>(polOrder)-1, polynomsPower_, 0);
      siLPB.addToBasis(-1.0, expc, 0, polynomsPower_, 0);
      return;
    } else {
      throw galerkMethodExcept("galerkinMethod take not valid"
          " parameter \"polynomsPower\"");
    }
}

//================================
// galerkMethod setResult
//================================
// method return indeces eigenvalues vector of ci(expc)
std::vector<size_t> galerkMethod::setResult(std::vector<compd> &veceigvs) {
  //  sequence of complex numbers(phase velocities) calculated with various
  //         expc ("ci" in method "GalerkMethodCalculate")

  if (veceigvs.empty())
    throw galerkMethodExcept("setResult take empty vector");
  std::vector<size_t> inds;
  size_t neigs = veceigvs.size();
  std::vector<size_t> indsnot;   // push duplicate indices
  for (size_t i = 0; i < neigs; ++i) {
      // if veceigvs[i] have duplicate in veceigvs[0] - veceigvs[i-1]
      if (std::find(indsnot.begin(), indsnot.end(), i) != inds.end())
        continue;
      // check duplicates for current element
      auto iter = veceigvs.begin()+i;
      while (1) {
          // compare with ~30% accuracy
          iter = std::find_if(iter + 1, veceigvs.end(),
              boost::bind(&cmplx_eq<double, double>,
                  *iter, _1, (*iter).real()*0.30));
          if (iter == veceigvs.end())
            break;
          indsnot.push_back(
              static_cast<size_t>(iter-veceigvs.begin()));
        }
      inds.push_back(i);
    }
  return inds;
}

//================================
// galerkMethod setTriangular
//================================
// overload of template method for "armadillo" matrices
bool galerkMethod::setTriangular(arma::cx_mat &A, const size_t n) {
  for (size_t i = 0; i < n; ++i)
    if (std::abs(A(i, i)) < 0.001)
      for (size_t z = 0; z < n; ++z) {
          if (z == i)
            continue;
          for (size_t k = 0; k < n; ++k)
            A(i, k) += A(z, k);
          if (std::abs(A(i, i)) > 0.001)
            break;
        }
  for (size_t i = 0; i < (n-1); ++i)
    for (size_t j = i+1; j < n; ++j) {
        std::complex<double> temp = A(j, i)/A(i, i);
        for (size_t k = i; k < n; ++k)
          A(j, k) -= A(i, k)*temp;
      }
  return false;
}

//================================
// galerkMethod fillMatrices
//================================
// calculate elements of matrices A and B - there is Orr-Sommerfeld equation
void galerkMethod::fillMatrices(arma::cx_mat &A, arma::cx_mat &B,
    const std::array<double, 4> &coef, const double wavnum, const double expc) {
  typedef laguerre_polynoms::Laguerre_polynomial_basis LPB_t;
  size_t nMat = countBasicFunc_,
         JMax = static_cast<size_t>(mds_->zerosEdge/mds_->dh);
  for (size_t i = 0; i < nMat; ++i) {
      LPB_t fiLPB;
      LPB_t siLPB;
      setBasisFunc(fiLPB, siLPB, expc, i+startOrder_, polynomsPower_);
      fiLPB.derivateLaguerrePolyns(4);
      siLPB.derivateLaguerrePolyns(2);
      for (size_t j = 0; j < nMat; ++j) {
          LPB_t phLPB;
          LPB_t tetLPB;
          setBasisFunc(phLPB, tetLPB, -1-expc, j+startOrder_, 0);
          A(i, j) = coef[1] * (
              LPB_t::integrateLaguerrePolynsBasis(fiLPB, phLPB, 4) -
              2.0*wavnum*wavnum *
              LPB_t::integrateLaguerrePolynsBasis(fiLPB, phLPB, 2) +
              std::pow(wavnum, 4.0) *
              LPB_t::integrateLaguerrePolynsBasis(fiLPB, phLPB, 0))
              / (imaginaryOne*wavnum) - coef[0] * (
              LPB_t::integrateLaguerrePolynsBasis(
                  fiLPB, phLPB, 2, &(mds_->F1[0]), mds_->dh, JMax) -
              LPB_t::integrateLaguerrePolynsBasis(
                  fiLPB, phLPB, 0, &(mds_->F3[0]), mds_->dh, JMax) -
              LPB_t::integrateLaguerrePolynsBasis(
                  fiLPB, phLPB, 0, &(mds_->F1[0]), mds_->dh, JMax)*
              wavnum*wavnum);

          A(i+nMat, j) = coef[2] *
              LPB_t::integrateLaguerrePolynsBasis(siLPB, phLPB, 1)
              /(imaginaryOne*wavnum);

          A(i+nMat, j+nMat) = coef[3] * (
              LPB_t::integrateLaguerrePolynsBasis(siLPB, tetLPB, 2) -
              wavnum*wavnum*
              LPB_t::integrateLaguerrePolynsBasis(siLPB, tetLPB, 0))
              / (imaginaryOne*wavnum) -
              coef[0] *
              LPB_t::integrateLaguerrePolynsBasis(
                  siLPB, tetLPB, 0, &(mds_->F1[0]), mds_->dh, JMax);

          A(i, j+nMat) = coef[0] *
              LPB_t::integrateLaguerrePolynsBasis(
                  fiLPB, tetLPB, 0, &(mds_->H1[0]), mds_->dh, JMax);
          B(i, j) = -LPB_t::integrateLaguerrePolynsBasis(fiLPB, phLPB, 2) +
              wavnum*wavnum *
              LPB_t::integrateLaguerrePolynsBasis(fiLPB, phLPB, 0);
          B(i+nMat, j+nMat) =
              -LPB_t::integrateLaguerrePolynsBasis(siLPB, tetLPB, 0);
        }
    }
}

//================================
// galerkMethod GalerkinMethodCalculate
//================================
// calculate eigenvalues and eigenvectors, select best (or valid) of thems
void galerkMethod::GalerkMethodCalculate(const  double Re, const double wavnum,
    std::vector<std::array<compd, 4>> &baseValues) {
  typedef std::array<compd, 4> params_t;
  const double Pr = mds_->Pr;

  std::vector<params_t> basicvaluesOfExpc,
                        tempResult;
  size_t nOrder  = countBasicFunc_,
         nMatrix = nOrder << 1;
  arma::cx_mat A (nMatrix, nMatrix),
               B (nMatrix, nMatrix, arma::fill::zeros),
               AA(nMatrix, nMatrix),
               BB(nMatrix, nMatrix),
               Q (nMatrix, nMatrix),
               Z (nMatrix, nMatrix);
  basicvaluesOfExpc.reserve(nMatrix);
  std::vector<compd> eigenValues(nMatrix);
  int niceInd;
  std::array<double, 4> coef 
  {{
    1.0, 1.0/Re, 1.0/Re, 1.0/Re/Pr
  }};
    // variable ci sets rate of decrease for basis functions
  double ci = 0.15;                           
    // basic functions f has a type:
    // a*x^px*Lag(^p,n)*exp(-ci*x)
    // or a*(Lag(^0,n+2)-n*Lag(^0,1)+(i-1)*Lag(^0,0))*exp(-ci*x)
  if ((polynomsPower_ == 0) && (startOrder_ < 2))
    startOrder_ = 2;
  std::vector<compd> phvelos;

  while (ci < 0.8) {
      ci += 0.05;
      fillMatrices(A, B, coef, wavnum, -ci);
      arma::qz(AA, BB, Q, Z, A, B);
      for (size_t i = 0; i < nMatrix; ++i) {
          eigenValues[i] = AA(i, i)/BB(i, i);
          if (!std::isfinite(eigenValues[i].imag())
              || !std::isfinite(eigenValues[i].real()))
            continue;
          basicvaluesOfExpc.push_back(
              {{
                std::complex<double>(0.0, 0.0),
                std::complex<double>(0.0, 0.0),
                std::complex<double>(0.0, 0.0),
                eigenValues[i]
              }});
          getBasicValues(A, B, basicvaluesOfExpc.back(), -ci);
        }
      niceInd = chooseEigenvalue(eigenValues);
      if (niceInd != -1) {
          tempResult.push_back(basicvaluesOfExpc[niceInd]);
          phvelos.push_back(eigenValues[niceInd]);
        }
      if (WIDE_OUTPUT) {
          std::cerr << "diagnostic message\n";
          outGalerkin(eigenValues, Re, wavnum, -ci);
        }
      std::vector<params_t>().swap(basicvaluesOfExpc);
    }
  std::vector<size_t> inds(setResult(phvelos));
  for (size_t i = 0; i < inds.size(); ++i )
    baseValues.push_back(tempResult[inds[i]]);
  if (baseValues.empty())
    throw galerkMethodExcept("Galerkin method return empty vector;");
}

//================================
// galerkMethod getBasicValues
//================================
// method return array {f''(0),f'''(0),s'(0)} for this eigenvalue
void galerkMethod::getBasicValues(arma::cx_mat &A, arma::cx_mat &B,
    std::array<compd, 4> &bv, const double expc) {
  typedef laguerre_polynoms::Laguerre_polynomial_basis LPB_t;
  compd phvel = bv[3],
        f2    = 0.0,
        f3    = 0.0,
        s1    = 0.0;
  size_t matOrder = countBasicFunc_ << 1;
  arma::cx_mat D(matOrder, matOrder);

  for (size_t i = 0; i < matOrder; i++) {
      for (size_t j = 0; j < matOrder; j++)
        D(i, j) = A(j, i)-phvel*B(j, i);
    }
  setTriangular(D, matOrder);
  arma::cx_vec y(matOrder, arma::fill::zeros);
  y(matOrder-1) = 1.0;
  D(matOrder-1, matOrder-1) = 1.0;
  arma::cx_vec coef;

// function "solve" can throw runtime error if matrix D is singular
  try {
      coef = arma::solve(D, y);
    } catch (std::runtime_error &e) {
      std::cerr << "Arma error for c " << bv[3] << "\n"
                << e.what() << std::endl;

  // in this case method return startcondition 0.0,0.0,0.0 - not good, not bad
  //  - perhaps, nachtsmethod will not provide result for this startcond
      return;
    }
  // define stream and temperature functions
  for (size_t i = 0; i < countBasicFunc_; ++i) {
      LPB_t LPBf;
      LPB_t LPBs;
      setBasisFunc(LPBf, LPBs, expc, i+startOrder_, polynomsPower_);
      LPBf.derivateLaguerrePolyns(3);
      LPBs.derivateLaguerrePolyns(1);
      LPB_t *tempLPB = LPBf.nextDerivate()->nextDerivate();
      f2 = tempLPB->getValue(coef[i], 0.0);
      tempLPB = tempLPB->nextDerivate();
      f3 = tempLPB->getValue(coef[i], 0.0);
      tempLPB = LPBs.nextDerivate();
      s1 = tempLPB->getValue(coef[i+countBasicFunc_], 0.0);
    }
  // set aproximation baseValues
  bv[0] = f2;
  bv[1] = f3;
  bv[2] = s1;
}

//================================
// galerkMethodInput struct
//================================

galerkInput::galerkInput(const size_t pp,
                         const size_t mo,
                         const size_t so)
  : polynomsPower(pp), countBasicFunc(mo), startOrder(so) {}

#ifdef UNIT_TEST

//================================
// galerkMethod::chooseEigenvalue_test
//================================

int galerkMethod::chooseEigenvalue_test(
    std::vector<compd> &vec_ev) const {
  return chooseEigenvalue(vec_ev);
}

//================================
// galerkMethod::setResult_test
//================================

std::vector<size_t> galerkMethod::setResult_test(
    std::vector<compd> &veceigvs) const {
  return setResult_test(veceigvs);
}

#endif
}  // namespace conv_flow
