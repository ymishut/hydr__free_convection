#include "nachts_method.h"

#include "maindatastorage.h"
#include "some_math.h"
#include "stabilitypoints.h"
#include "wide_output.h"

#include <armadillo>
#include <boost/numeric/odeint.hpp>

#include <algorithm>
#include <functional>
#include <iostream>
#include <mutex>
#include <utility>

#define imaginaryOne std::complex<double>(0.0, 1.0)

extern bool   WITH_IMPLICIT_METHODS;
extern bool   WIDE_OUTPUT;
extern int    ADAMS_METHOD_ORDER;
extern double ACCURACY;

std::mutex mtx;

namespace conv_flow {
//================================
// nachtsMethod::AdamsMethods
//================================

void nachtsMethod::AdamsMethods(
    flowparameters &fp, const double dx, size_t IMax,
    flowparameters *fc) {
  std::array<double, 3> fd;
  IMax -= 1;
  auto intStep = [this, &fp, fc, &fd, dx](adamFunc_nach f, int i) {
        for (int j = 6, k = 7; j >= 0; --j, --k) {
            if (j == 4)
              continue;
            fp.fisi[j][i+1] = f(dx, &fp.fisi[k][0], fp.fisi[j][i], i);
          }

        fp.fisi[4][i+1] = 2.0*wavnum_*wavnum_*fp.fisi[2][i+1] -
            std::pow(wavnum_, 4)*fp.fisi[0][i+1] -
            fp.fisi[6][i+1]*coef_[2]/coef_[1] + imaginaryOne*wavnum_*
            ((coef_[0]*fd[0] - phvel_)*(fp.fisi[2][i+1] -
            wavnum_*wavnum_*fp.fisi[0][i+1]) -
            coef_[0]*fd[1]*fp.fisi[0][i+1]) / coef_[1];
        fp.fisi[7][i+1] = wavnum_*wavnum_*fp.fisi[5][i+1] +
            imaginaryOne*wavnum_*((coef_[0]*fd[0] - phvel_)*fp.fisi[5][i+1] -
            coef_[0]*fd[2]*fp.fisi[0][i+1]) / coef_[3];

        if (fc != NULL) {
            fp.fisi[4][i+1] -= imaginaryOne*wavnum_*(fc->fisi[2][i+1] -
                wavnum_*wavnum_*fc->fisi[0][i+1])/coef_[1];
            fp.fisi[7][i+1] -= imaginaryOne*wavnum_*fc->fisi[5][i+1]/coef_[3];
          }
    };

  if (WITH_IMPLICIT_METHODS)
    for (size_t i = ADAMS_METHOD_ORDER; i < IMax; ++i) {
// basic flow was calculated with double accuracy
        fd[0] = mds_->F1[(i << 1)+2]; fd[1] = mds_->F3[(i << 1)+2]; 
        fd[2] = mds_->H1[(i << 1)+2];
        intStep(adamsBashforthMethod<compd>, i);
        intStep(adamsMoultonMethod<compd>, i);
      }
  else
    for (size_t i = ADAMS_METHOD_ORDER; i < IMax; ++i) {
        fd[0] = mds_->F1[(i << 1)+2]; fd[1] = mds_->F3[(i << 1)+2];
        fd[2] = mds_->H1[(i << 1)+2];
        intStep(adamsBashforthMethod<compd>, i);
      }
}

//================================
// nachtsMethod ctor
//================================

nachtsMethod::nachtsMethod(const MainDataStorage *mds, const bool uniqueResult)
  : mds_(mds), uniq_(uniqueResult) {
  for (auto &x : parametersABC_) {
      x = std::unique_ptr<flowparameters>(
          new flowparameters(mds_->JMax >> 1));
    }
}

//================================
// nachtsMethod::NachtsMethodCalculate
//================================

bool nachtsMethod::NachtsMethodCalculate(
    const nachtsInput &inp, const std::array<compd, 4> &baseValues,
    StabilityPoints &sp) {
  typedef boost::numeric::odeint::runge_kutta4
      <nachtsMethod::state_type_nach> stepper_type_nach;
  size_t IMax = static_cast<size_t>(static_cast<size_t>
      (mds_->zerosEdge/mds_->dh) >> 1);
  if (IMax > (mds_->JMax >> 1))
    IMax = mds_->JMax >> 1;
  const double dx = mds_->dh*2;
// for Newton method (calculating f'''[0], s'[0] and c)
  const size_t MAX_ITERATION = 100;

  parametersABC_[3]->fisi[2][0] = baseValues[0];
  parametersABC_[3]->fisi[3][0] = baseValues[1];
  parametersABC_[3]->fisi[6][0] = baseValues[2];

  std::unique_ptr<compd[]> corrections(new compd[3]);
  compd **Acorr = new compd*[3];
  for (size_t i = 0 ; i < 3; i++)
    Acorr[i] = new compd[4];

  // coefficients of stability equtions.
  coef_[0] = 1.0;
  coef_[1] = coef_[2] = compd(1.0, 0.0)/inp.Re;
  coef_[3] = coef_[1]/mds_->Pr;

  double wavnumMax = inp.wavnum_max,
         dwavnum   = inp.dwavnum,
         limitEdge = (ADAMS_METHOD_ORDER+1)*dx;

  compd   cMax,
          ctemp = baseValues[3];
  size_t  currentiter = 0,
          edgeL  = static_cast<size_t>(mds_->zerosEdge/dx)-1;
  stPoint tempPoint(inp.Re, inp.wavnum_min, baseValues[3],
      baseValues[1], baseValues[2]);

  compd fsiabcedge[32];
  state_type_nach x;
  for (size_t i = 0; i < 32; i++)
    fsiabcedge[i] = compd(0.0, 0.0);
  cMax = ctemp*1.2;
  wavnum_ = inp.wavnum_min;

// set integrate functor and observer
  nachtsheinSystem nachSys(false, this);
  pushresult_observer nachObs(3, this);

// MAIN LOOP BY WAVENUMBERS   wavenumbers was setted by client
// as well as interval phasevelocity
  while (wavnum_ <= wavnumMax) {
      ctemp  = 0.64*cMax;
      phvel_ = ctemp;
      double inaccuracy;

    // INNER LOOP BY PHASE VELOCITY  phase velocity - also by client,
    // or as Galerkin method results
      while (phvel_.real() < cMax.real()) {
          for (size_t i = 0; i < 3; i++)
            corrections[i] = 0.0;
          bool pointfound = false;
        // initialization of start data by previous result 
        // or (if it is first iteration) by input data
          parametersABC_[3]->fisi[3][0] = tempPoint.f3b;
          parametersABC_[3]->fisi[6][0] = tempPoint.s1b;

          // ANOTHER ONE LOOP
          // MAX_ITERATION  is  program configuration
          for ( ; currentiter < MAX_ITERATION; ++currentiter) {
              parametersABC_[3]->fisi[3][0] += corrections[0];
              parametersABC_[3]->fisi[6][0] += corrections[1];
              phvel_ += corrections[2];

              // integtation of main equations
              x[0] = 0.0; x[1] = 0.0;                   // f and f'
              x[2] = parametersABC_[3]->fisi[2][0];     // f''
              x[3] = parametersABC_[3]->fisi[3][0];     // f'''
              x[4] = 0.0; x[5] = parametersABC_[3]->fisi[6][0];  // s and s'
              flowparameters *fcref = 0;
              nachObs.nOrder = 3;
              boost::numeric::odeint::integrate_const(stepper_type_nach(),
                  nachSys, x, 0.0, limitEdge, dx, nachObs);
              AdamsMethods(*parametersABC_[3], dx, IMax, fcref);

              // integtation of system for caltulate corrections
              for (size_t i = 0; i <3 ; ++i) {
                  x[0] = 0.0; x[1] = 0.0;         // as above
                  x[2] = 0.0; x[3] = 0.0;
                  x[4] = 0.0; x[5] = 0.0;
                  if (i == 0) {
                      x[3] = 1.0;
                    } else if (i == 1) {
                      x[5] = 1.0;
                    } else {
                      nachSys.phVelocDerivated = true;
                      fcref = parametersABC_[3].get();
                    }
                  nachObs.nOrder = i;
                  boost::numeric::odeint::integrate_const(stepper_type_nach(),
                      nachSys, x, 0.0, limitEdge, dx, nachObs);
                  AdamsMethods(*parametersABC_[i], dx, IMax, fcref);
                }
              fcref = 0;
              nachSys.phVelocDerivated = false;

              for (size_t i = 0; i < 8; ++i)
                fsiabcedge[i] = parametersABC_[3]->fisi[i][edgeL];
              for (size_t j = 0; j < 3; ++j) {
                  size_t jtemp = j << 3;
                  for (size_t i = 0; i < 8; ++i)
                    fsiabcedge[8+i+jtemp] = parametersABC_[j]->fisi[i][edgeL];
                }

              //  variables from publication
              // have not meaning, but it is convenient
              compd beta  = std::sqrt(wavnum_*wavnum_-
                        imaginaryOne*wavnum_*phvel_/coef_[1]),
                    gamma = std::sqrt(wavnum_*wavnum_-
                        imaginaryOne*wavnum_*phvel_/coef_[3]);

              /* if ES is equations of Nachtsheim system, then :
                 Acorr 0 column - d(ES)/da   a = f'''(0)
                 Acorr 1 column - d(ES)/db   b = s'(0)
                 Acorr 2 column - d(ES)/dc   c = phase velocity   
                    - for this system we have slight differences equations
                 Acorr 3 column - ES */
              size_t jcorr;
              for (size_t j = 0; j < 4; ++j) {
                  if (j > 0)
                    jcorr = j-1;
                  else
                    jcorr = 3;
                  size_t jtemp = j << 3;
                  Acorr[0][jcorr] = fsiabcedge[6+jtemp] + 
                      gamma*fsiabcedge[5+jtemp];

                  Acorr[1][jcorr] = fsiabcedge[3+jtemp] -
                      wavnum_*wavnum_*fsiabcedge[1+jtemp] +
                      beta* (fsiabcedge[2+jtemp] -
                      wavnum_*wavnum_*fsiabcedge[0+jtemp]) +
                      coef_[2]*gamma*fsiabcedge[5+jtemp] /
                      (coef_[1]* (gamma+beta));

                  Acorr[2][jcorr] = fsiabcedge[3+jtemp] +
                      wavnum_*fsiabcedge[2+jtemp] - beta*beta *
                      (fsiabcedge[1+jtemp]+wavnum_*fsiabcedge[0+jtemp]) +
                      coef_[2]*gamma*fsiabcedge[5+jtemp] /
                      (coef_[1]* (gamma+wavnum_));
                }

              Acorr[0][2] += -0.5*imaginaryOne*wavnum_*fsiabcedge[5] /
                  (coef_[3]*gamma);

              Acorr[1][2] += -imaginaryOne*wavnum_*0.5 *
                  (fsiabcedge[2]-wavnum_*wavnum_*fsiabcedge[0]) /
                  (beta*coef_[1]) -
                  imaginaryOne*(mds_->Pr-1.0)*wavnum_*wavnum_*wavnum_ /
                  (coef_[1]*2.0*gamma*beta*std::pow(gamma+beta,2));

              Acorr[2][2] += imaginaryOne*wavnum_ * 
                  (fsiabcedge[1]+wavnum_*fsiabcedge[0] -
                  0.5*wavnum_*mds_->Pr*fsiabcedge[5] /
                  (gamma*std::pow(gamma+wavnum_,2.0)))/coef_[1];

              for (size_t j = 0; j < 3; ++j)
                Acorr[j][3] *= -1.0;

              bool isSingular = false;
              // check result
              inaccuracy = std::abs(Acorr[0][3])+std::abs(Acorr[1][3])
                  +std::abs(Acorr[2][3]);
              gaussMethod(Acorr, corrections.get(), 3, isSingular);
              // could not finish calculate, try next approximation for phase velocity
              if (std::abs(corrections[0]) > 10e6)
                break;
              if (ACCURACY > inaccuracy) {
                  pointfound = true;
                  break;
                }
            }  // for ( ; currentiter < MAX_ITERATION; ++currentiter) 
          if (pointfound) {
              tempPoint.wavnum = wavnum_;
              tempPoint.phvel = (mds_->flowID == FlowID::SIEGEL) ?
                  phvel_*coef_[1] : phvel_;
              tempPoint.f3b = parametersABC_[3]->fisi[3][0] /
                  parametersABC_[3]->fisi[2][0];
              tempPoint.s1b = parametersABC_[3]->fisi[6][0] / 
                  parametersABC_[3]->fisi[2][0];

              if (WIDE_OUTPUT) {
                  compd norma = fsiabcedge[0];
                  auto normalize = [norma](std::vector<compd> &vec) {
                    for_each(vec.begin(), vec.end(),
                        std::bind2nd(std::divides<compd>(), norma));
                  };
                  normalize(parametersABC_[3]->fisi[0]);
                  normalize(parametersABC_[3]->fisi[1]);
                  normalize(parametersABC_[3]->fisi[5]);
                  // output normalized parameters
                  std::lock_guard<std::mutex> lock(mtx);
                  outWaves(parametersABC_[3]->fisi[0], parametersABC_[3]->fisi[1],
                      parametersABC_[3]->fisi[5], 2.0*mds_->dh);
                }
              // calculate only first phase velocity
              //   for animation
              if (uniq_)
                break;
              sp.addPoint(tempPoint);
            }
          ctemp += std::real(cMax)*0.08;
          phvel_ = ctemp;
          currentiter = 0;
        }  // while (phvel_.real() < cMax.real()) 
      wavnum_ += dwavnum;
    }  // while (wavnum_ <= wavnumMax) {
  {
    std::lock_guard<std::mutex> lock(mtx);
    outStPoints(mds_);
  }
  for (size_t i = 0; i < 3; i++)
    delete[] Acorr[i];
  delete[] Acorr;
  return false;
}

//================================
// nachtsMethod::getFlowparameters
//================================

const std::vector<nachtsMethod::compd> &nachtsMethod::getFlowparameters(
    size_t i) const {
  if (i > 7) {
    std::cerr << " method getFlowparameters take not correct arg 'i'\n";
    i = 0;
  }
  return parametersABC_[3]->fisi[i];
}

//================================
// nachtsMethod::getPhasevelocity
//================================

nachtsMethod::compd nachtsMethod::getPhasevelocity() const {
  return phvel_;
}

//================================
// nachtsMethod::pushresult_observer struct
//================================

nachtsMethod::pushresult_observer::pushresult_observer(
    size_t nOrd, nachtsMethod *nm)
  : nOrder(nOrd), nm(nm) {}

void nachtsMethod::pushresult_observer::operator()(
    const state_type_nach &x, double t) {
  double wavnum_ = nm->wavnum_;
  size_t n   = static_cast<size_t>(t/nm->mds_->dh/2.0),
         iBF = n << 1;
  for (size_t i = 0; i < 4; ++i)
    nm->parametersABC_[nOrder]->fisi[i][n] = x[i];

  nm->parametersABC_[nOrder]->fisi[5][n] = x[4];
  nm->parametersABC_[nOrder]->fisi[6][n] = x[5];
  nm->parametersABC_[nOrder]->fisi[4][n] = 2.0*wavnum_*wavnum_*x[2] -
      pow(wavnum_, 4)*x[0] - x[5]*nm->coef_[2]/nm->coef_[1] +
      imaginaryOne * wavnum_* ((nm->coef_[0]*nm->mds_->F1[iBF] -
      nm->phvel_)*(x[2] - wavnum_*wavnum_*x[0]) -
      nm->coef_[0]*nm->mds_->F3[iBF]*x[0]) / nm->coef_[1];
  nm->parametersABC_[nOrder]->fisi[7][n] = wavnum_*wavnum_*x[4] +
      imaginaryOne*wavnum_*nm->mds_->Pr *
      ((nm->coef_[0]*nm->mds_->F1[iBF] -nm->phvel_)*x[4] -
      nm->coef_[0]*nm->mds_->H1[iBF]*x[0])/nm->coef_[1];
  if (nOrder == 2) {
      nm->parametersABC_[nOrder]->fisi[4][n] -= imaginaryOne*wavnum_*(
          nm->parametersABC_[3]->fisi[2][iBF] -
          wavnum_*wavnum_*nm->parametersABC_[3]->fisi[0][iBF]) / nm->coef_[1];
      nm->parametersABC_[nOrder]->fisi[7][n] -= imaginaryOne*wavnum_ *
        nm->parametersABC_[3]->fisi[6][iBF]/nm->coef_[3];
    }
}

//================================
// nachtsMethod::nachtsheinSystem struct
//================================

nachtsMethod::nachtsheinSystem::nachtsheinSystem(
    const bool pvd, nachtsMethod *nm)
  : phVelocDerivated(pvd), nm(nm) {}

void nachtsMethod::nachtsheinSystem::operator()(
    const nachtsMethod::state_type_nach &x,
    nachtsMethod::state_type_nach &dxdt, double t) {
  double wavnum_ = nm->wavnum_,
         dx      = nm->mds_->dh*2.0;
  size_t iBF = static_cast<size_t>(t/nm->mds_->dh);
  dxdt[0] = x[1];
  dxdt[1] = x[2]+0.5*x[1]*dx;
  dxdt[2] = x[3]+0.5*x[2]*dx+0.16666667*x[1]*dx*dx;
  dxdt[3] = 2.0*wavnum_*wavnum_*x[2] - std::pow(wavnum_, 4)*x[0] -
      x[5]*nm->coef_[2]/nm->coef_[1] + imaginaryOne*wavnum_ *
      ((nm->coef_[0]*nm->mds_->F1[iBF] - nm->phvel_)*(x[2] -
      wavnum_*wavnum_*x[0]) - nm->coef_[0]*nm->mds_->F3[iBF]*x[0]) /
      nm->coef_[1];
  dxdt[4] = x[5];
  dxdt[5] = wavnum_*wavnum_*x[4]+imaginaryOne*wavnum_ *
      ((nm->coef_[0]*nm->mds_->F1[iBF] - nm->phvel_)*x[4] -
      nm->coef_[0]*nm->mds_->H1[iBF]*x[0])/nm->coef_[3];
  if (phVelocDerivated) {
      dxdt[3] -= imaginaryOne*wavnum_* (nm->parametersABC_[3]->fisi[2][iBF]
          -wavnum_*wavnum_*nm->parametersABC_[3]->fisi[0][iBF]) / nm->coef_[1];
      dxdt[5] -= imaginaryOne*wavnum_* nm->parametersABC_[3]->fisi[6][iBF] /
          nm->coef_[3];
    }
}

//================================
// flowparameters struct
//================================

flowparameters::flowparameters(size_t n) {
  for (auto &x : fisi)
    x.assign(n, 0.0);
}

//================================
// nachtsInput struct
//================================

nachtsInput::nachtsInput(double Re, double REM, double dRe,
    double w_min, double w_max, double dw)
  : Re(Re), Re_max(REM), dRe(dRe),
    wavnum_min(w_min), wavnum_max(w_max), dwavnum(dw) {}
}  // namespace conv_flow
