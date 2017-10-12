#include "basicflow.h"

#include "galerk_method.h"
#include "nachts_method.h"
#include "stabilitypoints.h"
#include "some_math.h"
#include "wide_output.h"

#include <armadillo>
#include <boost/numeric/odeint.hpp>

#include <algorithm>
#include <cmath>
#include <exception>
#include <functional>
#include <numeric>
#include <iostream>
#include <tuple>

extern int    ADAMS_METHOD_ORDER;
extern bool   WITH_IMPLICIT_METHODS;
extern bool   WIDE_OUTPUT;
extern double ACCURACY;

typedef std::array<double, 5> state_type_bf;
typedef boost::numeric::odeint::runge_kutta4 <state_type_bf> stepper_type_bf;

namespace conv_flow {
//================================
// PolhausenFlow ctor
//================================

PolhausenFlow::PolhausenFlow(
    const double Pr, size_t Jm, const double zerosEdge, double dh)
  : FlowParameters(FlowID::POLHAUSEN, Pr, Jm, zerosEdge, dh) {
  // initialization error
  if (mdsp_ == nullptr)                     
    throw;
  try {
      calculate();
      // calculation error
    } catch (std::exception &e) {
      std::cerr << e.what()<< std::endl;
      throw;
    }
  if (WIDE_OUTPUT)
    outBasicflow(mdsp_.get());
}

//================================
// PolhausenFlow::calculate()
//================================

//  Хитрый метод из книжки Ц. На :
// 1) Вместо условия F=F'=H=0 на бесконечности
// он явно принимает значение границы х
// на которой, по мнению исследователя,
// эти условия уже соблюдены и сводит
// крвевую задачу к задаче Коши

// 2) Чтобы проверить корректонсть выбора границы
// метод расчитывает n точек за границей х -
// в идеале на них F=F'=H=0
void PolhausenFlow::calculate() {
  const double dh         = mdsp_->dh,
               Pr         = mdsp_->Pr;

  const size_t JMax       = size_t(mdsp_->zerosEdge/dh),
               N_MATRIX   = 5,
               K_MATRIX   = 8;

  arma::mat alf(N_MATRIX, N_MATRIX);
  arma::mat alfInv(N_MATRIX, N_MATRIX);

  float bw   [N_MATRIX];

  float ksi  [JMax][K_MATRIX],
        bet  [JMax][K_MATRIX],
        rArr [JMax][K_MATRIX],
        W    [JMax][N_MATRIX],
        D    [JMax][N_MATRIX],
        ans  [N_MATRIX][N_MATRIX];

  float A[JMax][N_MATRIX][N_MATRIX],
        B[JMax][N_MATRIX][N_MATRIX],
        C[JMax][N_MATRIX][N_MATRIX],
        G[JMax][N_MATRIX][N_MATRIX];

  for (size_t j = 0; j < JMax; ++j)
    for (size_t i = 0; i < N_MATRIX; ++i)
      D[j][i] = 0.0;

  for (size_t j = 0; j < JMax; ++j) {
      F[j]  = 0.1;
      H[j]  = 0.0;
      F1[j] = 0.0;
    }
  for (size_t j = 1; j < JMax; ++j) {
      F3[j] = 0.0;
      F2[j] = 0.1;
      H1[j] = 0.0;
    }
  F2[0] = 0.7;
  F3[0] = -1.0;
  H1[0] = -0.8;

  for (size_t z = 0; z < 8; ++z) {
      F[0]  = 0.0;
      F1[0] = 0.0;
      H[0]  = 1.0;
      F2[0] += D[1][0];
      H1[0] += D[1][1];
      for (size_t j = 1; j < JMax; ++j) {
          F[j]  += D[j][2];
          F2[j] += D[j][3];
          H1[j] += D[j][4];
          if (j < JMax-1) {
              F1[j] += D[j+1][0];
              H[j]  += D[j+1][1];
            }
          ksi[j][0] =  1.0+3.0*dh*(F[j]+F[j-1])/4.0;
          ksi[j][1] = -1.0+3.0*dh*(F[j]+F[j-1])/4.0;
          ksi[j][2] = ksi[j][3] = 3.0*dh*(F2[j]+F2[j-1])/4.0;
          ksi[j][4] = ksi[j][5] = -dh*(F1[j]+F1[j-1]);
          ksi[j][6] = ksi[j][7] = dh/2;
          bet[j][0] =  1.0+3.0*Pr*dh*(F[j]+F[j-1])/4.0;
          bet[j][1] = -1.0+3.0*Pr*dh*(F[j]+F[j-1])/4.0;
          bet[j][2] = bet[j][3] = 3.0*Pr*dh*(H1[j]+H1[j-1])/4.0;
          bet[j][4] = bet[j][5] = bet[j][6] = bet[j][7] = 0.0;

          rArr[j][0] =  F[j-1]- F[j]+dh*(F1[j]+F1[j-1])/2.0;
          rArr[j][1] = F1[j-1]-F1[j]+dh*(F2[j]+F2[j-1])/2.0;
          rArr[j][2] =  H[j-1]- H[j]+dh*(H1[j]+H1[j-1])/2.0;
          rArr[j][3] = F2[j-1]-F2[j]-dh*(H[j]+H[j-1])/2.0 -
              3.0*dh*(F[j]+F[j-1])*(F2[j]+F2[j-1])/4.0 +
              2.0*dh*std::pow((F1[j]+F1[j-1])/2.0, 2.0);
          rArr[j][4] = H1[j-1]-H1[j] - 3.0*Pr*dh*(F[j]+F[j-1]) *
              (H1[j]+H1[j-1])/4.0;

          if (j == 1) {
              A[j][0][0] = A[j][0][1] = A[j][0][3] = A[j][0][4] = A[j][1][1] =
                A[j][1][2] = A[j][1][4] = A[j][2][0] = A[j][2][2] =
                A[j][2][3] = A[j][3][1] = A[j][3][4] = A[j][4][0] = A[j][4][3] = 0.0;
              A[j][0][2] = 1.0;
              A[j][1][0] = A[j][1][3] = A[j][2][1] = A[j][2][4] = -dh/2.0;
              A[j][3][0] = ksi[j][1];
              A[j][3][2] = ksi[j][2];
              A[j][3][3] = ksi[j][0];
              A[j][4][1] = bet[j][1];
              A[j][4][2] = bet[j][2];
              A[j][4][4] = bet[j][0];
              for (size_t i = 0; i < N_MATRIX; ++i) {
                  for (size_t k = 0; k < N_MATRIX; ++k)
                    alf(i, k) = A[j][i][k];
                }

              alfInv = arma::inv(alf);
              for (size_t i = 0; i < N_MATRIX; ++i) {
                  W[j][i] = 0.0;
                  for (size_t n = 0; n < N_MATRIX; ++n)
                    W[j][i] += (alfInv(i, n) * rArr[j][n]);
                }
            }
          if (j > 1) {
              A[j][0][1] = A[j][0][3] = A[j][0][4] = A[j][1][1] =
                A[j][1][2] = A[j][1][4] = A[j][2][0] = A[j][2][2] =
                A[j][2][3] = A[j][3][4] = A[j][4][3] = 0.0;
              A[j][0][2] = 1.0;
              A[j][1][0] = A[j][2][1] = -1.0;
              A[j][0][0] = A[j][1][3] = A[j][2][4] = -dh/2.0;
              A[j][3][0] = ksi[j][5];
              A[j][3][2] = ksi[j][2];
              A[j][3][3] = ksi[j][0];
              A[j][3][1] = ksi[j][7];
              A[j][4][1] = bet[j][5];
              A[j][4][2] = bet[j][2];
              A[j][4][4] = bet[j][0];
              A[j][4][0] = bet[j][7];

              for (size_t i = 0; i < N_MATRIX; ++i)
                for (size_t k = 0; k < 2; ++k)
                  B[j][i][k] = 0.0;
              B[j][0][3] = B[j][0][4] = B[j][1][2] = B[j][1][4] =
                B[j][2][2] = B[j][2][3] = B[j][3][4] = B[j][4][3] = 0.0;
              B[j][0][2] = -1.0;
              B[j][1][3] = B[j][2][4] = -dh/2.0;
              B[j][3][2] = ksi[j][3];
              B[j][3][3] = ksi[j][1];
              B[j][4][2] = bet[j][3];
              B[j][4][4] = bet[j][1];
              for (size_t i = 0; i < N_MATRIX; ++i)
                for (size_t k = 0; k < N_MATRIX; ++k) {
                    ans[i][k] = 0.0;
                    for (size_t n = 0; n < N_MATRIX; ++n)
                      ans[i][k] += (B[j][i][n] * G[j-1][n][k]);
                  }

              for (size_t i = 0; i < N_MATRIX; ++i)
                for (size_t k = 0; k < N_MATRIX; ++k)
                  alf(i, k) = A[j][i][k]-ans[i][k];
              alfInv = arma::inv(alf);
              for (size_t i = 0; i < N_MATRIX; ++i) {
                  bw[i] = 0.0;
                  for (size_t n = 0; n < N_MATRIX; ++n)
                    bw[i] += (B[j][i][n] * W[j-1][n]);
                  bw[i] =- bw[i]+rArr[j][i];
                }
              for (size_t i = 0; i < N_MATRIX; ++i) {
                  W[j][i] = 0.0;
                  for (size_t n = 0; n < N_MATRIX; ++n)
                    W[j][i] += (alfInv(i, n) * bw[n]);
                }
            }
          if (j < JMax-1) {
              for (size_t i = 0; i < N_MATRIX; ++i)
                for (size_t k = 2; k < N_MATRIX; ++k)
                  C[j][i][k] = 0.0;
              C[j][0][1] = C[j][1][1] = C[j][2][0] = 0.0;
              C[j][0][0] = -dh/2.0;
              C[j][1][0] = C[j][2][1] = 1.0;
              C[j][3][0] = ksi[j][4];
              C[j][3][1] = ksi[j][6];
              C[j][4][0] = bet[j][6];
              C[j][4][1] = bet[j][4];
              for (size_t i = 0; i < N_MATRIX; ++i)
                for (size_t k = 0; k < N_MATRIX; ++k) {
                    G[j][i][k] = 0.0;
                    for (size_t n = 0; n < N_MATRIX; ++n)
                      G[j][i][k] += (alfInv(i, n) * C[j][n][k]);
                  }
            }
        }
      for (size_t j = JMax-1; j > 0; --j) {
          if (j == JMax-1)
            for (size_t i = 0; i < N_MATRIX; ++i)
              D[j][i] = W[j][i];
          if (j < JMax-1) { 
              for (size_t i = 0; i < N_MATRIX; ++i) {
                  bw[i] = 0.0;
                  for (size_t n = 0; n < N_MATRIX; ++n)
                    bw[i] += (G[j][i][n] * D[j+1][n]);
                }
              for (size_t i = 0; i < N_MATRIX; ++i)
                D[j][i] = W[j][i]-bw[i];
            }
        }
    }
  RungeKutta();
  size_t mdsJMax = mdsp_->JMax - 1;

  if (WITH_IMPLICIT_METHODS)
    for (size_t j = ADAMS_METHOD_ORDER; j < mdsJMax; ++j) {
        AdamsIntegrate(adamsBashforthMethod<double>, j);
        AdamsIntegrate(adamsMoultonMethod<double>, j);
      }
  else
    for (size_t j = ADAMS_METHOD_ORDER; j < mdsJMax; ++j)
      AdamsIntegrate(adamsBashforthMethod<double>, j);
  checkResult();
}

//================================
// PolhausenFlow::RungeKutta()
//================================

void PolhausenFlow::RungeKutta() {
  state_type_bf x;
  x[0] = F[0];
  x[1] = F1[0];
  x[2] = F2[0];
  x[3] = H[0];
  x[4] = H1[0];
  boost::numeric::odeint::integrate_n_steps(
      stepper_type_bf(),
      polhausSystem(mdsp_), x, 0.0,
      mdsp_->dh, ADAMS_METHOD_ORDER +1,
      pushresult_observer(mdsp_));
}

//================================
// PolhausenFlow::AdamsIntegrate(...)
//================================

void PolhausenFlow::AdamsIntegrate(
    FlowParameters::adamFunc_bf f, const size_t j) {
  const double dh = mdsp_->dh;
  F2[j+1] = f(dh, &(F3[0]), F2[j], j);
  F1[j+1] = f(dh, &(F2[0]), F1[j], j);
  F[j+1]  = f(dh, &(F1[0]), F[j], j);
  H1[j+1] = f(dh, &(H2[0]), H1[j], j);
  H[j+1]  = f(dh, &(H1[0]), H[j], j);
  F3[j+1] = -3.0*F[j+1]*F2[j+1]+2.0*F1[j+1]*F1[j+1]-H[j+1];
  H2[j+1] = -3.0*mdsp_->Pr*H1[j+1]*F[j+1];
}

//================================
// SiegelFlow ctor
//================================

SiegelFlow::SiegelFlow(
    const double Pr, size_t Jm, const double zerosEdge, double dh)
  : FlowParameters(FlowID::SIEGEL, Pr, Jm, zerosEdge, dh) {
  if (mdsp_ == nullptr)                     // initialization error
    throw;
  try {
      calculate();
      if (WIDE_OUTPUT)
        outBasicflow(mdsp_.get());
    } catch (std::exception &e) {              // calculation error
      std::cerr << e.what() << std::endl;
    }
}

//================================
// SiegelFlow::calculate()
//================================

void SiegelFlow::calculate() {
  double h    = 0,
         dh   = mdsp_->dh,
         Pr   = mdsp_->Pr;
  size_t JMax = mdsp_->JMax;
  if (std::abs(Pr-1.0) > 0.001)
    for (size_t i = 0; i < JMax; ++i) {
        F1[i]  = ((1.0+2.0*h*h/Pr)*erfc(h/sqrt(Pr)) -
             2.0*h*exp(-h*h/Pr)/sqrt(M_PI*Pr) -
             (1.0+2.0*h*h)*erfc(h) +
             2.0*h*exp(-h*h)/sqrt(M_PI)) / (1.0-Pr);
        F2[i] = (4.0*h*erfc(h/sqrt(Pr))/Pr +
             4.0*exp(-h*h)/sqrt(M_PI) -
             4.0*h*erfc(h) -
             4.0*exp(-h*h/Pr)/sqrt(M_PI*Pr)) / (1.0-Pr);
        F3[i] = 4.0*(erfc(h/sqrt(Pr))/Pr - erfc(h)) / (1.0-Pr);
        F[i]  = 0.0;  // 8.0*(exp(-h*h)/sqrt(M_PI) -
                      //  exp(-h*h/Pr)/(sqrt(M_PI)*pow(Pr,1.5)))/
                      //  (1.0-Pr);           //it is F''''
        h += dh;
      }
  else
    for (size_t i = 0; i < JMax; ++i) {
        F1[i] = 2.0*h*(h*erf(h)-h+exp(-h*h)/sqrt(M_PI));
        F2[i] = 4.0*h*(-erfc(h)) + 2.0*exp(-h*h)/sqrt(M_PI);
        F3[i] = 4.0*(-erfc(h)+h*exp(-h*h)/sqrt(M_PI));
        F[i]  =  0.0;  // (12.0*exp(-h*h) -
                       //  8.0*h*h*exp(-h*h))/
                       //  sqrt(M_PI); //it is F''''
        h+=dh;
      }
  h = 0.0;
  for (size_t i = 0; i < JMax; ++i) {
      H[i]  = erfc(h);
      H1[i] = -2.0*exp(-h*h)/sqrt(M_PI);
      H2[i] = 4.0*h*exp(-h*h)/sqrt(M_PI);
      h += dh;
    }
  checkResult();
}

//================================
// PolhausenFlow::pushresult_observer functor
//================================

PolhausenFlow::pushresult_observer::pushresult_observer(
    const std::unique_ptr<MainDataStorage> &mdsp_)
  : mds(mdsp_) {}

void PolhausenFlow::pushresult_observer::operator() (
    const state_type_bf &x, double t) {
  size_t n = static_cast<size_t>(t/mds->dh);
  mds->F[n]  = x[0];
  mds->F1[n] = x[1];
  mds->F2[n] = x[2];
  mds->F3[n] = -3.0*x[0]*x[2] + 2.0*x[1]*x[1] - x[3];
  mds->H[n]  = x[3];
  mds->H1[n] = x[4];
  mds->H2[n] = -3.0*mds->Pr*x[0]*x[4];
}

//================================
// PolhausenFlow::polhausSystem functor
//================================

PolhausenFlow::polhausSystem::polhausSystem(
    const std::unique_ptr<MainDataStorage> &mdsp_)
  : mds(mdsp_) {}

void PolhausenFlow::polhausSystem::operator() (
    const state_type_bf &x, state_type_bf &dxdt, double) {
  dxdt[0] = x[1];
  dxdt[1] = x[2];
  dxdt[2] = -3.0*x[0]*x[2]+2.0*x[1]*x[1]-x[3];
  dxdt[3] = x[4];
  dxdt[4] = -3.0*mds->Pr*x[0]*x[4];
}

//================================
// FlowParameters::useNachts
//================================

void FlowParameters::useNachts(const nachtsInput &inp,
    const std::array<FlowParameters::compd, 4> &baseValues) {
  StabilityPoints calculatedpoints(mdsp_.get());
  nachtsMethod nc(mdsp_.get());
  try {
      nc.NachtsMethodCalculate(inp, baseValues, calculatedpoints);
      calculatedpoints.pushToStorage();
    } catch (std::exception &e) {
      std::cerr << "Nachtsheim method error: " << e.what()  << std::endl;
    }
}

//================================
// FlowParameters::useGaler
//================================

// get approximate values for Nachtsheim method
void FlowParameters::useGaler(
    const double Re, const double wavnum, galerkInput igp,
    std::vector<std::array<FlowParameters::compd, 4>> &baseValues) {
  galerkMethod gl(igp, mdsp_.get());
  try {
      gl.GalerkMethodCalculate(Re, wavnum, baseValues);
    } catch (std::exception &e) {
      std::cerr << "Galerkin method error: " << e.what() << std::endl;
    }
}

const std::vector<double> &FlowParameters::getFlowParameters(size_t i) const {
  if (i > 6) {
      std::cerr << " method getFlowParameters take not correct arg 'i'\n";
      i = 0;
    }
  std::array<std::vector<double> *, 7> vecs =
  {{
    &mdsp_->F,
    &mdsp_->F1,
    &mdsp_->F2,
    &mdsp_->F3,
    &mdsp_->H,
    &mdsp_->H1,
    &mdsp_->H2
  }};
  return *vecs[i];
}

//================================
// FlowParameters ctor
//================================

FlowParameters::FlowParameters(
    FlowID fid, const double Pr, size_t Jm, const double zerosEdge, double dh) {
  try {
      mdsp_ = std::unique_ptr <MainDataStorage>(
          new MainDataStorage(fid, Pr, Jm, zerosEdge, dh));
      F  = &(mdsp_->F [0]);
      F1 = &(mdsp_->F1[0]);
      F2 = &(mdsp_->F2[0]);
      F3 = &(mdsp_->F3[0]);
      H  = &(mdsp_->H [0]);
      H1 = &(mdsp_->H1[0]);
      H2 = &(mdsp_->H2[0]);
    } catch (MainDataStorage::DataStorageExcept &e) {
      std::cerr << e.what() << std::endl;
      throw;
    }
}

//================================
// FlowParameters::checkResult
//================================

void FlowParameters::checkResult() {
  size_t tail = mdsp_->JMax >> 7;  // ~1% JMax
  auto tailaverage = [tail](std::vector<double> &vec) {
      return  std::accumulate(vec.rbegin(), vec.rbegin() + tail, 0.0 ,
                  [](double a, double b) {
                      return a + std::abs(b);
                    })
                  / tail;
    };
  auto abs_max_element =
    [](std::vector<double> &vec) {
        double max_pos = *(std::max_element(vec.begin(), vec.end())),
               max_neg = std::abs(*(std::min_element(vec.begin(), vec.end())));
        return (max_pos > max_neg) ? max_pos : max_neg;
      };
  bool isCorrect = ((abs_max_element(mdsp_->F1)*0.001) >
      tailaverage(mdsp_->F1));
       isCorrect &=  ((abs_max_element(mdsp_->F3)*0.001) >
      tailaverage(mdsp_->F3));
       isCorrect &=  ((abs_max_element(mdsp_->H1)*0.001) >
      tailaverage(mdsp_->H1));
  if (!isCorrect)
    std::cerr << " Please, check result of calculation of a"
              << " basicflow for Pr = " << mdsp_->Pr
              << " there may be a mistake \n";
  else
    std::cerr << " Basic flow parameters valid (doesn't mean correct')\n";
}

//================================
// FlowParameters::getMDS
//================================

const MainDataStorage *FlowParameters::getMDS() const {
  return mdsp_.get();
}
}  // namespace conv_flow
