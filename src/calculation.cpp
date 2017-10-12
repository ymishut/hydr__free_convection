#include "calculation.h"

#include "basicflow.h"
#include "some_math.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <complex>
#include <functional>
#include <future>
#include <vector>
#include <utility>

extern bool WITH_IMPLICIT_METHODS;
extern bool WIDE_OUTPUT;

namespace conv_flow {
//================================
// SelectOptions() ctor
//================================

Calculation::SelectOptions::SelectOptions(
    FlowID fid, const  double Pr, size_t Jm, const double zerosEdge, double dh) {
  try {
      switch (fid) {
        case(FlowID::POLHAUSEN):
          flop = std::unique_ptr<FlowParameters>(
              new PolhausenFlow(Pr, Jm, zerosEdge, dh));
          break;
        case(FlowID::SIEGEL):
          flop = std::unique_ptr<FlowParameters>(
              new SiegelFlow(Pr, Jm, zerosEdge, dh));
          break;
        default:
          std::cerr << "Unknowing type of flow ";
          flop = nullptr;
        }
    } catch (std::exception &e) {                 //  calculate error
      std::cerr << e.what() << std::endl;
      flop = nullptr;
    }
}

//================================
// SelectOptions::useGalerkin(...)
//================================

Calculation::SelectOptions &
Calculation::SelectOptions::useGalerkin(const size_t polynomPower,
                                        const size_t matOrder,
                                        const size_t startOrder) {
  gin = std::unique_ptr<galerkInput>(
      new galerkInput(polynomPower, matOrder, startOrder));
  return *this;
}

Calculation::SelectOptions &
Calculation::SelectOptions::useGalerkin(galerkInput &&gi) {
  gin = std::unique_ptr<galerkInput> (
      new galerkInput(gi));
  return *this;
}

//================================
// SelectOptions::useNachtseim(...)
//================================

Calculation::SelectOptions &
Calculation::SelectOptions::useNachtsheim(const double Re_min,
                                          double Re_max,
                                          const double dRe,
                                          const double w_min,
                                          const double w_max,
                                          const double dw) {
  nin = std::unique_ptr<nachtsInput>(
      new nachtsInput(Re_min, Re_max, dRe, w_min, w_max, dw));
  return *this;
}

//================================
// SelectOptions::build()
//================================

Calculation Calculation::SelectOptions::build() {
  if ((gin == nullptr) && (nin == nullptr))
    throw CalculationExcept("Builder require at least one option use* ");
  return Calculation(flop, std::move(gin), std::move(nin));
}

//================================
// Calculation::calculating(...)
//================================

void Calculation::calculating() {
  if (nin_ == nullptr) {
      std::cerr << " Calculation object has not any"
                << " parameters for calculate ";
      return;
    }
  WIDE_OUTPUT = false;
  std::vector<std::array<compd, 4>> gal_result;

  if (gin_ == nullptr) {
    // one shoot for get start condition
      std::cerr << "Galerkin method use default parameters"
                << " (pow = 2, order = 7, start = 0)\n";
      flop_->useGaler(nin_->Re, nin_->wavnum_min,
          galerkInput(2, 7, 0), gal_result);
    } else {
      flop_->useGaler(nin_->Re, nin_->wavnum_min, *gin_, gal_result);
    }
  if (gal_result.empty()) {
      std::cerr << "Please try call \"calculation\" overload"
                << " with aditional arguments,"
                << " or use Galerkin method for take them \n";
      return;
    }

  if (gin_ == nullptr) {
      concurrNachts(gal_result[0]);
    } else {
      if ((nin_->Re + nin_->dRe) > nin_->Re_max) {
          auto nip = std::unique_ptr<nachtsInput>(
              new nachtsInput(nin_->Re, nin_->Re, 1, nin_->wavnum_min,
                  nin_->wavnum_max, nin_->dwavnum));
          concurrNachts(nip.get(), gal_result);
        } else {
            // call gal_resulr.size() times concurr_nachts(array<compd,4>)
            for_each(gal_result.begin(), gal_result.end(),
                [this](std::array<compd, 4> &stcond) {
                    concurrNachts(stcond);
                  });
          }
    }
}

//================================
// Calculation::calculating(...)
//================================

void Calculation::calculating(const double Re, const double wan,
    const std::array<Calculation::compd, 4> &stcond) {
  bool corr = (Re > 0.0) && (wan > 0.0) &&
    (std::isfinite(Re)) && (std::isfinite(wan));
  for (auto &x : stcond)
    corr &= std::isfinite(std::abs(x));
  if (!corr) {
      std::cerr << "calculating(Re, wan) get not valid parameters";
      return;
    }

  WITH_IMPLICIT_METHODS = true;
  WIDE_OUTPUT = true;

  if (gin_ != nullptr)
    std::cerr << "Calculation object has enough parameters for"
              << " calculation. Extra galerkInput* will be ignored."
              << std::endl;
  if (nin_!= nullptr)
    std::cerr << "Calculation object has enough parameters for"
              << " calculation. Extra nachtsInput* will be ignored."
              << std::endl;
  nin_ = std::unique_ptr<nachtsInput>(
      new nachtsInput(Re, Re, 1, wan, wan, 1));
  flop_->useNachts(*nin_, stcond);
  nin_ = nullptr;
}

//================================
// Calculation::calculating(...)
//================================

void Calculation::calculating(const double Re, const double wan) {
  bool corr = (Re > 0.0) && (wan > 0.0)
    && (std::isfinite(Re)) && (std::isfinite(wan));
  if (!corr) {
      std::cerr << " calculating(Re, wan) get not valid parameters";
      return;
    }

  WITH_IMPLICIT_METHODS = true;
  WIDE_OUTPUT = true;

  std::vector<std::array<compd, 4>> gal_result;

  if (gin_ == nullptr)
    // one shoot for get start condition
    flop_->useGaler(Re, wan, galerkInput(2, 7, 0), gal_result);
  else
    flop_->useGaler(Re, wan, *gin_, gal_result);
  if (gal_result.empty()) {
      std::cerr << "Please try call \"calculation\" overload"
                << " with aditional arguments, or use"
                << " Galerkin method for take them \n";
      return;
    }
  if (nin_ != nullptr)
    std::cerr << "Calculation object has enough parameters for"
              << " calculation. Extra nachtsInput*"
              << "  will be ignored." << std::endl;

  if (gin_ == nullptr) {
      calculating(Re, wan, gal_result[0]);
    } else {
      auto nip = std::unique_ptr<nachtsInput>(
          new nachtsInput(Re, Re, 1, wan, wan, 1));
      concurrNachts(nip.get(), gal_result);
    }
}

//================================
// Calculation() setGalerkinInput
//================================

void Calculation::setGalerkinInput(
    std::unique_ptr<galerkInput> newgin) {
  gin_.swap(newgin);
}

//================================
// Calculation() setNachtsInput
//================================

void Calculation::setNacgtsInput(
    std::unique_ptr<nachtsInput> newnin) {
  nin_.swap(newnin);
}

//================================
// Calculation() copysctors
//================================

Calculation::Calculation(const Calculation &clc) {
  flop_ = clc.flop_;
  gin_  = std::unique_ptr<galerkInput>(new galerkInput(*clc.gin_));
  nin_  = std::unique_ptr<nachtsInput>(new nachtsInput(*clc.nin_));
}

Calculation &Calculation::operator=(const Calculation &clc) {
  if (this != &clc) {
      flop_ = clc.flop_;
      gin_  = std::unique_ptr<galerkInput>(new galerkInput(*clc.gin_));
      nin_  = std::unique_ptr<nachtsInput>(new nachtsInput(*clc.nin_));
    }
  return *this;
}


//================================
// Calculation() ctor
//================================

Calculation::Calculation(std::shared_ptr<FlowParameters> flop,
                         std::unique_ptr<galerkInput> gin,
                         std::unique_ptr<nachtsInput> nin)
  : flop_(flop), gin_(std::move(gin)), nin_(std::move(nin)) {}

//================================
// Calculation::concurrNachts(...)
//================================

// async calls
void Calculation::concurrNachts(
    const std::array<Calculation::compd, 4> &stcond) {
  std::vector<std::future<void>> threads;
  std::vector<nachtsInput> vnip { *nin_ };
  int vnipsize = (nin_->Re_max - nin_->Re) / nin_->dRe + 1;
  if (vnipsize < 1)
    return;
  vnip.reserve(vnipsize);
  threads.reserve(vnipsize);

  double Re = nin_->Re + nin_->dRe;
  while (Re < nin_->Re_max) {
      vnip.push_back(
          { Re, nin_->Re_max, nin_->dRe, nin_->wavnum_min,
              nin_->wavnum_max, nin_->dwavnum
          });
      Re += nin_->Re + nin_->dRe;
    }

  for_each(vnip.begin(), vnip.end(),
      [this, &stcond, &threads](const nachtsInput &in) {
        threads.push_back(
            std::async(&FlowParameters::useNachts, flop_.get(),
                in, std::cref(stcond)));
        });

  for_each(threads.begin(), threads.end(),
      [](std::future<void> &t) { t.get(); });
}

void Calculation::concurrNachts(nachtsInput *nip,
    std::vector<std::array<Calculation::compd, 4>> &stcond) {
  std::vector<std::future<void>> threads;
  for (size_t i = 0; ; ++i) {
      if (stcond.empty())
        break;
      threads.push_back(std::async(&FlowParameters::useNachts, flop_.get(),
          std::cref(*nip), std::ref(*(stcond.rbegin() )) ));
      stcond.pop_back();
    }
  for_each(threads.begin(), threads.end(),
      [](std::future<void> &t) { t.get(); });
}
}  // namespace conv_flow
