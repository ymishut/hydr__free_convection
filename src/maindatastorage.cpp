#include "maindatastorage.h"

#include "some_math.h"

#include <algorithm>
#include <functional>
#include <fstream>
#include <tuple>
#include <vector>

// with Adams - Moulton method
bool WITH_IMPLICIT_METHODS = true;
// write all calculated parameters in files
bool WIDE_OUTPUT           = true;
// for change this - realize Adams - Bashforth - Moulton 
// methods for another order
int  ADAMS_METHOD_ORDER    = 7;
// calculation accuracy, not machinery, but if you use another float -
// reduce or increase this
double ACCURACY            = 10e-6;
// number of points in output files (for basic flow and "perturbation" flow)
size_t POINT_COUNT         = 100;

namespace conv_flow {
//================================
// MainDataStorage ctor
//================================

MainDataStorage::MainDataStorage(
    FlowID fid, const double _Pr, const size_t _JMax,
    const double _zerosEdge, const double _dh)
  : flowID(fid), JMax(_JMax), Pr(_Pr), zerosEdge(_zerosEdge), dh(_dh) {
  if (zerosEdge > dh*JMax)
    JMax = zerosEdge/dh + 1;
  auto isOK = [](double val) {
      return (val > 0.0 || std::isfinite(val));
    };
  bool correctInput =
      isOK(Pr) && isOK(zerosEdge) && isOK(zerosEdge) && isOK(dh) &&
          (JMax < 100000) && (JMax > 500);
  if (!correctInput)
    throw DataStorageExcept("Bad input for MainDataStorage");

  // legacy code for Polhauzen flow
  //   expect it are pointers
  F.assign(JMax, 0.0);  F1.assign(JMax, 0.0);
  F2.assign(JMax, 0.0); F3.assign(JMax, 0.0);
  H.assign(JMax, 0.0);  H1.assign(JMax, 0.0); H2.assign(JMax, 0.0);
}

//================================
// MainDataStorage getMutex
//================================

std::mutex &MainDataStorage::getMutex() {
  return mtx;
}

//================================
// MainDataStorage dtor
//================================

MainDataStorage::~MainDataStorage() {
  try {
    std::ofstream fout("log_points.txt");
    for_each(points.begin(), points.end(),
        [this, &fout](const stPoint &p) { fout << p; });
    } catch (std::exception &e) {}
}

//================================
// MainDataStorage::DataStorageException class
//================================

MainDataStorage::DataStorageExcept::DataStorageExcept(
    const std::__cxx11::string &message)
  : message_(message) {}

MainDataStorage::DataStorageExcept::DataStorageExcept(
    std::__cxx11::string &&message) noexcept
  : message_(message) {}

const char *MainDataStorage::DataStorageExcept::what() const noexcept {
  return message_.c_str();
}

//================================
// stPoint ctor
//================================

stPoint::stPoint(double R, double wn, compd cb, compd f3, compd s1)
  : Re(R), wavnum(wn), phvel(cb), f3b(f3), s1b(s1) {}

//================================
// stPoint operator<
//================================

bool operator<(const stPoint &lsp, const stPoint &rsp) {
  double phr  = rsp.phvel.real(),
         phl  = lsp.phvel.real();
  return (std::tie(lsp.Re, lsp.wavnum, phl) <
          std::tuple<double, double, double>(rsp.Re + ACCURACY,
            rsp.wavnum + ACCURACY, rsp.phvel.real() + ACCURACY) &&
          std::tie(rsp.Re, rsp.wavnum, phr) >
          std::tuple<double, double, double>(lsp.Re + ACCURACY,
            lsp.wavnum + ACCURACY, lsp.phvel.real() + ACCURACY));
}

//================================
// stPoint operator<<
//================================

std::ostream &operator<<(std::ostream &out, const stPoint &rsp) {
  out << rsp.Re << "  \t" << rsp.wavnum << "  \t" << rsp.f3b
      << "  \t" << rsp.s1b << "  \t" << rsp.phvel << std::endl;
  return out;
}
}  // namespace conv_flow
