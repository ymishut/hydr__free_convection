#ifndef SRC_MAINDATASTORAGE_H_
#define SRC_MAINDATASTORAGE_H_

#include <boost/noncopyable.hpp>

#include <cstddef>

#include <complex>
#include <iostream>
#include <mutex>
#include <set>
#include <string>
#include <vector>

namespace conv_flow {
//================================
// enum FlowID
//================================

enum class FlowID {POLHAUSEN = 0, SIEGEL};

//================================
// stPoint struct
//================================

struct stPoint {
 private:
  typedef std::complex<double> compd;
 public:
  double Re,
         wavnum;
  compd  phvel,
         f3b,
         s1b;
  stPoint(double R, double wn, compd cb, compd f3, compd s1);
};

bool operator<(const stPoint &lsp, const stPoint &rsp);

std::ostream &operator<<(std::ostream& out, const stPoint &rsp);

//================================
// MainDataStorage
//================================

class MainDataStorage : private boost::noncopyable {
  typedef std::complex<double> compd;

 public:
  MainDataStorage(FlowID fid, const double Pr, const size_t JMax, 
                  const double zerosEdge, const double dh);
  std::mutex &getMutex();
  ~MainDataStorage();

 public:
  std::vector<double> F, F1, F2, F3, H, H1, H2;
  std::set<stPoint> points;
  const FlowID flowID;
  size_t JMax;
  double Pr,
         zerosEdge,
         dh;

 public:
  class DataStorageExcept;

 private:
  std::mutex mtx;
};

//================================
// DataStorageException
//================================

class MainDataStorage::DataStorageExcept : public std::exception {
 public:
  DataStorageExcept(const std::string &message);
  DataStorageExcept(std::string &&message) noexcept;
  const char *what() const noexcept;
  ~DataStorageExcept() noexcept {}

 private:
  std::string message_;
};
}  // namespace conv_flow
#endif  // SRC_MAINDATASTORAGE_H_
