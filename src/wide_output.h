#ifndef SRC_WIDE_OUTPUT_H_
#define SRC_WIDE_OUTPUT_H_

#include <complex>
#include <iomanip>
#include <ostream>
#include <vector>

//================================
// formatted_output
//================================

class formatted_output {
 private:
  int width     = 11,
      precision = 8;
  std::ostream &stream_obj;

 public:
  explicit formatted_output(std::ostream &obj)
    : stream_obj(obj) {}

  template<class T>
  formatted_output &operator<<(const T &output) {
    stream_obj << std::setprecision(precision)
               << std::setw(width) << output << " ";
    return *this;
  }

  inline formatted_output &operator<<(std::ostream &(*func)(std::ostream &)) {
    func(stream_obj);
    return *this;
  }
};

namespace conv_flow {
  class MainDataStorage;
  void outWaves(const std::vector<std::complex<double>> &fi,
                const std::vector<std::complex<double>> &fi1,
                const std::vector<std::complex<double>> &si,
                                           const double dh);
  void outBasicflow(const MainDataStorage *mds);
  void outGalerkin(const std::vector<std::complex<double>> &eigenValues,
                   const double Re, const double wavnum, const double expc);
  void outStPoints(const MainDataStorage *mds);
}  // namespace conv_flow
#endif  // SRC_WIDE_OUTPUT_H_
