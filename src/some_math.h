#ifndef SOME_MATH_H_
#define SOME_MATH_H_

#include <cstddef>

#include <complex>
#include <cmath>
#include <type_traits>

template<class T>
inline T adamsBashforthMethod(const double dh, const T *arr, 
                              const T prev, const size_t cur) {
  return prev +
    dh * arr[cur] * 434241.0 / 120960.0 - dh * arr[cur-1] * 1152169.0 / 120960.0 +
    dh * arr[cur-2] * 2183877.0 / 120960.0 - dh * arr[cur-3] * 2664477.0 / 120960.0 +
    dh * arr[cur-4] * 2102243.0 / 120960.0 - dh * arr[cur-5] * 1041723.0 / 120960.0 +
    dh * arr[cur-6] * 295767.0 / 120960.0 - dh * arr[cur-7] * 36799.0 / 120960.0;
}

template<class T>
inline T adamsMoultonMethod(const double dh, const T *arr,
                            const T prev, const size_t cur) {
  return prev +
    dh * arr[cur+1] * 36799.0 / 120960.0 + dh * arr[cur] * 139849.0 / 120960.0 -
    dh * arr[cur-1] * 121797.0 / 120960.0 + dh * arr[cur-2] * 123133.0 / 120960.0 -
    dh * arr[cur-3] * 88547.0 / 120960.0 + dh * arr[cur-4] * 41499.0 / 120960.0 -
    dh * arr[cur-5] * 11351.0 / 120960.0 + dh * arr[cur-6] * 1375.0 / 120960.0;
}

template<class T, class Acc_T,
      class = typename std::enable_if<std::is_floating_point <T>::value>::type>
bool cmplx_eq(const std::complex<T> &l, const std::complex<T> &r,
              Acc_T accuracy = 0.0001) {
  bool    isEq = (((l.real() + accuracy) > r.real()) ==
       ((r.real() + accuracy) > l.real()));
  return isEq && (((l.imag() + accuracy) > r.imag()) ==
      ((r.imag() + accuracy) > l.imag()));
}

template<class T>
void setTriangular(T **Mat, const size_t nrows, const size_t ncols) {
  for (size_t i = 0; i < nrows; ++i)
    if (std::abs(Mat[i][i]) < 0.001)
      for (size_t z = 0; z < nrows; ++z) {
          if (z == i)
            continue;
          for (size_t k = 0; k < ncols; ++k)
            Mat[i][k] += Mat[z][k];
          if (std::abs(Mat[i][i]) > 0.001)
            break;
        }
  for (size_t i = 0; i < (nrows-1); ++i)
    for (size_t j = i+1; j < nrows; ++j) {
        T temp = Mat[j][i]/Mat[i][i];
        for (size_t k = i; k < ncols; ++k)
          Mat[j][k] -= Mat[i][k]*temp;
      }
}

template<class T>
void gaussMethod(T **Mat, T *ans, const size_t nrows, bool &isSingular) {
  setTriangular(Mat, nrows, nrows+1);
  size_t ndiag_zeroes = 0;
  size_t i = nrows-1;
  while (i < nrows) {
      if (std::abs(Mat[i][i]) < 0.001) {
          ++ndiag_zeroes;
          ans[i] = 1.0;
        } else {
          ans[i] = Mat[i][nrows];
          for (size_t j = i+1; j < nrows; j++)
            ans[i] -= Mat[i][j]*ans[j];
          ans[i] /= Mat[i][i];
        }
      --i;
    }
  if (ndiag_zeroes != 0)
    isSingular = true;
}
#endif  // SOME_MATH_H_
