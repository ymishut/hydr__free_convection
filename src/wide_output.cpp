#include "wide_output.h"

#include "maindatastorage.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>

extern size_t POINT_COUNT;
typedef std::complex<double> compd;

namespace conv_flow {
//================================
// outWaves
//================================

void outWaves(const std::vector<compd> &fi,
              const std::vector<compd> &fi1,
              const std::vector<compd> &si, const double dh) {
  std::string file("parameters/nachts_results.txt");
  std::ofstream of(file, std::ios_base::out);
  if (!of.is_open()) {
      std::cerr << "Cannot open output file: " << file << std::endl;
      return;
    }
  formatted_output fout(of);
  size_t leng  = size_t{fi.size()/POINT_COUNT},
         temp  = 0;
  double x = 0.0;
  fout << "#x          |fr         |fi         |fr'" 
       << "        |fi'        |sr         |si\n";
  for (size_t i = 0; i < POINT_COUNT; ++i) {
      x = temp * dh;
      fout << x << fi[temp].real() << fi[temp].imag() << fi1[temp].real()
           << fi1[temp].imag() << si[temp].real() << si[temp].imag() << "\n";
      temp += leng;
    }
}

//================================
// outBasicflow
//================================

void outBasicflow(const MainDataStorage *mds) {
  std::string file("parameters/basicflow.txt");
  std::ofstream of(file, std::ios_base::out);
  if (!of.is_open()) {
      std::cerr << "Cannot open output file: " << file << std::endl;
      return;
    }
  formatted_output fout(of);
  size_t leng  = size_t{mds->JMax/POINT_COUNT},
         temp  = 0;
  double x = 0.0;
  fout << "#x          |F''          |F'''        |H'" << std::endl;
  for (size_t i = 0; i < POINT_COUNT; ++i) {
      x = temp * mds->dh;
      fout << x << mds->F1[temp] << mds->F3[temp] << mds->H1[temp] << "\n";
      temp += leng;
    }
}

//================================
// outGalerkin
//================================

void outGalerkin(const std::vector<compd> &eigenValues,
    const double Re, const double wavnum, const double expc) {
  std::string file("parameters/galerkin_results.txt");
  std::ofstream of(file, std::ios_base::app);
  if (!of.is_open()) {
      std::cerr << "Cannot open output file: " << file << std::endl;
      return;
    }
  formatted_output fout(of);
  fout << "#Re _1, wavenumber _2 and c _3 in exp(-cx) :" << Re
       << wavnum << expc << "\n";
  for_each(eigenValues.begin(), eigenValues.end(),
      [&fout](const compd &c) {
          fout << c.real() << c.imag() << "\n";
        });
}

//================================
// outStPoints
//================================

void outStPoints(const MainDataStorage *mds) {
  std::string file("parameters/points.txt", std::ios_base::app);
  std::ofstream of(file);
  if (!of.is_open()) {
      std::cerr << "Cannot open output file: " << file << std::endl;
      return;
    }
  formatted_output fout(of);
  fout << "#Re        |wav        |f3       |s1        |ph  \n";
  for_each(mds->points.begin(), mds->points.end(),
      [&fout](const stPoint &p) {
          fout << p.Re << p.wavnum << p.f3b << p.s1b << p.phvel;
        });
}
}  // namespace conv_flow
