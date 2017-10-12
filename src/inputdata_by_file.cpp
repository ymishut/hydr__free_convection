#include "inputdata_by_file.h"

#include "galerk_method.h"
#include "filereading.h"
#include "nachts_method.h"

#include <boost/lexical_cast.hpp>

#include <exception>
#include <fstream>
#include <iostream>

namespace conv_flow {
//================================
// InputData ctor
//================================

InputData::InputData(const std::__cxx11::string &filename) {
  std::ifstream infile_(filename, std::ios_base::in);
  if (!infile_.is_open()) {
      throw InputExcept("Cannot open file to read init data: "
          + filename + "\n put its in same directory with .exe or .out ");
    }
  try {
      ReadFile d;
      inputParameters_ = (d.parseFile(infile_));
    } catch (ReadFile::ReadFileExcept &e) {
      std::cerr << e.what() << std::endl;
      throw;
    }
  if (inputParameters_.size() != 13)
    throw InputExcept(" Count of input parameters is incorrect. Must be 13\n");
}

//================================
// InputData getCalcMethods
//================================

std::__cxx11::string InputData::getCalcMethods() const {
  return inputParameters_.at(0);
}

//================================
// InputData getFlow
//================================

std::pair<FlowID, double> InputData::getFlow() const {
  double pr_num;
  FlowID fid;
  try {
      std::string flowidstr = inputParameters_.at(1);
      fid = (flowidstr == "POLHAUSEN") ? FlowID::POLHAUSEN :
          (flowidstr == "SIEGEL") ? FlowID::SIEGEL :
              throw InputExcept("unknowing flow type in input file");
      pr_num = boost::lexical_cast<double>(inputParameters_.at(2));
    } catch (std::exception &e) {
      std::cerr << e.what() << std::endl;
      throw;
    }
  return {fid, pr_num};
}

//================================
// InputData castArray
//================================

std::array<double, 3> InputData::castArray(size_t i, std::string *str) const {
  std::array<double, 3> arr;
  try {
      arr[0] = boost::lexical_cast<double>(str[i]);
      arr[1] = boost::lexical_cast<double>(str[i+1]);
      arr[2] = boost::lexical_cast<double>(str[i+2]);
    } catch (std::exception &e) {
      std::cerr << e.what() << std::endl;
    }
  return arr;
}

//================================
// InputData getCalcParams
//================================

std::array<double, 3> InputData::getCalcParams() {
  return castArray(3, &inputParameters_[0]);
}

//================================
// InputData getGalerkParams
//================================

galerkInput InputData::getGalerkParams() {
  std::array<double, 3> arr = castArray(6, &inputParameters_[0]);
  return galerkInput(arr[0], arr[1], arr[2]);
}

//================================
// InputData castBoundary
//================================

std::array<double, 3> InputData::castBoundary(size_t i) const {
  size_t j = inputParameters_[i].find(',');
  std::string str[] =
  {
    inputParameters_[i].substr(0, j),
    inputParameters_[i].substr(j+1),
    inputParameters_[i+2]
  };
  return castArray(0, str);
}

//================================
// InputData getReBoundaries
//================================

std::array<double, 3> InputData::getReBoundaries() const {
  return castBoundary(9);
}

//================================
// InputData getWavnumBoundaries
//================================

std::array<double, 3> InputData::getWavnumBoundaries() const {
  return castBoundary(10);
}
}  // namespace conv_flow

