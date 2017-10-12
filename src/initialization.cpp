#include "initialization.h"

#include "inputdata_by_file.h"
#include "maindatastorage.h"

namespace conv_flow {
//================================
// InitializationByFile::funclist  map
//================================

std::map<std::string, InitializationByFile::pf>
  InitializationByFile::funclist =
  {
    {"NACHTS",   &InitializationByFile::setNachts},
    {"GALERKIN", &InitializationByFile::setGalerkin},
    {"BOTH",     &InitializationByFile::setBoth}
  };

//================================
// InitializationByFile ctor
//================================

InitializationByFile::InitializationByFile(
    const std::__cxx11::string &file) {
  reset(file);
}

//================================
// Initialization reset
//================================

void InitializationByFile::reset(const std::__cxx11::string &file) {
  try {
      inp = std::unique_ptr<InputData>(new InputData(file));
      setclc();
    } catch (InputData::InputExcept &e) {
      std::cerr << e.what() << std::endl;
      inp = nullptr;
      clc = nullptr;
      return;
    }
}

//================================
// Initialization setclc
//================================

void InitializationByFile::setclc() {
  try {
      pf method = InitializationByFile::funclist[inp->getCalcMethods()];
      setbldr();
      (this->*method)();
    } catch (InputData::InputExcept &e) {
      std::cerr << e.what() << std::endl;
    }
  return;
}

//================================
// Initialization setbldr
//================================
// set builder
void InitializationByFile::setbldr() {
  auto pr  = inp->getFlow();
  auto arr = inp->getCalcParams();
  bldr = std::unique_ptr<Calculation::SelectOptions>(
      new Calculation::SelectOptions(
          pr.first, pr.second, arr[1], arr[0], arr[2]));
}

//================================
// Initialization setNachts
//================================

void InitializationByFile::setNachts() {
  auto Re  = inp->getReBoundaries();
  auto Wav = inp->getWavnumBoundaries();
  clc = std::unique_ptr<Calculation>(
      new Calculation(bldr->useNachtsheim(
          Re[0], Re[1], Re[2],
          Wav[0], Wav[1], Wav[2]).build()));
}

//================================
// Initialization setGalerkin
//================================

void InitializationByFile::setGalerkin() {
  clc = std::unique_ptr<Calculation>(
      new Calculation(bldr->useGalerkin(inp->getGalerkParams()).build()));
}

//================================
// Initialization setBoth
//================================

void InitializationByFile::setBoth() {
  auto Re  = inp->getReBoundaries();
  auto Wav = inp->getWavnumBoundaries();
  clc = std::unique_ptr<Calculation>(
      new Calculation(bldr->useNachtsheim(
          Re[0], Re[1], Re[2], Wav[0], Wav[1], Wav[2]).
              useGalerkin(inp->getGalerkParams()).build()));
}

//================================
// Initialization getCalculator
//================================

Calculation *Initialization::getCalculator() {
   return (clc == nullptr) ? nullptr : clc.get();
}
}  // namespace conv_flow
