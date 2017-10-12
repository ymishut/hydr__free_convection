#ifndef SRC_INITIALIZATION_H_
#define SRC_INITIALIZATION_H_

#include "calculation.h"
#include "inputdata_by_file.h"

#include <map>
#include <memory>
#include <string>

namespace conv_flow {
//================================
// Initialization
//================================

class Initialization {
 public:
  Calculation *getCalculator();

 protected:
  std::unique_ptr<Calculation> clc;
  std::unique_ptr<Calculation::SelectOptions> bldr;

 protected:
  virtual void setNachts()   = 0;
  virtual void setGalerkin() = 0;
  virtual void setBoth()     = 0;

  virtual ~Initialization() {}
};

//================================
// InitializationByFile
//================================

class InitializationByFile final : public Initialization {
  typedef void (InitializationByFile::*pf)();

 public:
  explicit InitializationByFile(const std::string &file);
  void reset(const std::string &file);

 private:
  std::unique_ptr<InputData> inp;
  static std::map<std::string, pf> funclist;

 private:
  void setclc();
  void setbldr();
  void setNachts() override;
  void setGalerkin() override;
  void setBoth() override;
};
}  // namespace conv_flow
#endif  // SRC_INITIALIZATION_H_
