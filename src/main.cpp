#include "animation.h"
#include "basicflow.h"
#include "calculation.h"
#include "initialization.h"
#include "maindatastorage.h"

#include <string>

#define ByFileCheck
#define AnotherCheck
#define JFF

int main(int argc, char **argv) {
#ifdef ByFileCheck
  std::string file = "input.txt";
  conv_flow::InitializationByFile a(file);
  conv_flow::Calculation *pc = a.getCalculator();
  if (pc != nullptr)
    pc->calculating();
#endif  // ByFileCheck

  // GALERKIN METHOD CAN GET ONLY 0 and 2 POLYNOM POWER 
  // (it's specifical for boundary conditions)

#ifdef AnotherCheck
  double const Pr = 6.7;
  double ZEROS_EDGE = 5.2;
  size_t JMax = 6001;
  double dh = 0.0025;
  conv_flow::Calculation::SelectOptions fg(
      conv_flow::FlowID::POLHAUSEN, Pr, JMax, ZEROS_EDGE, dh);
  conv_flow::Calculation cl = fg.useGalerkin(0, 7, 2).
      useNachtsheim(34, 35, 5, 0.45, 0.46, 0.2).build();
  cl.calculating();
#endif  // AnotherCheck

#ifdef JFF
  conv_flow::animation::Drawing::StartDrawing(argc, argv);
#endif  // JFF

  return 0;
}
