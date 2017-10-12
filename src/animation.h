#ifndef SRC_ANIMATION_H_
#define SRC_ANIMATION_H_

#include <boost/noncopyable.hpp>
#include <GL/glut.h>

#include <complex>
#include <memory>
#include <vector>

namespace conv_flow {
class FlowParameters;

namespace animation {
class Drawing;
//================================
// point struct
//================================

struct point {
GLfloat y,
        t;    // temperature
size_t  x_ind;
};

//================================
// scales struct
//================================

struct scales {
std::vector<GLfloat> Velocity,
                     Length,
                     yZerosEdge,
                     PolhausenU,
                     PolhausenV;
};

//================================
// perturb struct
//================================

struct perturb {
 private:
  typedef std::vector<std::complex<GLfloat>> veccompc;

 public:
  veccompc fi,
           fi1,
           si;
};

//================================
// basicflow struct
//================================

struct basicflow {
  std::vector<GLfloat> F,
                       F1,
                       H;
};

//================================
// Animation
//================================

class Animation : private boost::noncopyable {
  typedef std::complex<GLfloat> compf;
  friend class Drawing;

 private:
  std::unique_ptr<FlowParameters> flop_;
  std::unique_ptr<scales>    scales_;
  std::unique_ptr<perturb>   perturb_;
  std::unique_ptr<basicflow> basicflow_;
  std::vector<point> points_;
  std::vector<GLfloat> xCoords_;
  std::vector<compf> exp_ksi_;
  static const size_t npoints_  = 1000,
                      pointsinrow_ = 2;
  size_t  nlines_;
  GLfloat xbottom_,
          xtop_,
          dx_,
          globalscale_,
          globaltime_;
  compf   phvel_;

 private:
  Animation();
  void calcBasicflow();
  void calcPerturb();
  void setExp_ksi();
  void setScales();
  void setPoints();
  void updatePoints();
};

//================================
// Drawing
//================================

class Drawing {
 public:
  static void StartDrawing(int argc, char**argv);

 private:
  static std::unique_ptr<Animation> anim;

 private:
  Drawing();
  static void timer(int = 0);
  static void drawPoints();
  static void draw();
  static void initialization();
  static void ReshapeFunction(GLsizei w, GLsizei h);
};
}  // namespace animation
}  // namespace conv_flow
#endif  // SRC_ANIMATION_H_
