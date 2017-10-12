#include "animation.h"

#include "basicflow.h"
#include "nachts_method.h"
#include "stabilitypoints.h"

#include <algorithm>
#include <functional>
#include <iterator>

#define imaginaryOne std::complex<GLfloat>(0.0, 1.0)

namespace conv_flow {
std::unique_ptr<animation::Animation>
    animation::Drawing::anim = nullptr;

namespace {
  // flow parameters
  const GLfloat Pr             = 6.7,
                dh             = 0.0025,
                wavnum         = 0.45,
                KIN_VISC       = 0.0001,
              //  COEF_THERM_EXP = 0.000207,
              //  g              = 9.81,
              //  gb             = g*COEF_THERM_EXP,    // 0.00203067
                gb             = 0.00203,
                ZEROS_EDGE     = 5.5;

  const int      Re             = 40,
              //  T_WALL         = 100,
              //  T_INF          = 20,
              //  dT             = T_WALL-T_INF,        // 80
                dT             = 80,
                JMax           = 2200;
  // animation parameters
  const GLfloat GLLength       = 10.0,
                dtau           = 40,
  // 10 ^ (-3) * 0.05  //^-3- millisec to sec and 0.5*10 ^-1 - time delay
                timeScale      = 0.00005,
                yBegin         = 5.0;
}  // namespace

//================================
// Animation ctor
//================================

animation::Animation::Animation() {
  basicflow_ = std::unique_ptr<basicflow>(new basicflow());
  perturb_   = std::unique_ptr<perturb>(new perturb());
  scales_    = std::unique_ptr<scales>(new scales());
}

//================================
// calcBasicflow
//================================

void animation::Animation::calcBasicflow() {
  flop_ = std::unique_ptr<FlowParameters>(
      new PolhausenFlow(Pr, JMax, ZEROS_EDGE, dh));
}

//================================
// calcPerturb
//================================

void animation::Animation::calcPerturb() {
  typedef std::complex<double> compd;
  typedef std::vector<compf> vcf;
  typedef std::vector<compd> vcd;
  nachtsInput nin(Re, Re, 1, wavnum, wavnum, 1);
  std::array<compd, 4> stcond =
  {{
    compd{1.0, 0.0},
    compd{0.447, -1.33},
    compd{4.526, -3.167},
    compd{0.20, 0.0}
  }};
  nachtsMethod nm(flop_->getMDS(), true);
  StabilityPoints sp(const_cast<MainDataStorage *>(flop_->getMDS()));
  nm.NachtsMethodCalculate(nin, stcond, sp);
  phvel_ = nm.getPhasevelocity();

  auto copy_fun = [](const vcd &inp, vcf &out) {
    out.reserve(inp.size());
    for (size_t i = 0; i < inp.size(); ++i)
      out.push_back(compf(inp[i].real(), inp[i].imag()));
  };
  copy_fun(nm.getFlowparameters(0), perturb_->fi);
  copy_fun(nm.getFlowparameters(1), perturb_->fi1);
  copy_fun(nm.getFlowparameters(5), perturb_->si);
  basicflow_->F.reserve(JMax);
  basicflow_->F1.reserve(JMax);
  basicflow_->H.reserve(JMax);
  std::copy_n(flop_->getFlowParameters(0).cbegin(), JMax,
      std::back_inserter(basicflow_->F));
  std::copy_n(flop_->getFlowParameters(1).cbegin(), JMax,
      std::back_inserter(basicflow_->F1));
  std::copy_n(flop_->getFlowParameters(4).cbegin(), JMax,
      std::back_inserter(basicflow_->H));
}

//================================
// setScales
//================================

void animation::Animation::setScales() {
  nlines_ = size_t{npoints_ / pointsinrow_};
  size_t halfLines = nlines_ >> 1;
  // calculate parameters with dimension [meter]
  // 64.0 and 0.3333 Polhauzen flow constants
  GLfloat xmiddle  = std::pow(Re * KIN_VISC * KIN_VISC /
              (64.0f * gb * dT), 0.33333f);
          dx_      = xmiddle / (nlines_ - 1);
          xtop_    = xmiddle + (halfLines - 1) * dx_;
          xbottom_ = xmiddle - halfLines * dx_;

  dx_ = (xtop_-xbottom_) / nlines_;

  // transition to dimensionless parameters
  globaltime_  = 0.0f;
  globalscale_ = GLLength / (xtop_-xbottom_);
  xbottom_     = 5.0f;
  dx_          = GLLength / nlines_;
  xtop_        = (GLLength + 5.0f) - dx_;
  xCoords_.push_back(xbottom_);
  std::generate_n(std::back_inserter(xCoords_), nlines_ - 1,
      [this]() {
          return xCoords_.back() + dx_;
        });

  GLfloat SQRTGRx, 
          // x with dimension (meter)
          x = xbottom_/globalscale_,       
          dxnotglobal = dx_ / globalscale_;
  scales_->PolhausenU.assign(nlines_, 0.0f);
  scales_->PolhausenV.assign(nlines_, 0.0f);
  scales_->Length.assign(nlines_, 0.0f);
  scales_->yZerosEdge.assign(nlines_, 0.0f);
  scales_->Velocity.assign(nlines_, 0.0f);
  const GLfloat SQRT2 = std::sqrt(2.0f);
  for (size_t i = 0; i < nlines_; ++i) {
      SQRTGRx = std::sqrt(gb * dT * x*x*x / (KIN_VISC*KIN_VISC));
      scales_->PolhausenU[i] = globalscale_ *
          SQRTGRx *8.0f*KIN_VISC * std::sqrt(Pr) / x;
      scales_->PolhausenV[i] = globalscale_ * std::sqrt(SQRTGRx) *
          SQRT2 * KIN_VISC * std::pow(Pr, 0.75) / x;
      scales_->Length[i]     = SQRT2 * x / std::sqrt(SQRTGRx);
      scales_->yZerosEdge[i] = globalscale_ * ZEROS_EDGE *
          scales_->Length.back();
      scales_->Velocity[i]   = globalscale_ *2.0f*KIN_VISC * SQRTGRx / x;
      x += dxnotglobal;
    }
}

//================================
// setPoints
//================================

void animation::Animation::setPoints() {
  std::srand(unsigned(std::time(0)));
  size_t j = 0,
         x_ind;
  GLfloat y;
  for (size_t i = 0; i < nlines_; ++i)
    for (size_t k = 0; k < pointsinrow_; ++k) {
      x_ind = i;
      if (scales_->yZerosEdge[i] > yBegin)
        y = static_cast<GLfloat>(std::rand() & 1023)
            * 0.00096f * yBegin;
      else
        y = static_cast<GLfloat>(std::rand() & 1023)
            * 0.00096f * scales_->yZerosEdge[i];
      ++j;
      points_.push_back({y, 0.0f, x_ind});
    }
  // this variables manage points color
  // 0.7 is basic red color intensive for "clear" flow -
  // without perturbasions.
  GLfloat temperConst = 0.7f,
          temper1     = 1.0f - temperConst;
  auto norm = 1.0f / std::abs(
      *std::max_element(perturb_->si.begin(), perturb_->si.end(),
          [](const compf a, const compf b) {
              return (std::abs(a) - std::abs(b)) < -0.0001;
            }));
  std::transform(perturb_->si.begin(), perturb_->si.end(),
                 perturb_->si.begin(), std::bind2nd(
      std::multiplies<std::complex<GLfloat>>(), temper1 * norm));
  std::transform(basicflow_->H.begin(), basicflow_->H.end(),
                 basicflow_->H.begin(), std::bind2nd(
      std::multiplies<GLfloat>(), temperConst));
  exp_ksi_.assign(nlines_, 0.0);
}

//================================
// setExp_ksi
//================================

void animation::Animation::setExp_ksi() {
  for (size_t i = 0; i < nlines_; ++i)
    exp_ksi_[i] = std::exp(imaginaryOne * wavnum *
        (xCoords_[i] - phvel_ * scales_->Velocity[i] *
        globaltime_)) / scales_->Length[i] / globalscale_;
}

//================================
// updatePoints
//================================

void animation::Animation::updatePoints() {
  setExp_ksi();
  size_t temphight  = 0,
         halfLength = nlines_ >> 1;
// u is vertical velocity
// v is perpendicular to a wall velocity
// h is length (hight)
  for (size_t i = 0; i < npoints_; ++i) {
      GLfloat v,
              h = points_[i].y;
      size_t nhw = static_cast<size_t>(h / dh),
             nh  = nhw >> 1;
      if (nh > JMax)
        nh = 1;
      GLfloat u = scales_->Velocity[points_[i].x_ind] *
                  std::real(perturb_->fi1[nh] * exp_ksi_[points_[i].x_ind]) +
                  scales_->PolhausenU[points_[i].x_ind] * basicflow_->F1[nhw],
              x = xCoords_[points_[i].x_ind] + u * dtau * timeScale;
      if (x <= xtop_) {
          v = scales_->Velocity[points_[i].x_ind] * std::real(-imaginaryOne *
              wavnum * perturb_->fi[nh]*exp_ksi_[points_[i].x_ind]) +
              scales_->PolhausenV[points_[i].x_ind] *
              (h * basicflow_->F1[nhw] - 3.0f*basicflow_->F[nhw]);
          points_[i].y += v * dtau * timeScale;
          // if nearly to a wall
          if (points_[i].y <= 0.02f) {
              points_[i].y = yBegin;
              points_[i].x_ind = nlines_ >> 2;
              points_[i].t = 0.0f;
              continue;
          }
          size_t ind = static_cast<size_t>((x - xbottom_) / dx_);
          points_[i].x_ind = ind % nlines_;
        } else {
          // set few points in certain places
          // it will look like lines of stream
          if (temphight < halfLength) {
              points_[i].y = yBegin;
              points_[i].x_ind = temphight;
              temphight += 20;
            } else {
              points_[i].x_ind = 0;
              points_[i].y = static_cast<GLfloat>(std::rand() & 1023) *
                  0.00096f * yBegin;
            }
        }
      // temperature function
      h   = points_[i].y;
      nhw = static_cast<size_t>(h / dh);
      nh  = nhw >> 1;
      points_[i].t = basicflow_->H[nhw] +
          std::real(perturb_->si[nh]*exp_ksi_[points_[i].x_ind]);
      if (points_[i].t > 1.0f)
        points_[i].t = 1.0f;
    }
}

//================================
// Drawing StartDrawing
//================================

void animation::Drawing::StartDrawing(int argc, char **argv) {
  anim = std::unique_ptr<Animation> (new Animation());
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
  glutInitWindowSize(1000, 800);
  glutInitWindowPosition(100, 100);
  glutCreateWindow("Polhauzen flow stability");
  initialization();
  glutReshapeFunc(ReshapeFunction);
  glutDisplayFunc(draw);
  timer();
  glutMainLoop();
}

//================================
// Drawing ctor
//================================

animation::Drawing::Drawing() {}

//================================
// Drawing timer
//================================

void animation::Drawing::timer(int) {
  draw();
  anim->globaltime_ += dtau * timeScale;
  anim->updatePoints();
  glutTimerFunc(dtau, timer, 0);
}

//================================
// Drawing drawPoints
//================================

void animation::Drawing::drawPoints() {
  glBegin(GL_POINTS);
  for (size_t i = 0; i < anim->npoints_; ++i) {
      GLfloat temp = 1.0f - anim->points_[i].t;
      glColor3f(anim->points_[i].t, temp, temp);
      glVertex3f(0.0f, anim->points_[i].y,
          anim->xCoords_[anim->points_[i].x_ind] - 5.0f);
    }
  glEnd();
}

//================================
// Drawing draw
//================================

void animation::Drawing::draw() {
  glClear(GL_COLOR_BUFFER_BIT);
  glColor3f(1.0, 1.0, 1.0);
  glBegin(GL_LINES);
  glVertex3f(0.0, 0.0, 0.0);
  glVertex3f(0.0, 0.0, GLLength);
  glVertex3f(0.0, 0.0, 0.0);
  glVertex3f(0.0, yBegin, 0.0);
  glVertex3f(0.0, 0.0, 0.0);
  glVertex3f(5.0, 0.0, 0.0);
  glEnd();

  drawPoints();
  glutSwapBuffers();
}

//================================
// Drawing initialization
//================================

void animation::Drawing::initialization() {
  anim->setScales();
  anim->calcBasicflow();
  anim->calcPerturb();
  anim->setPoints();
}

//================================
// Drawing ResharpeFunction
//================================

void animation::Drawing::ReshapeFunction(GLsizei w, GLsizei h) {
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60, (GLfloat)w / (GLfloat)h, 1, 60);
  glTranslatef(0.0, 0.0, -30);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glRotatef(-80, 1, 0, 0);
  glRotatef(-100, 0, 0, 1);
  glTranslatef(0.0, 0.0, -5.0);
  glScalef(1.0, 1.9, 1.9);
}
}  // namespace conv_flow;
