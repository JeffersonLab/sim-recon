#include "DLine.h"
using namespace std;

#define DELTA_R = 10000.0

// default constructor line along z-axis
DLine::DLine() : r0(0.0, 0.0, 0.0), r1(0.0, 0.0, 1.0), debug_level(0) {}

DLine::DLine(Hep3Vector r0_in, Hep3Vector r1_in, int level): r0(r0_in), r1(r1_in), debug_level(level) {};

DLine::DLine(double x, double y, double z, double theta, double phi, int level) : debug_level(level) {
  r0(0) = x;
  r0(1) = y;
  r0(2) = z;
  double sinth = sin(theta);
  r1(0) = x + sinth*cos(phi);
  r1(1) = y + sinth*sin(phi);
  r1(2) = z + cos(theta);
  if (debug_level > 2) cout << "DLine::DLine: line constructed: " << r0 << r1 << endl;
  return;
}

double DLine::doca(HepVector point) {
  Hep3Vector point3(point(1), point(2), point(3));
  return doca(point3);
}

double DLine::doca(Hep3Vector point) {
  Hep3Vector num, diff;
  diff = r1 - r0;
  num = diff.cross(r0 - point);
  double doca = num.mag()/diff.mag();
  if (debug_level >= 4) cout << "DLine:doca: num = " << num.mag() << " diff = " << diff.mag() << " doca = " << doca << endl;
  return doca;
};

HepVector DLine::poca(HepVector point) {
  Hep3Vector point3(point(1), point(2), point(3)), poca3;
  poca3 = poca(point3);
  HepVector poca(3);
  poca(1) = poca3(0);
  poca(2) = poca3(1);
  poca(3) = poca3(2);
  return poca;
}

Hep3Vector DLine::poca(Hep3Vector point) {
  double t;
  Hep3Vector diff;
  diff = r1 - r0;
  t = (point - r0)*diff/diff.mag2();
  Hep3Vector poca = r0 + diff*t;
  // cout << point << r0 << r1 << diff << t << poca << endl;
  return poca;
};

