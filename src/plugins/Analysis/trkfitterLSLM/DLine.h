#ifndef _DLINE_H_
#define _DLINE_H_

#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Vector/ThreeVector.h>
using namespace CLHEP;

class DLine {
 public:
  DLine(); // default constructor
  DLine(double x, double y, double z, double theta, double phi, int level = 0); // constructor with explicit point and direction
  DLine(Hep3Vector r0_in, Hep3Vector r1_in, int level = 0); // constructor with 3-vector positions of two points
  double doca(HepVector point);
  double doca(Hep3Vector point);
  HepVector poca(HepVector point);
  Hep3Vector poca(Hep3Vector point);
 private:
  Hep3Vector r0, r1;
  int debug_level;
};

#endif // _DLINE_H_
