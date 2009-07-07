#ifndef _MYTRAJECTORYGRKUTA_H_
#define _MYTRAJECTORYGRKUTA_H_

#include "MyTrajectoryBfield.h"

class MyTrajectoryGrkuta : public MyTrajectoryBfield {
 public:
  MyTrajectoryGrkuta(const DMagneticFieldMap *bfield, int level = 1);
  MyTrajectoryGrkuta(const HepVector B_in, int level = 1);
  void swim(const HepVector &param); // cut and paste
  void swim(double startingXprime, double startingZ, double startingTheta, double startingPhi, double ptinv);
 private:
  int debug_level;
};

#endif // _MYTRAJECTORYGRKUTA_H_
