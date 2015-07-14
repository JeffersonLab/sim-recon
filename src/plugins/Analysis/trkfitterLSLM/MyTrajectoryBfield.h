#ifndef _MYTRAJECTORYBFIELD_H_
#define _MYTRAJECTORYBFIELD_H_

#include "HDGEOMETRY/DMagneticFieldMap.h"
#include "MyTrajectory.h"

class MyTrajectoryBfield : public MyTrajectory {
 public:
  MyTrajectoryBfield(const DMagneticFieldMap *bfield, int level = 1);
  MyTrajectoryBfield(const HepVector B_in, int level = 1);
  void swim(const HepVector &param);
  void swim(double startingXprime, double startingZ, double startingTheta, double startingPhi, double ptinv);
  unsigned int getNumberOfParams();
  virtual vector<double> getDelta() {
    return delta;
  }
  HepVector getParams();
  HepVector getField(HepVector& r);
 protected:
  unsigned int nparams;
  HepVector params;
 private:
  HepVector BConst;
  vector<double> delta;
  const DMagneticFieldMap *bfield;
  int debug_level;
};

#endif // _MYTRAJECTORYBFIELD_H_
