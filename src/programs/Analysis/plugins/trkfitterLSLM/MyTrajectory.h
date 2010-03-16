#ifndef _MYTRAJECTORY_H_
#define _MYTRAJECTORY_H_

#include <iostream>
#include <iomanip>
#include <CLHEP/Vector/LorentzVector.h>
#include "TRACKING/DMCTrackHit.h"
#include "DLine.h"

#define DIST_BIG 10000.0
#define MAX_ITERATIONS 100000

using namespace std;

class MyTrajectory {
 public:
  MyTrajectory(int level = 1);
  virtual ~MyTrajectory();
  void clear();
  virtual void swim(HepVector startingPoint, double theta, double phi);
  virtual void swim(const HepVector &startingVector);
  void swim(const vector<double> &startingStdVector);
  void swimMC(vector<const DMCTrackHit*> &mctrackhits);
  void print();
  vector<HepLorentzVector*>* getTrajectory();
  template<class C> double doca(C& spaceObject, HepLorentzVector &poca) const;
  void checkClear();
  virtual unsigned int getNumberOfParams();
  double dist(HepVector& point, int trajIndex) const;
  double dist(DLine& line, int trajIndex) const;
  void para_min(double yMinus, double yZero, double yPlus, double &xMinFrac,
		double &yMin) const;
  virtual vector<double> getDelta() {
    return delta;
  }
  int getXYT(double z, double &x, double &y, double &ct) const;
  void dump_ascii(ostream *trajFile, int tag);

 protected:
  vector<HepLorentzVector*> traj;

 private:
  unsigned int nparams; // number of parameters for trajectory
  vector<double> delta;
  int debug_level;
};

template<class C> double MyTrajectory::doca(C& spaceObject, HepLorentzVector &poca) const {
  unsigned int ilo, ihi, imid;
  ilo = 0;
  ihi = traj.size() - 1;
  double distlo, disthi, distmid, distdown, distup;
  distlo = dist(spaceObject, ilo);
  disthi = dist(spaceObject, ihi);
  if (distlo > DIST_BIG || disthi > DIST_BIG) {
    cout << "MyTrajectory::doca: end point of trajectory too far away " << distlo << " " << disthi << endl;
    int ierror = 3;
    throw ierror;
  }
  if (distlo < disthi) {
    imid = ilo + 1;
  } else {
    imid = ihi - 1;
  }
  distmid = dist(spaceObject, imid);
  if (debug_level >= 4) cout << setprecision(14) << "MyTrajectory::doca: initialize " << ilo << " " << imid << " " << ihi << " " << distlo << " " << distmid << " " << disthi << endl;
  if (isnan(distmid)) {
    cout << "MyTrajectory::doca: distmid is not a number: " << ilo << " " << imid << " " << ihi << " " << distlo << " " << distmid << " " << disthi << endl;
    int ierror = 2;
    throw ierror;
  }
  if (distmid >= distlo || distmid >= disthi) {
    cout << "MyTrajectory::doca: bad initialization of doca search: " << ilo << " " << imid << " " << ihi << " " << distlo << " " << distmid << " " << disthi << endl;
    int ierror = 1;
    throw ierror;
  }
  double x1, x2, y1, y2, xnew, distnew = 0;
  unsigned int inew;
  int iterations = 0;
  while (imid - ilo > 1 || ihi - imid > 1) {
    x1 = (double)(imid - ilo);
    x2 = -(double)(ihi - imid); // unsigned int's cannot be negative
    y1 = distmid - distlo;
    y2 = distmid - disthi;
    xnew = (double)imid
      - 0.5*(
	     (x1*x1*y2 - x2*x2*y1)
	     /
	     (x1*y2 - x2*y1)
	     );
    inew = (unsigned int)(xnew + 0.5);
    if (inew < ilo || inew > ihi) { // jumped out of the bracket
      cout << "MyTrajectory::doca: jumped out of the miminum bracket: " << ilo << " " << inew << " " << ihi << endl;
      int ierror = 4;
      throw ierror;
    }
    distnew = dist(spaceObject, inew);
    if (debug_level >= 4) cout << "MyTrajectory::doca: xnew = " << xnew << " inew = " << inew << " distnew = " << distnew << endl;
    if (inew == imid) {
      // Look on either side for minimum bracket
      distup = dist(spaceObject, imid + 1);
      distdown = dist(spaceObject, imid - 1);
      if (distup > distmid && distdown > distmid) {
	ihi = imid + 1;
	disthi = distup;
	ilo = imid - 1;
	distlo = distdown;
      } else if (distup > distmid) { // higher on the up side
	if (ihi - imid > 1) { // room on the up side
	  ihi = imid + 1;
	  disthi = dist(spaceObject, ihi);
	} else { // no room on the up side
	  imid--;
	  distmid = dist(spaceObject, imid);
	}
      } else { // higher on the down side
	if (imid - ilo > 1) { // room on the down side
	  ilo = imid - 1;
	  distlo = dist(spaceObject, ilo);
	} else { // no room on the down side
	  imid++;
	  distmid = dist(spaceObject, imid);
	}
      }
    } else {
      
      if (distnew < distmid) { // new point is the new lowest point
	if (inew < imid) { // if new index less than old mid point index
	  ihi = imid; // new high point is the old mid point
	  disthi = distmid;
	} else { // new index is greater than the old mid point index
	  ilo = imid; // new low point is the old mid point
	  distlo = distmid;
	}
	imid = inew; // new point is the new mid point
	distmid = distnew;
      } else { // old mid-point is still the lowest point
	if (inew < imid) { // if new index less than old mid point index
	  ilo = inew; // new low point is the new point
	  distlo = distnew;
	} else { // new index is greater than the old mid point index
	  ihi = inew; // new high point is the new point
	  disthi = distnew;
	}
      }

    }
    if (debug_level >= 4) cout << "MyTrajectory::doca, after housekeeping: " << ilo << " " << imid << " " << ihi << " " << distlo << " " << distmid << " " << disthi << endl;
    iterations++;
    if (iterations > MAX_ITERATIONS) {
    cout << "MyTrajectory::doca: too many iterations " << ilo << " " << imid << " " << ihi << " " << distlo << " " << distmid << " " << disthi << endl;
    int ierror = 4;
    throw ierror;
    }
  }
  double dlon, dmidn, dhin, dinterp;
  dlon = dist(spaceObject, ilo);
  dmidn = dist(spaceObject, imid);
  dhin = dist(spaceObject, ihi);
  double xMinFrac, yMin;
  para_min(dlon*dlon, dmidn*dmidn, dhin*dhin, xMinFrac, yMin);
  dinterp = sqrt(yMin);
  if (xMinFrac > 0.0) {
    poca = xMinFrac*(*traj[ihi]) + (1.0 - xMinFrac)*(*traj[imid]);
  } else {
    poca = -xMinFrac*(*traj[ilo]) + (1.0 + xMinFrac)*(*traj[imid]);
  }
  if (debug_level >= 4) cout << "doca final calc: ilo,imid,ihi = " << ilo << ',' << imid << ',' << ihi << " dists = " << dlon << "," << dmidn << "," << dhin << " dinterp = " << dinterp << endl;
  return dinterp;

}

#endif //  _MYTRAJECTORY_H_
