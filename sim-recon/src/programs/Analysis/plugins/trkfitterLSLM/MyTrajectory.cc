#include <iostream>
#include <vector>

using namespace std;

#include "MyTrajectory.h"

MyTrajectory::MyTrajectory(int level) : nparams(4), delta(4, 0.001),
					debug_level(level) {
  // explicit default constructor, set number of params for straight line
  //  cout << "MyTrajectory constructor called\n";
  //  cout << "nparams = " << nparams << endl;
  return;
}

MyTrajectory::~MyTrajectory() {
  clear();
  return;
}

unsigned int MyTrajectory::getNumberOfParams() {
  return nparams;
}

void MyTrajectory::clear() {
  if (traj.size() != 0) {
    for (vector<HepLorentzVector*>::iterator iVector = traj.begin();
	 iVector != traj.end();
	 iVector++) {
      delete *iVector;
    }
  }
  traj.clear();
  return;
}

void MyTrajectory::swim(HepVector startingPoint, double theta,
			double phi)
// default swim, straight line
{
  checkClear();
  double stepLength = 1.0;
  double t0 = 0.0;
  HepLorentzVector step;
  step.setX(stepLength*sin(theta)*cos(phi));
  step.setY(stepLength*sin(theta)*sin(phi));
  step.setZ(stepLength*cos(theta));
  step.setT(stepLength);
  HepLorentzVector* thisVector;
  thisVector = new HepLorentzVector(startingPoint(1), startingPoint(2),
				    startingPoint(3), t0);
  traj.push_back(thisVector);
  HepLorentzVector lastVector(*thisVector);
  for (int i = 0; i < 600; i++) {
    thisVector = new HepLorentzVector();
    *thisVector = lastVector + step;
    traj.push_back(thisVector);
    lastVector = *thisVector;
  }
  return;
}

void MyTrajectory::swim(const HepVector& startingVector) {
  HepVector startingPoint(3);
  startingPoint(1) = startingVector(1);
  startingPoint(2) = startingVector(2);
  startingPoint(3) = 0.0; // start at z = 0 every time
  double theta = startingVector(3);
  double phi = startingVector(4);
  swim(startingPoint, theta, phi);
  return;
}

void MyTrajectory::swim(const vector<double> &params) {
  int size = params.size();
  HepVector hepParams(size);
  for (int i = 0; i < size; i++) {
    hepParams(i + 1) = params[i];
  }
  swim(hepParams);
  return;
}

void MyTrajectory::swimMC(vector<const DMCTrackHit*> &mctrackhits) {
  checkClear();
  const DMCTrackHit* mchit;
  HepLorentzVector *point;
  for (vector<const DMCTrackHit*>::iterator imchit = mctrackhits.begin();
       imchit != mctrackhits.end();
       imchit++) {
    mchit = *imchit;
    if (mchit->system & (SYS_CDC | SYS_FDC) && mchit->primary == 1) {
      double r = mchit->r, phi = mchit->phi, z = mchit->z;
      point = new HepLorentzVector(r*cos(phi), r*sin(phi), z);
      traj.push_back(point);
    }
  }
}

void MyTrajectory::print()
{
  cout << "### MyTrajectory print out begin ###" << endl;
  for (vector<HepLorentzVector*>::iterator iVector = traj.begin();
       iVector != traj.end();
       iVector++) {
    cout << (**iVector).x() << " " << (**iVector).y() << " " << (**iVector).z() << " " << (**iVector).t()
    	 << endl;
  }
  cout << "### MyTrajectory print out end ###" << endl;
}

vector<HepLorentzVector*>* MyTrajectory::getTrajectory() {
  return &traj;
}

double MyTrajectory::dist(HepVector& point, int trajIndex) const {
  Hep3Vector delta, point3(point(1), point(2), point(3)), trajPoint;
  trajPoint = traj[trajIndex]->getV();
  delta = point3 - trajPoint;
  if (debug_level >= 4) cout << "point3 = " << point3
			     << "traj point = " << trajPoint
			     << "delta = " << delta
			     << "delta.mag = " << delta.mag() << endl;
  return delta.mag();
}

double MyTrajectory::dist(DLine& line, int trajIndex) const {
  return line.doca(*traj[trajIndex]);
}

int MyTrajectory::getXYT(double z, double &x, double &y, double &t) const {
  int iBefore = 0, iAfter = traj.size() - 1, iTry;
  double zBefore = traj[iBefore]->z();
  double zAfter = traj[iAfter]->z();
  double zTry;
  if (z < zBefore || z > zAfter) {
    return 1;
  }
  while (iAfter - iBefore > 1) {
    iTry = iBefore
      + (int)((double)(iAfter - iBefore)*(z - zBefore)/(zAfter - zBefore)
	      + 0.5);
    if (debug_level > 3) cout << iBefore << ' ' << iTry
			      << ' ' << iAfter << endl;
    if (iBefore == iTry) iTry++;
    if (iAfter == iTry) iTry--;
    if (debug_level > 3) cout << iBefore << ' ' << iTry
			      << ' ' << iAfter << endl;
    zTry = traj[iTry]->z();
    if (debug_level > 3) cout << "zTry = " << zTry << endl;
    if (z < zTry) {
      iAfter = iTry;
      zAfter = traj[iAfter]->z();
    } else if (z > zTry) {
      iBefore = iTry;
      zBefore = traj[iBefore]->z();
    } else {
      iBefore = iTry;
      zBefore = traj[iBefore]->z();
      iAfter = iTry + 1;
      zAfter = traj[iAfter]->z();
    }
    if (debug_level > 3) cout << z << ' ' << zBefore << ' ' << zTry
			      << ' ' << zAfter << endl;
  }
  double frac, otherfrac;
  frac = (z - zBefore)/(zAfter - zBefore);
  otherfrac = 1.0 - frac;
  if (debug_level > 3) cout << frac << ' ' << otherfrac << endl;
  if (debug_level > 3) cout << "x before, after " << traj[iBefore]->x()
			    << ' '<< traj[iAfter]->x() << endl;
  if (debug_level > 3) cout << "y before, after " << traj[iBefore]->y()
			    << ' '<< traj[iAfter]->y() << endl;
  x = frac*traj[iAfter]->x() + otherfrac*traj[iBefore]->x();
  y = frac*traj[iAfter]->y() + otherfrac*traj[iBefore]->y();
  t = frac*traj[iAfter]->t() + otherfrac*traj[iBefore]->t();
  if (debug_level > 3) cout << "MyTrajectory::getXYT, x, y, t = " << x << " " << y
			    << " " << t << endl;
  return 0;
}

void MyTrajectory::para_min(double yMinus, double yZero, double yPlus,
			    double &xMinFrac, double &yMin) const {
  double a, b, c;
  a = 0.5*(yPlus - 2.0*yZero + yMinus);
  b = 0.5*(yPlus - yMinus);
  c = yZero;
  yMin = -b*b/(4.0*a) + c;
  if (yMin < 0.0) yMin = 0.0;
  xMinFrac = -b/(2.0*a);
  return;
}

void MyTrajectory::checkClear() {
  if (traj.size() != 0) {
    int ierror = 2;
    throw ierror;
  }
  return;
}

void MyTrajectory::dump_ascii(ostream *trajFile, int tag) {
  HepVector trajPoint(3);
  for (unsigned int i = 0; i < traj.size(); i++) {
    trajPoint = *(traj[i]);
    *trajFile << tag << " " << i + 1 << " " << trajPoint(1) << " "
	      << trajPoint(2) << " " << trajPoint(3) << endl;
  }
}
