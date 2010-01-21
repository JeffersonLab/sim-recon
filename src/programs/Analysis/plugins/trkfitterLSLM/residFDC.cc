#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/Vector.h>

using namespace std;

#include "FDC/DFDCPseudo_factory.h"
#include "MyTrajectory.h"
#include "HDGEOMETRY/DLorentzDeflections.h"
#include "residFDC.h"

const double velDrift = 55.e-4; // cm/ns
const double c = 29.9792548; // speed of light, cm/ns

// constructor, takes pointer to vector of pseudopoints and a pointer
// to a trajectory
residFDC::residFDC(vector<const DFDCPseudo*> *pseudopoints,
		   const MyTrajectory *trajectory,
		   const DLorentzDeflections *lorentz_def_in, int level) : 
  n_fdc(pseudopoints->size()), ppPtr(pseudopoints),
  trajPtr(trajectory), debug_level(level), lorentz_def(lorentz_def_in), errorFDC(0.046) {}

void residFDC::calcResids() {
  double docaThis, errorThis, residThis;
  HepVector pointThis(3);
  HepLorentzVector pocaThis;
  const DFDCPseudo* ppointPtr;
  doca.clear();
  poca.clear();
  error.clear();
  resid.clear();
  for (unsigned int i = 0; i < n_fdc; i++) {
    ppointPtr = (*ppPtr)[i];
    pointThis = pseudo2HepVector(*ppointPtr);
    point.push_back(pointThis);
    docaThis = trajPtr->doca(pointThis, pocaThis);
    errorThis = errorFDC;
    residThis = docaThis/errorThis;
    doca.push_back(docaThis);
    poca.push_back(pocaThis);
    resid.push_back(residThis);
    error.push_back(errorThis);
    if (debug_level > 2) cout << "residFDC: i = " << i << " doca = " << docaThis << " poca xyzt = " << pocaThis.getX() << ' ' << pocaThis.getY() << ' ' << pocaThis.getZ() << ' ' << pocaThis.getT()/c << " resid = " << residThis << endl;

  }
}

HepVector residFDC::pseudo2HepVector(const DFDCPseudo &ppoint) {
  double x;
  double y;
  double ct;
  double z = ppoint.wire->origin(2);
  trajPtr->getXYT(z, x, y, ct); // on trajectory
  double delta_x = 0.0, delta_y = 0.0;
  getCorrectionValue(ppoint, x, y, z, ct, delta_x, delta_y);
  bool ispos = getCorrectionSign(ppoint, x, y, delta_x, delta_y);
  HepVector point(3);
  if (ispos) {
    point(1) = ppoint.x + delta_x;
    point(2) = ppoint.y + delta_y;
  } else {
    point(1) = ppoint.x - delta_x;
    point(2) = ppoint.y - delta_y;
  }
  point(3) = z; 
  if (debug_level >= 4) {
    cout << "residFDC::pseudo2HepVector, x = " << x << " y = " << y << " z = " << z << " ct = " << ct << " delta_x = " << delta_x << " delta_y = " << delta_y << " ispos = " << ispos << " point = " << point << endl;
  }
  return point;
}

bool residFDC::getCorrectionSign(const DFDCPseudo &ppoint, double x, double y, double delta_x, double delta_y) {
  double xWire = ppoint.wire->udir(0);
  double yWire = ppoint.wire->udir(1);
  double wireCrossTraj = xWire*(y - ppoint.y) - yWire*(x - ppoint.x);
  double wireCrossDelta = xWire*delta_y - yWire*delta_x;
  bool isposTraj = wireCrossTraj > 0?true:false;
  bool isposDelta = wireCrossDelta > 0?true:false;
  bool ispos = !(isposTraj ^ isposDelta);
  if (debug_level > 3) cout << setprecision(14)
			    << "residFDC::getCorrectionSign,"
			    << " x = " << x
			    << " y = " << y
			    << " ppx = " << ppoint.x
			    << " ppy = " << ppoint.y
			    << " dx = " << x - ppoint.x
			    << " dy = " << y - ppoint.y
			    << " delta_x = " << delta_x
			    << " delta_y = " << delta_y
			    << " xWire = " << xWire
			    << " yWire = " << yWire
			    << " wireCrossTraj = " << wireCrossTraj
			    << " wireCrossDelta = " << wireCrossDelta
			    << " isposTraj = " << isposTraj
			    << " isposDelta = " << isposDelta
			    << " ispos = " << ispos
			    << endl;
  return ispos;
}

void residFDC::getCorrectionValue(const DFDCPseudo &ppoint, double x, double y, double z, double ct, double &delta_x, double &delta_y) {
  double driftDist = (ppoint.time - ct/c)*DRIFT_VELOCITY;
  //double driftDist = ppoint.time*DRIFT_VELOCITY;
  double cosangle = ppoint.wire->udir(1);
  double sinangle= ppoint.wire->udir(0);
  double alpha = 0.0; // for now
  if (debug_level > 3) cout << "x, y, z, ct = " << x << ' ' << y << ' ' << z << ' ' << ct << " time, ct/c = " << ppoint.time << " " << ct/c << " driftdist = " << driftDist << endl;
  double lorentzShift = lorentz_def->GetLorentzCorrection(x, y, z, alpha, driftDist);
  double ds = -lorentzShift;
  double dw = driftDist;
  if (debug_level > 3) cout << "ds, dw" << ds << " " << dw << endl;
  delta_x = dw*cosangle + ds*sinangle;
  delta_y = -dw*sinangle + ds*cosangle;
  return;
}

void residFDC::getResids(vector<double> &residRef) {
  residRef = resid;
  return;
}

void residFDC::getDetails(vector<HepVector> &pointsRef, vector<double> &docasRef, vector<double> &errorsRef,
			  vector<HepLorentzVector> &pocasRef) {
  pointsRef = point;
  docasRef = doca;
  errorsRef = error;
  pocasRef = poca;
  return;
}

// end of C++ source
