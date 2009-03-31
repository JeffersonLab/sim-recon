#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/Vector.h>

#define SPEED_OF_LIGHT 29.9792548

using namespace std;

#include "FDC/DFDCPseudo_factory.h"
#include "MyTrajectory.h"
#include "residFunc.h"
#include "combinedResidFunc.h"
#include "HDGEOMETRY/DLorentzDeflections.h"

const double velDrift = 55.e-4; // cm/ns
const double c = 29.9792548; // speed of light, cm/ns

// constructor, takes pointer to vector of pseudopoints and a pointer
// to a trajectory
combinedResidFunc::combinedResidFunc(vector<const DFDCPseudo*> *pseudopoints,
				     vector<const DCDCTrackHit*> *trackHits,
				     MyTrajectory *trajectory, const DLorentzDeflections *lorentz_def_in, int level) : 
  n_fdc(pseudopoints->size()), n_cdc(trackHits->size()), ppPtr(pseudopoints),
  trkhitPtr(trackHits), trajPtr(trajectory), delta(trajPtr->getDelta()),
  debug_level(level), lorentz_def(lorentz_def_in), storeDetails(false),
  innerResidFrac(1.0)
{}

void combinedResidFunc::resid(const HepVector *x, void *data, HepVector *f){
// input parameters
//   x: pointer to a vector of fit parameters
//   data: ???
// output parameters
//   f: pointer to vector of residuals
  double doca;
  if (debug_level > 2) {
    cout << "combinedResidFunc::resid: resid called\n";
    cout << "                          params: " << *x;
  }
  // do a swim with the input parameters
  trajPtr->swim(*x);
  // populate f vector with residuals
  HepVector point(3);
  HepLorentzVector poca;
  FDCHitDetails *FDCHitDetailsPtr;
  for (unsigned int i = 0; i < n_fdc; i++) {
    point = pseudo2HepVector(*((*ppPtr)[i]));
    (*f)(i + 1) = trajPtr->doca(point, poca)/ERROR_FDC;
    if (storeDetails) {
      FDCHitDetailsPtr = new FDCHitDetails();
      *FDCHitDetailsPtr = getDetails((*ppPtr)[i], point);
      FDCDetails.push_back(FDCHitDetailsPtr);
    }
  }
  DLine line;
  double dist;
  CDCHitDetails *CDCHitDetailsPtr;
  double thisChiSquared = 0.0;
  double thisResid;
  for (unsigned int j = 0; j < n_cdc; j++) {
    if (debug_level > 2) cout << "working on cdc hit " << j << endl;
    line = trackhit2line(*((*trkhitPtr)[j]));
    doca = trajPtr->doca(line, poca);
    dist = velDrift*((*trkhitPtr)[j]->tdrift - poca.getT()/c);
    if (debug_level > 2) cout << "resid, cdc: j = " << j << " dist = " << dist << " doca = " << doca << " poca xyzt = " << poca.getX() << ' ' << poca.getY() << ' ' << poca.getZ() << ' ' << poca.getT()/c << " resid = " << dist - doca << endl;
    if (isnan(dist)) {
      thisResid = 0.0;
    } else {
      if (doca > dist) {
	thisResid = (dist - doca)/ERROR_CDC;
      } else {
	thisResid = innerResidFrac*(dist - doca)/ERROR_CDC;
      }
    }
    (*f)(n_fdc + j + 1) = thisResid;
    if (storeDetails) {
      thisChiSquared += thisResid*thisResid;
      CDCHitDetailsPtr = new CDCHitDetails();
      *CDCHitDetailsPtr = getDetails((*trkhitPtr)[j], line);
      CDCDetails.push_back(CDCHitDetailsPtr);
    }
  }
  if (debug_level > 2) cout << "combinedResidFunc::resid: resids:" << *f;
  if (storeDetails) chiSquared = thisChiSquared;
  trajPtr->clear();
};

void combinedResidFunc::deriv(const HepVector *x, void *data, HepMatrix *J){
  if (debug_level > 2) {
    cout << "combinedResidFunc::deriv: deriv called\n";
    cout << "                          params: " << *x << endl;
  }
  HepLorentzVector poca(3);
  // save base parameters
  HepVector xBase = *x;
  // do central swim
  trajPtr->swim(xBase);
  // store pseudo points as three vectors
  vector<HepVector *>pPoints;
  HepVector *thisPointPtr;
  for (unsigned int i = 0; i < n_fdc; i++) {
    thisPointPtr = new HepVector(3);
    *thisPointPtr = pseudo2HepVector(*((*ppPtr)[i]));
    pPoints.push_back(thisPointPtr);
  }
  // store track hits as lines
  vector<DLine *> linePtrs;
  DLine *thisLinePtr;
  for (unsigned int j = 0; j < n_cdc; j++) {
    thisLinePtr = new DLine();
    *thisLinePtr = trackhit2line(*((*trkhitPtr)[j]));
    linePtrs.push_back(thisLinePtr);
  }
  // store base residuals
  HepVector residBase(n_fdc + n_cdc);
  // base resids for FDC
  for (unsigned int i = 0; i < n_fdc; i++) {
    residBase(i + 1) = trajPtr->doca(*(pPoints[i]), poca)/ERROR_FDC;
  }
  // base resids for CDC
  double docaThis, distThis;
  for (unsigned int j = 0; j < n_cdc; j++) {
    distThis = velDrift*((*trkhitPtr)[j]->tdrift - poca.getT()/c);
    docaThis = trajPtr->doca(*(linePtrs[j]), poca);
    if (docaThis > distThis) {
      residBase(n_fdc + j + 1) = (distThis - docaThis)/ERROR_CDC;
    } else {
      residBase(n_fdc + j + 1) = innerResidFrac*(distThis - docaThis)/ERROR_CDC;
    }
  }
  if (debug_level > 2) cout << "base resids:" << residBase;
  trajPtr->clear();
  // calculate Jacobian
  unsigned int p = trajPtr->getNumberOfParams();
  HepVector xThis(p);
  int iHep, jHep; // index for HepVector () notation
  for (unsigned int i = 0; i < p; i++) {
    iHep = i + 1;
    xThis = xBase; // set params back to base
    xThis(iHep) = xBase(iHep) + delta[i];
    if (debug_level > 2) cout << "perturbed params: iHep = " << iHep << ", values:" << xThis << endl;
    // do the perturbed swim
    trajPtr->swim(xThis);
    // calculate derivatives for FDC points
    for (unsigned int j = 0; j < n_fdc; j++) {
      jHep = j + 1;
      docaThis = trajPtr->doca(*(pPoints[j]), poca)/ERROR_FDC;
      if (debug_level > 2) cout << "resid " << j << " = " << docaThis << endl;
      (*J)(jHep, iHep) = (docaThis - residBase(jHep))/delta[i];
    }
    // calculate derivatives for CDC points
    for (unsigned int j = 0; j < n_cdc; j++) {
      jHep = n_fdc + j + 1;
      distThis = velDrift*((*trkhitPtr)[j]->tdrift - poca.getT()/c);
      docaThis = trajPtr->doca(*(linePtrs[j]), poca);
      if (debug_level > 2) cout << j << " dist = " << distThis << " doca = " << docaThis << " resid  = " << distThis - docaThis << endl;
      if (isnan(distThis)) {
	(*J)(jHep, iHep) = 0;
      } else {
	if (docaThis > distThis) {
	  (*J)(jHep, iHep) = ((distThis - docaThis)/ERROR_CDC - residBase(jHep))/delta[i];
	} else {
	  (*J)(jHep, iHep) = (innerResidFrac*(distThis - docaThis)/ERROR_CDC - residBase(jHep))/delta[i];
	}
      }
    }
    trajPtr->clear();
  }
  if (debug_level >= 3) {
    for (unsigned int i = 0; i < p; i++) {
      iHep = i + 1;
      cout << iHep;
      for (unsigned int j = 0; j < n_fdc + n_cdc; j++) {
	jHep = j + 1;
	cout << ' ' << (*J)(jHep, iHep);
      }
      cout << endl;
    }
  }
  for (unsigned int i = 0; i < n_fdc; i++) {
    delete pPoints[i];
  }
  for (unsigned int j = 0; j < n_cdc; j++) {
    delete linePtrs[j];
  }
};

void combinedResidFunc::residAndDeriv(const HepVector *x, void *data, HepVector *f,
		   HepMatrix *J){};

HepVector combinedResidFunc::pseudo2HepVector(const DFDCPseudo &ppoint) {
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
  return point;
}

bool combinedResidFunc::getCorrectionSign(const DFDCPseudo &ppoint, double x, double y, double delta_x, double delta_y) {
  double xWire = ppoint.wire->udir(0);
  double yWire = ppoint.wire->udir(1);
  double wireCrossTraj = xWire*(y - ppoint.y) - yWire*(x - ppoint.x);
  double wireCrossDelta = xWire*delta_y - yWire*delta_x;
  bool isposTraj = wireCrossTraj > 0?true:false;
  bool isposDelta = wireCrossDelta > 0?true:false;
  bool ispos = !(isposTraj ^ isposDelta);
  if (debug_level > 3) cout << " ppx = " << ppoint.x
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

void combinedResidFunc::getCorrectionValue(const DFDCPseudo &ppoint, double x, double y, double z, double ct, double &delta_x, double &delta_y) {
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

DLine combinedResidFunc::trackhit2line(const DCDCTrackHit &trkhit) {
  const DCDCWire *wire = trkhit.wire;
  double x = wire->origin.X();
  double y = wire->origin.Y();
  double z = wire->origin.Z();
  //  double phi_wire = atan2(y, x);
  //  double phi_naive = phi_wire + PIOVER2;
  //  double theta_naive = wire->stereo;
  double theta = acos(wire->udir.z());
  double phi = atan2(wire->udir.y(), wire->udir.x());
  /*
  cout << "theta " << theta_naive << " " << theta << " phi " << phi_naive
       << " "<< phi << endl;
  */
  DLine line(x, y, z, theta, phi);
  return line;
}

FDCHitDetails combinedResidFunc::getDetails(const DFDCPseudo *ppoint, HepVector point) {
  FDCHitDetails details;
  details.doca = trajPtr->doca(point, details.poca);
  details.rCorr = point;
  return details;
}

CDCHitDetails combinedResidFunc::getDetails(const DCDCTrackHit *trackhit, DLine line) {CDCHitDetails details;
  HepLorentzVector poca;
  details.doca = trajPtr->doca(line, poca);
  details.poca = poca;
  details.dist = velDrift*(trackhit->tdrift - poca.getT()/c);
  details.posWire = line.poca(poca);
  return details;
}

void combinedResidFunc::setStoreDetails(bool value) {
  storeDetails = value;
  return;
}

void combinedResidFunc::clearDetails() {
  for (unsigned int i = 0; i < FDCDetails.size(); i++) {delete FDCDetails[i];}
  FDCDetails.clear();
  for (unsigned int i = 0; i < CDCDetails.size(); i++) {delete CDCDetails[i];}
  CDCDetails.clear();
}

void combinedResidFunc::setInnerResidFrac(double innerResidFracIn) {
  innerResidFrac = innerResidFracIn;
}
