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
				     MyTrajectory *trajectory,
				     const DLorentzDeflections *lorentz_def_in,
				     int level) : 
  n_fdc(pseudopoints->size()), n_cdc(trackHits->size()), ppPtr(pseudopoints),
  trkhitPtr(trackHits), trajPtr(trajectory), delta(trajPtr->getDelta()),
  debug_level(level), lorentz_def(lorentz_def_in), storeDetails(false),
  innerResidFrac(1.0), rCDC(trackHits, trajectory, level),
  rFDC(pseudopoints, trajectory, lorentz_def_in, level)
{}

void combinedResidFunc::resid(const HepVector *x, void *data, HepVector *f){
// input parameters
//   x: pointer to a vector of fit parameters
//   data: ???
// output parameters
//   f: pointer to vector of residuals
  if (debug_level > 2) {
    cout << "combinedResidFunc::resid: resid called\n";
    cout << "                          params: " << *x;
  }
  // do a swim with the input parameters
  trajPtr->swim(*x);
  // populate f vector with residuals

  double thisChiSquared = 0.0;
  double thisResid;

  // get info from residFDC class

  rFDC.calcResids();
  vector<double> residsF;
  rFDC.getResids(residsF);
  FDCHitDetails *FDCHitDetailsPtr;
  vector<HepVector> pointF;
  vector<double> docasF, errorsF;
  vector<HepLorentzVector> pocasF;
  rFDC.getDetails(pointF, docasF, errorsF, pocasF);
  for (unsigned int ir = 0; ir < n_fdc; ir++) {
    thisResid = residsF[ir];
    (*f)(ir + 1) = thisResid;
    if (storeDetails) {
      thisChiSquared += thisResid*thisResid;
      FDCHitDetailsPtr = new FDCHitDetails();
      FDCHitDetailsPtr->doca = docasF[ir];
      FDCHitDetailsPtr->poca = pocasF[ir];
      FDCHitDetailsPtr->rCorr = pointF[ir];
      FDCDetails.push_back(FDCHitDetailsPtr);
    }
  }

  // get info from CDC residual class

  rCDC.setInnerResidFrac(innerResidFrac);
  rCDC.calcResids();
  vector<double> residsC;
  rCDC.getResids(residsC);
  CDCHitDetails *CDCHitDetailsPtr;
  vector<double> docasC, distsC, errorsC;
  vector<HepLorentzVector> pocasC;
  vector<HepVector> posWiresC;
  rCDC.getDetails(docasC, distsC, errorsC, pocasC, posWiresC);
  for (unsigned int ir = 0; ir < n_cdc; ir++) {
    thisResid = residsC[ir];
    (*f)(n_fdc + ir + 1) = thisResid;
    if (storeDetails) {
      thisChiSquared += thisResid*thisResid;
      CDCHitDetailsPtr = new CDCHitDetails();
      CDCHitDetailsPtr->doca = docasC[ir];
      CDCHitDetailsPtr->poca = pocasC[ir];
      CDCHitDetailsPtr->dist = distsC[ir];
      CDCHitDetailsPtr->posWire = posWiresC[ir];
      CDCDetails.push_back(CDCHitDetailsPtr);
    }
  }

  if (debug_level > 2) cout << "combinedResidFunc::resid: resids:" << *f;
  if (storeDetails) chiSquared = thisChiSquared;

  // clear the trajectory
  trajPtr->clear();

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
  if (debug_level >= 4) {
    cout << "combinedResidFunc::pseudo2HepVector, x = " << x << " y = " << y << " z = " << z << " ct = " << ct << " delta_x = " << delta_x << " delta_y = " << delta_y << " ispos = " << ispos << " point = " << point << endl;
  }
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
  if (debug_level > 3) cout << setprecision(14)
			    << "combinedResidFunc::getCorrectionSign,"
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

void combinedResidFunc::getResidsBoth(vector<double> &residsBoth) {

  rCDC.setInnerResidFrac(innerResidFrac);
  rCDC.calcResids();
  vector<double> residsC;
  rCDC.getResids(residsC);

  rFDC.calcResids();
  vector<double> residsF;
  rFDC.getResids(residsF);

  residsBoth.reserve(n_fdc + n_cdc); // reserve slots for FDC and CDC residuals

  for (unsigned int i = 0; i < n_fdc; i++) {
    residsBoth[i] = residsF[i];
  }

  for (unsigned int i = 0; i < n_cdc; i++) {
    residsBoth[n_fdc + i] = residsC[i];
  }

  return;

}

void combinedResidFunc::deriv(const HepVector *params, void *data, HepMatrix *Jacobian) {
  HepVector paramsCentral = *params, paramsThis;
  int nResids = n_fdc + n_cdc;
  unsigned int nParams = trajPtr->getNumberOfParams();
  int iHep = 0, jHep = 0;
  if (debug_level > 2){
    for (unsigned int j = 0; j < nParams; j++) {
      jHep = j + 1;
      cout << "central params: jHep = " << jHep << ", values:" << paramsCentral(jHep) << endl;
    }
  }
  // do central swim
  trajPtr->swim(paramsCentral);
  // save central residuals
  vector<double> residsCentral, residsThis;
  getResidsBoth(residsCentral);
  if (debug_level >=3) {
    cout << "combinedResidFunc::deriv2: central resids, ";
    for (int k = 0; k < nResids; k++) {
      cout << k << "=" << residsCentral[k] << " ";
    }
    cout << endl;
  }
  trajPtr->clear();
  // calculate the Jacobian matrix
  for (unsigned int j = 0; j < nParams; j++) {
    jHep = j + 1;
    // prepare perturbed parameters
    paramsThis = paramsCentral;
    paramsThis(jHep) = paramsCentral(jHep) + delta[j];
    if (debug_level > 2) cout << "perturbed params: jHep = " << jHep << ", values:" << paramsThis << endl;
    // do the perturbed swim
    trajPtr->swim(paramsThis);
    // get the perturbed residuals
    getResidsBoth(residsThis);
    if (debug_level >=3) {
      cout << "combinedResidFunc::deriv2: resids, ";
      for (int k = 0; k < nResids; k++) {
	cout << k << "=" << residsThis[k] << " ";
      }
      cout << endl;
    }
    // calculate the derivatives
    for (int i = 0; i < nResids; i++) {
      iHep = i + 1;
      (*Jacobian)(iHep, jHep) = (residsThis[i] - residsCentral[i])/delta[j];
    } 
    trajPtr->clear();
  }
  return;
}
