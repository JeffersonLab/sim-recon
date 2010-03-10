// residCDC: class to calculate drift distance, distance of closest
// approach, point of closest approach, residual and error on residual
// for the CDC for a given set of track parameters

#include "residCDC.h"

using namespace std;

const double velDrift = 55.e-4; // cm/ns
const double c = 29.9792548; // speed of light, cm/ns

// constructor, takes pointer to vector of trackhits and a pointer
// to a trajectory
residCDC::residCDC(vector<const DCDCTrackHit*> *trackHits,
				     const MyTrajectory *trajectory, int level) : 
  n_cdc(trackHits->size()), trkhitVectorPtr(trackHits), trajPtr(trajectory),
  debug_level(level), errorCDC(0.018)
{}

void residCDC::calcResids(){
  const DCDCTrackHit* trkhitPtr;
  DLine line;
  double docaThis, distThis, residThis, errorThis;
  HepLorentzVector pocaThis;
  HepVector posWireThis;
  doca.clear();
  dist.clear();
  poca.clear();
  posWire.clear();
  error.clear();
  resid.clear();
  for (unsigned int j = 0; j < n_cdc; j++) {
    trkhitPtr = (*trkhitVectorPtr)[j];
    line = trackhit2line(*trkhitPtr);
    docaThis = trajPtr->doca(line, pocaThis);
    distThis = velDrift*(trkhitPtr->tdrift - pocaThis.getT()/c);
    if (docaThis > distThis) {
      residThis = (distThis - docaThis)/errorCDC;
    } else {
      residThis = innerResidFrac*(distThis - docaThis)/errorCDC;
    }
    posWireThis = line.poca(pocaThis);
    if (debug_level > 2) cout << "residCDC: j = " << j << " dist = " << distThis << " doca = " << docaThis << " poca xyzt = " << pocaThis.getX() << ' ' << pocaThis.getY() << ' ' << pocaThis.getZ() << ' ' << pocaThis.getT()/c << " resid = " << residThis << endl;
    errorThis = errorCDC;
    doca.push_back(docaThis);
    dist.push_back(distThis);
    poca.push_back(pocaThis);
    posWire.push_back(posWireThis);
    resid.push_back(residThis);
    error.push_back(errorThis);
  }
};

DLine residCDC::trackhit2line(const DCDCTrackHit &trkhit) {
  const DCDCWire *wire = trkhit.wire;
  double x = wire->origin.X();
  double y = wire->origin.Y();
  double z = wire->origin.Z();
  double theta = acos(wire->udir.z());
  double phi = atan2(wire->udir.y(), wire->udir.x());
  if (debug_level >= 4) cout << setprecision(14) << "residCDC::trackhit2line: x = " << x << " y = " << y << " z = " << z << " theta " << theta << " phi " << phi << endl;
  DLine line(x, y, z, theta, phi, debug_level);
  return line;
}

void residCDC::setInnerResidFrac(double innerResidFracIn) {
  innerResidFrac = innerResidFracIn;
}

void residCDC::getResids(vector<double> &residsRef) {
  residsRef = resid;
  return;
}

void residCDC::getDetails(vector<double> &docasRef, vector<double> &distsRef,
			  vector<double> &errorsRef,
			  vector<HepLorentzVector> &pocasRef,
			  vector<HepVector> &posWiresRef) {
  docasRef = doca;
  distsRef = dist;
  errorsRef = error;
  pocasRef = poca;
  posWiresRef = posWire;
  return;
}
