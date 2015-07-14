// residFDCCathode: class to calculate drift distance, distance of closest
// approach, point of closest approach, residual and error on residual
// for the FDC anode wires for a given set of track parameters

#include "residFDCCathode.h"

using namespace std;

const double velDrift = 55.e-4; // cm/ns
const double c = 29.9792548; // speed of light, cm/ns

// constructor, takes pointer to vector of trackhits and a pointer
// to a trajectory
residFDCCathode::residFDCCathode(vector<const DFDCPseudo*> *pseudopoints,
				     const MyTrajectory *trajectory, int level) : 
  n_fdca(pseudopoints->size()), pseudopointVectorPtr(pseudopoints), trajPtr(trajectory),
  debug_level(level), errorFDCA(0.018)
{}

void residFDCCathode::calcResids(){
  const DFDCPseudo* pseudopointPtr;
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
  for (unsigned int j = 0; j < n_fdca; j++) {
    pseudopointPtr = (*pseudopointVectorPtr)[j];
    line = pseudopoint2line(*pseudopointPtr);
    docaThis = trajPtr->doca(line, pocaThis);
    distThis = velDrift*(pseudopointPtr->time - pocaThis.getT()/c);
    if (docaThis > distThis) {
      residThis = (distThis - docaThis)/errorFDCA;
    } else {
      residThis = innerResidFrac*(distThis - docaThis)/errorFDCA;
    }
    posWireThis = line.poca(pocaThis);
    if (debug_level > 2) cout << "residFDCCathode: j = " << j
			      << " dist = " << distThis
			      << " doca = " << docaThis
			      << " poca xyzt = " << pocaThis.getX()
			      << ' ' << pocaThis.getY()
			      << ' ' << pocaThis.getZ()
			      << ' ' << pocaThis.getT()/c
			      << " resid = " << residThis << endl;
    errorThis = errorFDCA;
    doca.push_back(docaThis);
    dist.push_back(distThis);
    poca.push_back(pocaThis);
    posWire.push_back(posWireThis);
    resid.push_back(residThis);
    error.push_back(errorThis);
  }
};

DLine residFDCCathode::pseudopoint2line(const DFDCPseudo &pseudopoint) {
  const DFDCWire *wire = pseudopoint.wire;
  double x = wire->origin.X();
  double y = wire->origin.Y();
  double z = wire->origin.Z();
  double theta = acos(wire->udir.z());
  double phi = atan2(wire->udir.y(), wire->udir.x());
  if (debug_level >= 4) cout << setprecision(14) << "residFDCCathode::pseudopoint2line: x = " << x << " y = " << y << " z = " << z << " theta " << theta << " phi " << phi << endl;
  DLine line(x, y, z, theta, phi, debug_level);
  return line;
}

void residFDCCathode::setInnerResidFrac(double innerResidFracIn) {
  innerResidFrac = innerResidFracIn;
}

void residFDCCathode::getResids(vector<double> &residsRef) {
  residsRef = resid;
  return;
}

void residFDCCathode::getDetails(vector<double> &docasRef, vector<double> &distsRef,
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
