#ifndef _COMBINEDRESIDFUNC_H_
#define _COMBINEDRESIDFUNC_H_

#define PIOVER2 1.570796327
#define PIOVER4 0.785398163
#include "FDC/DFDCPseudo.h"
#include "CDC/DCDCTrackHit.h"
#include "hitDetails.h"
#include "MyTrajectory.h"

#define BIG_DOUBLE 1.0e12
#define DRIFT_VELOCITY 55e-4
#define ERROR_FDC 0.0250
#define ERROR_CDC 0.0180

class combinedResidFunc:public residFunc {
 public:
  combinedResidFunc(vector<const DFDCPseudo*> *pseudopoints,
		    vector<const DCDCTrackHit*> *trackhits,
		    MyTrajectory *trajectory, const DLorentzDeflections *lorentz_def,
		    int level = 1);
  void resid(const HepVector *x, void *data, HepVector *f);
  void deriv(const HepVector *x, void *data, HepMatrix *J);
  void residAndDeriv(const HepVector *x, void *data, HepVector *f,
		     HepMatrix *J);
  inline unsigned int getN() {return n_fdc + n_cdc;};
  inline unsigned int getP() {return trajPtr->getNumberOfParams();};
  void setStoreDetails(bool storeDetailsValue);
  void clearDetails();
  inline vector<CDCHitDetails*> *getCDCDetails(){return &CDCDetails;};
  inline vector<FDCHitDetails*> *getFDCDetails(){return &FDCDetails;};
  void setInnerResidFrac(double innerResidFracIn);
 private:
  unsigned int n_fdc, n_cdc;
  vector<const DFDCPseudo*> *ppPtr;
  vector<const DCDCTrackHit*> *trkhitPtr;
  MyTrajectory *trajPtr;
  HepVector pseudo2HepVector(const DFDCPseudo &pseudopoint);
  DLine trackhit2line(const DCDCTrackHit &trackhit);
  vector<double> delta;
  int debug_level;
  bool getCorrectionSign(const DFDCPseudo &pseudopoint, double x, double y, double deltaX, double deltaY);
  void getCorrectionValue(const DFDCPseudo &pseudopoint, double x, double y, double z, double t, double &delta_x, double &delta_y);
  const DLorentzDeflections *lorentz_def;
  bool storeDetails;
  FDCHitDetails getDetails(const DFDCPseudo *ppoint, HepVector point);
  vector<FDCHitDetails*> FDCDetails;
  CDCHitDetails getDetails(const DCDCTrackHit *trackhit, DLine line);
  vector<CDCHitDetails*> CDCDetails;
  double innerResidFrac;
};

#endif // _COMBINEDRESIDFUNC_H_

// end of C++ source
