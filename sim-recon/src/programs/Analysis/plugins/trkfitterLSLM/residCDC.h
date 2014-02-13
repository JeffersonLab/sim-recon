#ifndef _RESIDCDC_H_
#define _RESIDCDC_H_

#include "CDC/DCDCTrackHit.h"
#include "MyTrajectory.h"
#include "DLine.h"

class residCDC {
 public:
  residCDC(vector<const DCDCTrackHit*> *trackhits, const MyTrajectory *trajectory,
		    int level = 1);
  void calcResids();
  void setInnerResidFrac(double innerResidFracIn);
  void getResids(vector<double> &residsRef);
  void getDetails(vector<double> &docasRef, vector<double> &distsRef,
		  vector<double> &errorsRef,
		  vector<HepLorentzVector> &pocasRef,
		  vector<HepVector> &posWiresRef);
 private:
  unsigned int n_cdc;
  vector<const DCDCTrackHit*> *trkhitVectorPtr;
  const MyTrajectory *trajPtr;
  DLine trackhit2line(const DCDCTrackHit &trackhit);
  int debug_level;
  double innerResidFrac;
  vector<double> doca, dist, resid, error;
  vector<HepLorentzVector> poca;
  vector<HepVector> posWire;
  double errorCDC;
};

#endif // _RESIDCDC_H_

// end of C++ source
