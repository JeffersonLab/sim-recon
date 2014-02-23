#ifndef _RESIDCDC_H_
#define _RESIDCDC_H_

#include "FDC/DFDCPseudo.h"
#include "MyTrajectory.h"
#include "DLine.h"

class residFDCCathode {
 public:
  residFDCCathode(vector<const DFDCPseudo*> *pseudopoints, const MyTrajectory *trajectory,
		    int level = 1);
  void calcResids();
  void setInnerResidFrac(double innerResidFracIn);
  void getResids(vector<double> &residsRef);
  void getDetails(vector<double> &docasRef, vector<double> &distsRef,
		  vector<double> &errorsRef,
		  vector<HepLorentzVector> &pocasRef,
		  vector<HepVector> &posWiresRef);
 private:
  unsigned int n_fdca;
  vector<const DFDCPseudo*> *pseudopointVectorPtr;
  const MyTrajectory *trajPtr;
  DLine pseudopoint2line(const DFDCPseudo &pseudopoint);
  int debug_level;
  double innerResidFrac;
  vector<double> doca, dist, resid, error;
  vector<HepLorentzVector> poca;
  vector<HepVector> posWire;
  double errorFDCA;
};

#endif // _RESIDCDC_H_

// end of C++ source
