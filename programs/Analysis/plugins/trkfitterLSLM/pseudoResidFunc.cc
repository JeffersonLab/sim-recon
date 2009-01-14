#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/Vector.h>

using namespace std;

#include "FDC/DFDCPseudo_factory.h"
#include "MyTrajectory.h"
#include "residFunc.h"
#include "pseudoResidFunc.h"

// constructor, takes pointer to vector of pseudopoints and a pointer
// to a trajectory
pseudoResidFunc::pseudoResidFunc(vector<const DFDCPseudo*> *pseudopoints,
				 MyTrajectory *trajectory) : 
  n(pseudopoints->size()), ppPtr(pseudopoints), trajPtr(trajectory),
  delta(trajPtr->getDelta())
{
  return;
}

void pseudoResidFunc::resid(const HepVector *x, void *data, HepVector *f){
// input parameters
//   x: pointer to a vector of fit parameters
//   data: ???
// output parameters
//   f: pointer to vector of residuals
  //  cout << "resid called\n";
  //  cout << "params: " << *x;
  // do a swim with the input parameters
  trajPtr->swim(*x);
  // populate f vector with residuals
  HepVector point(3);
  HepLorentzVector poca;
  for (unsigned int i = 0; i < n; i++) {
    point(1) = ((*ppPtr)[i])->x;
    point(2) = ((*ppPtr)[i])->y;
    point(3) = ((*ppPtr)[i])->wire->origin(2);
    (*f)(i + 1) = trajPtr->doca(point, poca);
  }
  //  cout << "resids:" << *f;
  trajPtr->clear();
};

void pseudoResidFunc::deriv(const HepVector *x, void *data, HepMatrix *J){
  //  cout << "deriv called\n";
  //  cout << "params: " << *x << endl;
  HepLorentzVector poca;
  // save base parameters
  HepVector xBase = *x;
  // do central swim
  trajPtr->swim(xBase);
  // store pseudo points as three vectors
  vector<HepVector *>pPoints;
  HepVector *thisPointPtr;
  for (unsigned int i = 0; i < n; i++) {
    thisPointPtr = new HepVector(3);
    (*thisPointPtr)(1) = ((*ppPtr)[i])->x;
    (*thisPointPtr)(2) = ((*ppPtr)[i])->y;
    (*thisPointPtr)(3) = ((*ppPtr)[i])->wire->origin(2);
    pPoints.push_back(thisPointPtr);
    //    cout << "stored pp " << i << *(pPoints[i]);
  }
  // store base residuals
  HepVector docaBase(n);
  for (unsigned int i = 0; i < n; i++) {
    //    cout << "storing base resid for pp "<< i << *(pPoints[i]);
    docaBase(i + 1) = trajPtr->doca(*(pPoints[i]), poca);
  }
  //  cout << "base resids:" << docaBase;
  trajPtr->clear();
  // calculate Jacobian
  unsigned int p = trajPtr->getNumberOfParams();
  //  cout << "pseudoResidFunc::deriv p = " << p << endl;
  HepVector xThis(p);
  double docaThis;
  int iHep, jHep; // index for HepVector () notation
  for (unsigned int i = 0; i < p; i++) {
    iHep = i + 1;
    xThis = xBase; // set params back to base
    xThis(iHep) = xBase(iHep) + delta[i];
    //    cout << "perturbed params: iHep = " << iHep << " delta = " << delta[i] << " values = " << xThis << endl;
    // do the perturbed swim
    trajPtr->swim(xThis);
    for (unsigned int j = 0; j < n; j++) {
      jHep = j + 1;
      docaThis = trajPtr->doca(*(pPoints[j]), poca);
      //      cout << "resid " << j << " = " << docaThis << endl;
      (*J)(jHep, iHep) = (docaThis - docaBase(jHep))/delta[i];
    }
    trajPtr->clear();
  }
  //  cout << "Jacobian:" << *J;
  for (unsigned int i = 0; i < n; i++) {
    delete pPoints[i];
  }
};

void pseudoResidFunc::residAndDeriv(const HepVector *x, void *data, HepVector *f,
		   HepMatrix *J){};
