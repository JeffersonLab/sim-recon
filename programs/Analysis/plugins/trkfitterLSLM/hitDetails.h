#ifndef _HITDETAILS_H_
#define _HITDETAILS_H_

#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Vector/LorentzVector.h>
using namespace CLHEP;

class FDCHitDetails {
  static int instanceCount;
 public:
  FDCHitDetails();
  ~FDCHitDetails();
  double doca;
  HepVector rCorr;
  HepLorentzVector poca;
  static int getInstanceCount();
}; 

class CDCHitDetails {
  static int instanceCount;
 public:
  CDCHitDetails();
  ~CDCHitDetails();
  double doca;
  double dist;
  HepLorentzVector poca;
  HepVector posWire;
  static int getInstanceCount();
};

#endif /* ifndef _HITDETAILS_H_ */
