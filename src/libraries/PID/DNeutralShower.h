// $Id$
//
//    File: DNeutralShower.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DNeutralShower_
#define _DNeutralShower_

#include <vector>
#include <utility>
#include <string>
#include <memory>

#include "JANA/JObject.h"

#include "GlueX.h"
#include "DLorentzVector.h"
#include "TMatrixFSym.h"

using namespace std;

class DNeutralShower : public jana::JObject
{
 public:
  JOBJECT_PUBLIC(DNeutralShower);

  oid_t dShowerID;
  DetectorSystem_t dDetectorSystem;

  DLorentzVector dSpacetimeVertex;
  double dEnergy;
  shared_ptr<TMatrixFSym> dCovarianceMatrix; //E, x, y, z, t

  // Quality is a number between 0 and 1 that is closer to
  // 1 for more photon-like showers;  currently only implemented
  // in the FCAL.  Quality = 1 for all BCAL showers.

  double dQuality;

  const JObject* dBCALFCALShower; //is either DBCALShower or DFCALShower: dynamic_cast as appropriate (based on dDetectorSystem)

  void toStrings(vector<pair<string,string> > &items) const{
    AddString(items, "E", "%3.5f", dEnergy);
    AddString(items, "x", "%3.2f", dSpacetimeVertex.X());
    AddString(items, "y", "%3.2f", dSpacetimeVertex.Y());
    AddString(items, "z", "%3.2f", dSpacetimeVertex.Z());
    AddString(items, "t", "%3.2f", dSpacetimeVertex.T());
  }
};

#endif // _DNeutralShower_

