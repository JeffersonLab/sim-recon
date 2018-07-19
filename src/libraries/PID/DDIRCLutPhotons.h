// $Id$
//
//    File: DDIRCLutPhotons.h
//

#ifndef _DDIRCLutPhotons_
#define _DDIRCLutPhotons_

#include <vector>
#include <utility>
#include <string>
#include <memory>

#include "JANA/JObject.h"

#include "GlueX.h"
#include "DLorentzVector.h"

using namespace std;

class DDIRCLutPhotons : public jana::JObject
{
 public:
  JOBJECT_PUBLIC(DDIRCLutPhotons);

  vector< pair<double,double> > dPhoton; // thetaC, deltaT

  void toStrings(vector<pair<string,string> > &items) const{
    //AddString(items, "E", "%3.5f", dEnergy);
    //AddString(items, "x", "%3.2f", dSpacetimeVertex.X());
    //AddString(items, "y", "%3.2f", dSpacetimeVertex.Y());
    //AddString(items, "z", "%3.2f", dSpacetimeVertex.Z());
    //AddString(items, "t", "%3.2f", dSpacetimeVertex.T());
  }
};

#endif // _DDIRCLutPhotons_

