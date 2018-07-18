// $Id$
//
//    File: DDIRCLut.h
//

#ifndef _DDIRCLut_
#define _DDIRCLut_

#include <vector>
#include <utility>
#include <string>
#include <memory>

#include "JANA/JObject.h"

#include "GlueX.h"
#include "DLorentzVector.h"

using namespace std;

class DDIRCLut : public jana::JObject
{
 public:
  JOBJECT_PUBLIC(DDIRCLut);

  vector< pair<double,double> > dPhoton; // thetaC, deltaT

  void toStrings(vector<pair<string,string> > &items) const{
    //AddString(items, "E", "%3.5f", dEnergy);
    //AddString(items, "x", "%3.2f", dSpacetimeVertex.X());
    //AddString(items, "y", "%3.2f", dSpacetimeVertex.Y());
    //AddString(items, "z", "%3.2f", dSpacetimeVertex.Z());
    //AddString(items, "t", "%3.2f", dSpacetimeVertex.T());
  }
};

#endif // _DDIRCLut_

