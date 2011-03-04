// $Id$
//
//    File: DVertex.h
// Created: Tue Apr  6 17:01:54 EDT 2010
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//

#ifndef _DVertex_
#define _DVertex_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

#include <TRACKING/DTrackTimeBased.h>

class DVertex:public jana::JObject{
 public:
  JOBJECT_PUBLIC(DVertex);
  
  DLorentzVector x; // vertex position in cm + vertex time in ns
  DMatrix cov;	// covariance matrix
  bool beamline_used;
  
  typedef struct{
    const DTrackTimeBased *track;
    double FOM;
    double tprojected;
  }track_info_t;

  vector<vector<track_info_t> >hypotheses;

  // Objects used to calculate this added as Associated Objects
  void toStrings(vector<pair<string,string> > &items)const{
    
    AddString(items, "x", "%3.2f", x.X());
    AddString(items, "y", "%3.2f", x.Y());
    AddString(items, "z", "%3.2f", x.Z());
    AddString(items, "t", "%3.2f", x.T());
    AddString(items, "beamline_used", "%d", beamline_used);
    AddString(items, "Ntracks", "%d", hypotheses.size());
  }
};

#endif // _DVertex_

