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
#include <BCAL/DBCALShower.h>
//#include <FCAL/DFCALCluster.h>
#include <FCAL/DFCALShower.h>
#include <DMatrix.h>

class DVertex:public jana::JObject{
 public:
  JOBJECT_PUBLIC(DVertex);
  
  DLorentzVector x; // vertex position in cm + vertex time in ns
  DMatrix cov;	// covariance matrix
  bool beamline_used;
  double t_sigma;
  
  typedef struct{
    const DTrackTimeBased *track;
    double FOM;
    double tprojected;
  }track_info_t;

  vector<vector<track_info_t> >hypotheses;

  class shower_info_t:public DKinematicData{
  public:
    shower_info_t(const DBCALShower *bcal,const DFCALShower *fcal){
      this->bcal=bcal;
      this->fcal=fcal;
      this->matched_track=NULL;
    }
    shower_info_t(const DBCALShower *bcal,const DTrackTimeBased *track){
      this->bcal=bcal;
      this->matched_track=track;
      this->fcal=NULL;
    }
    shower_info_t(const DFCALShower *fcal,const DTrackTimeBased *track){
      this->fcal=fcal;
      this->matched_track=track;
      this->bcal=NULL;      
    }
    const DBCALShower *bcal;
    const DFCALShower *fcal;
    const DTrackTimeBased *matched_track;
  };
  vector<shower_info_t>showers;

  // Objects used to calculate this added as Associated Objects
  void toStrings(vector<pair<string,string> > &items)const{
    
    AddString(items, "x", "%3.2f", x.X());
    AddString(items, "y", "%3.2f", x.Y());
    AddString(items, "z", "%3.2f", x.Z());
    AddString(items, "t", "%3.2f", x.T());
    AddString(items, "beamline_used", "%d", beamline_used);
    AddString(items, "Ntracks", "%d", hypotheses.size());
    AddString(items, "Nphotons","%d", showers.size());
  }
};

#endif // _DVertex_

