// $Id$
//
//    File: DTrackFinder.h
// Created: Fri Aug 15 09:43:08 EDT 2014
// Creator: staylor (on Linux ifarm1102 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#ifndef _DTrackFinder_
#define _DTrackFinder_

#include <JANA/jerror.h>
#include "DVector3.h"
#include "CDC/DCDCTrackHit.h"
#include "CDC/DCDCWire.h"
#include "DMatrixSIMD.h"

#include <vector>

class DTrackFinder:public jana::JObject{
 public:
  JOBJECT_PUBLIC(DTrackFinder);
  
  DTrackFinder(JEventLoop *loop);	// require JEventLoop in constructor
  virtual ~DTrackFinder();

  enum state_vector{
    state_x,
    state_y,
    state_tx,
    state_ty,
  };

  class cdc_hit_t{
  public:
  cdc_hit_t(const DCDCTrackHit *hit=NULL,bool used=false):
    hit(hit),used(used){}
    const DCDCTrackHit *hit;
    bool used;
  };
  
  class cdc_segment_t{
  public:
    cdc_segment_t(vector<const DCDCTrackHit *>&input_hits,const DVector3 &dir,
		  bool matched=false){
      for (unsigned int i=0;i<input_hits.size();i++){
	this->hits.push_back(input_hits[i]);
      }
      this->dir=dir;
      this->matched=matched;
    };
    ~cdc_segment_t(){};

    bool matched;
    DVector3 dir;
    vector<const DCDCTrackHit *>hits;
    
  };

  class cdc_track_t{
  public:
    cdc_track_t(vector<const DCDCTrackHit *>myhits){
      for (unsigned int i=0;i<myhits.size();i++){
	if (myhits[i]->is_stereo)this->stereo_hits.push_back(myhits[i]);
	else this->axial_hits.push_back(myhits[i]);
      }
    };  
    ~cdc_track_t(){};

    jerror_t FindStateVector(void);

    vector<const DCDCTrackHit *>axial_hits; 
    vector<const DCDCTrackHit *>stereo_hits;
    DVector3 dir;
    double z; // z-position at which S is reported
    DMatrix4x1 S;

  };

  
  void Reset(void);
  void AddHit(const DCDCTrackHit *hit);
  bool FindAxialSegments(void);
  void LinkCDCSegments(void);
  bool MatchCDCHit(const DVector3 &vhat,const DVector3 &pos0,
		   const DCDCTrackHit *hit);
  const vector<cdc_track_t>&GetTracks(void) const {return cdc_tracks;};
  
  double FindDoca(double z,const DMatrix4x1 &S,const DVector3 &wdir,
		  const DVector3 &origin) const;

 private:
  // Prohibit default constructor
  DTrackFinder();
   
  vector<cdc_hit_t>axial_hits;
  vector<cdc_hit_t>stereo_hits;
  vector<cdc_segment_t>axial_segments;
  vector<cdc_track_t>cdc_tracks;
};

#endif // _DTrackFinder_

