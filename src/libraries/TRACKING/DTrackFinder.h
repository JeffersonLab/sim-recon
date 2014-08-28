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
#include "FDC/DFDCPseudo.h"
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

  class fdc_hit_t{
  public:
  fdc_hit_t(const DFDCPseudo *hit=NULL,bool used=false):
    hit(hit),used(used){}
    const DFDCPseudo *hit;
    bool used;
  };
  
  class fdc_segment_t{
  public:
    fdc_segment_t(vector<const DFDCPseudo *>&input_hits,bool matched=false){
      for (unsigned int i=0;i<input_hits.size();i++){
	this->hits.push_back(input_hits[i]);
      }
      this->S=FindStateVector();
      this->matched=matched;
    };
    ~fdc_segment_t(){};

    DMatrix4x1 FindStateVector(void) const;

    bool matched;
    DMatrix4x1 S;
    vector<const DFDCPseudo *>hits;
    
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
  void AddHit(const DFDCPseudo *hit);
  bool FindFDCSegments(void);
  bool LinkFDCSegments(void);
  bool FindAxialSegments(void);
  bool LinkCDCSegments(void);
  bool MatchCDCHit(const DVector3 &vhat,const DVector3 &pos0,
		   const DCDCTrackHit *hit);
  const vector<cdc_track_t>&GetCDCTracks(void) const {return cdc_tracks;};
  const vector<fdc_segment_t>&GetFDCTracks(void) const {return fdc_tracks;};
  
  double FindDoca(double z,const DMatrix4x1 &S,const DVector3 &wdir,
		  const DVector3 &origin,DVector3 *poca=NULL) const;

 private:
  // Prohibit default constructor
  DTrackFinder();
   
  vector<cdc_hit_t>axial_hits;
  vector<cdc_hit_t>stereo_hits;
  vector<cdc_segment_t>axial_segments;
  vector<cdc_track_t>cdc_tracks;

  vector<fdc_hit_t>fdc_hits;
  vector<fdc_segment_t>fdc_segments[4];
  vector<fdc_segment_t>fdc_tracks;
};

#endif // _DTrackFinder_

