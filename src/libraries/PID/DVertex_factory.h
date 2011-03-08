// $Id$
//
//    File: DVertex_factory.h
// Created: Tue Apr  6 17:01:54 EDT 2010
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//

#ifndef _DVertex_factory_
#define _DVertex_factory_

#include <JANA/JFactory.h>
#include "DVertex.h"
#include "DParticleID.h"
#include <TRACKING/DHoughFind.h>
#include <PID/DPhoton.h>

#include <TH1.h>
#include <TH2.h>


class DVertex_factory:public jana::JFactory<DVertex>{
 public:
  DVertex_factory(){};
  ~DVertex_factory(){};

  class vertexInfo_t : public DHoughFind {
  public:
    bool is_in_group;
    vector<const DTrackTimeBased *>*hypotheses;
    const DPhoton *photon;
    double t;
    double sigmat;
    double z;
    double sigmaz;
    
    void Reset(void){
      ResetHist(); // (from DHoughFind)
      is_in_group = false;
      hypotheses = NULL;
      photon = NULL;
    }
  };
  jerror_t MakeVertices(vector<vector<const DTrackTimeBased*> >&tracks_by_candidate);
  void FillVertexInfoChargedTrack(DVertex_factory::vertexInfo_t *pi, 
				  vector<const DTrackTimeBased *>*hypotheses);
  bool AllInGroups(vector<vertexInfo_t*> &vertices);
  
 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.
  
  double target_length;
  double target_z;
  
  // Pool of memory heavy partInfo_t objects
  unsigned int MAX_VERTEXINFOS;
  vector<vertexInfo_t*> vertexInfos_pool;
  
  // Values to define the histo limits
  unsigned int Nbinst;
  double tmin;
  double tmax;
  unsigned int Nbinsz;
  double zmin;
  double zmax;

  int eventnumber;
 
  DParticleID *pid_algorithm;

  bool DEBUG_HISTS;
  TH2F *fcal_match;
  TH2F *fcal_dt;
};

#endif // _DVertex_factory_

