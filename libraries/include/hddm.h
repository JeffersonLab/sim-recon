
#include "DContainer.h"
#include <TVector3.h>
#include <TLorentzVector.h>

#ifndef _HDDM_H_
#define _HDDM_H_

//------------------------------- FCALcluster -------------------------------
typedef struct{ 
	int nhits;				///< number of hits
	int *hit_index;		///< indicies to rows in FCALhits container 
	TVector3 position;	///< reconstructed hit position
	float time;				///< time of hit
	float width;			///< transverse spread of cluster
	float energy;			///< energy of cluster 
}FCALcluster_t; 

class FCALclusters_t: public DContainer 
{ 
	public: 
		FCALclusters_t(void):DContainer((void**)&FCALcluster,sizeof(FCALcluster_t), "FCALclusters"){} 
		FCALcluster_t *FCALcluster; 
}; 

//------------------------------- CDChit -------------------------------
typedef struct{ 
	TVector3 pos;	///< Coordinates of hit
	float t;			///< time of hit
	int  track;		///< track id associated with hit
}CDChit_t; 

class CDChits_t: public DContainer 
{ 
	public: 
		CDChits_t(void):DContainer((void**)&CDChit,sizeof(CDChit_t), "CDChits"){} 
		CDChit_t *CDChit; 
};

//------------------------------- CDCtrack -------------------------------
typedef struct{
	TLorentzVector p;	///< 4-momentum of track
	TVector3 dir;		///< Unit vector in direction of track
	float q;				///< Charge of particle
	int track;			///< track id
	float x0;			///< x-coord. center of track in X/Y plane
	float y0;			///< y-coord. center of track in X/Y plane
}CDCtrack_t;

class CDCtracks_t: public DContainer 
{ 
	public: 
		CDCtracks_t(void):DContainer((void**)&CDCtrack,sizeof(CDCtrack_t), "CDCtracks"){} 
		CDCtrack_t *CDCtrack; 
};

//----------------------------------------------------------------------------
//------------------------------- hddm_containers_t -------------------------------
//----------------------------------------------------------------------------
typedef struct {

		/// This struct should consist ONLY of classes derived from DContainer
		FCALclusters_t 	*FCALclusters;
		CDChits_t 			*CDChits;
		CDCtracks_t			*CDCtracks;
}hddm_containers_t;

// in DANA/hddm_containers.cc
derror_t init_hddm_containers_t(hddm_containers_t *hddm);
derror_t delete_hddm_containers_t(hddm_containers_t *hddm);

#endif // _HDDM_H_
