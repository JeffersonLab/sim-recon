
#include "DBank.h"
#include <TVector3.h>
#include <TLorentzVector.h>

#ifndef _HDDM_H_
#define _HDDM_H_

//------------------------------- FCALcluster -------------------------------
typedef struct{ 
	int nhits;				///< number of hits
	int *hit_index;		///< indicies to rows in FCALhits bank 
	TVector3 position;	///< reconstructed hit position
	float time;				///< time of hit
	float width;			///< transverse spread of cluster
	float energy;			///< energy of cluster 
}FCALcluster_t; 

class FCALclusters_t: public DBank 
{ 
	public: 
		FCALclusters_t(void):DBank((void**)&FCALcluster,sizeof(FCALcluster_t), "FCALclusters"){} 
		FCALcluster_t *FCALcluster; 
}; 

//------------------------------- CDChit -------------------------------
typedef struct{ 
	TVector3 pos;	///< Coordinates of hit
	float t;			///< time of hit
	int  track;		///< track id associated with hit
}CDChit_t; 

class CDChits_t: public DBank 
{ 
	public: 
		CDChits_t(void):DBank((void**)&CDChit,sizeof(CDChit_t), "CDChits"){} 
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

class CDCtracks_t: public DBank 
{ 
	public: 
		CDCtracks_t(void):DBank((void**)&CDCtrack,sizeof(CDCtrack_t), "CDCtracks"){} 
		CDCtrack_t *CDCtrack; 
};

//----------------------------------------------------------------------------
//------------------------------- hddm_banks_t -------------------------------
//----------------------------------------------------------------------------
typedef struct {

		/// This struct should consist ONLY of classes derived from DBank
		FCALclusters_t 	*FCALclusters;
		CDChits_t 			*CDChits;
		CDCtracks_t			*CDCtracks;
}hddm_banks_t;

// in DANA/hddm_banks.cc
derror_t init_hddm_banks_t(hddm_banks_t *hddm);
derror_t delete_hddm_banks_t(hddm_banks_t *hddm);

#endif // _HDDM_H_
