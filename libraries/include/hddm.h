
#include <DBank.h>
#include <TVector3.h>

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
	float time;		///< time of hit
	int  track;		///< transverse spread of cluster
	float energy;	///< energy of cluster 
}CDChit_t; 

class CDChits_t: public DBank 
{ 
	public: 
		CDChits_t(void):DBank((void**)&CDChit,sizeof(CDChit_t), "CDChits"){} 
		CDChit_t *CDChit; 
}; 


//----------------------------------------------------------------------------
//------------------------------- hddm_banks_t -------------------------------
//----------------------------------------------------------------------------
class  hddm_banks_t
{
	public:
		FCALclusters_t 	*FCALclusters;
		CDChits_t 			*CDChits;

		hddm_banks_t(){
			FCALclusters	= new FCALclusters_t();
			CDChits			= new CDChits_t();
		}
		
		~hddm_banks_t(){
			delete FCALclusters;
			delete CDChits;
		}
		
};

#endif // _HDDM_H_
