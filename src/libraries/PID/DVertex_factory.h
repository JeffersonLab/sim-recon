// $Id$
//
//    File: DVertex_factory.h
// Created: Tue Apr  6 17:01:54 EDT 2010
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//

#ifndef _DVertex_factory_
#define _DVertex_factory_

#include <JANA/JFactory.h>
#include <DVertex.h>
#include <TRACKING/DHoughFind.h>
#include <PID/DChargedTrack.h>
#include "HDGEOMETRY/DRootGeom.h"
#include <PID/DNeutralShower.h>
#include <ANALYSIS/DAnalysisUtilities.h>
#include <TH1F.h>

class DVertex_factory : public jana::JFactory<DVertex>{
	public:
		DVertex_factory(){};
		~DVertex_factory(){
		  for (unsigned int i=0;i<dVertexInfoPool.size();i++){
		    delete dVertexInfoPool[i];
		  }
		};

	       class vertex_t{
	       public:
		 unsigned int trackbits;
		 DVector3 pos;
		 DVector3 weight;
		 double t0,t0_weight;

		 vertex_t(unsigned int trackbits,DVector3 &pos,DVector3 &weight,double t0,double t0_weight)
		   :trackbits(trackbits),pos(pos),weight(weight),t0(t0),t0_weight(t0_weight){};
	       };

		class vertexInfo_t : public DHoughFind {
			public:
				bool is_in_group;
				bool is_matched_to_vertex;
				const DChargedTrack* dChargedTrack;
				double t;
				double sigmat;
				double z;
				double sigmaz;
    
				void Reset(void){
					ResetHist(); // (from DHoughFind)
					is_in_group = false;
					is_matched_to_vertex=false;
					dChargedTrack = NULL;
				}
		};

		virtual jerror_t MakeVertices(vector<const DChargedTrack*> &locChargedTracks);
		void FillVertexInfoChargedTrack(DVertex_factory::vertexInfo_t *locVertexInfo, const DChargedTrack *locChargedTrack);
		virtual void AssignParticlesToGroups(vector<vertexInfo_t*> &locVertexInfos, vector< vector<vertexInfo_t *> > &locVertexInfoGroups);
		bool AllInGroups(vector<vertexInfo_t*> &locVertexInfos);
		


	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		const DAnalysisUtilities* dAnalysisUtilities;

		vector<vertexInfo_t*> dVertexInfoPool;

		float GROUP_NUM_SIGMAS_TIME;
		float GROUP_NUM_SIGMAS_Z;
		double dTargetCenter_Z;  

		// Pool of memory heavy vertexInfo_t objects
		unsigned int MAX_VERTEXINFOS;
  
		// Values to define the histo limits
		unsigned int Nbinst;
		double tmin;
		double tmax;
		unsigned int Nbinsz;
		double zmin;
		double zmax;

		bool DEBUG_HISTS;
		TH1F *Nsigmas_t_particles;
		TH1F *Nsigmas_z_particles;
};

#endif // _DVertex_factory_

