// $Id: DEventProcessor_candidate_tree.h 5685 2009-11-05 19:57:00Z staylor $
//
//    File: DEventProcessor_candidate_tree.h
// Created: Fri Aug 20 13:22:33 EDT 2010
// Creator: davidl (on Darwin harriet.jlab.org 9.8.0)
//

#ifndef _DEventProcessor_candidate_tree_
#define _DEventProcessor_candidate_tree_

#include <pthread.h>
#include <map>
using std::map;

#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include <JANA/JFactory.h>
#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>

#include <PID/DKinematicData.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <TRACKING/DMCTrackHit.h>
#include <TRACKING/DMCThrown.h>
#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCHit.h>
#include <FDC/DFDCWire.h>
#include <HDGEOMETRY/DLorentzDeflections.h>

#include "dchit.h"
#include "trackpar.h"

class DEventProcessor_candidate_tree:public JEventProcessor{

	public:
		DEventProcessor_candidate_tree();
		~DEventProcessor_candidate_tree();

		TTree *fdchits, *cdchits, *ttrack;
		dchit cdchit, fdchit;
		dchit *cdchit_ptr, *fdchit_ptr;
		trackpar trk;
		trackpar *trk_ptr;
		
		DCoordinateSystem target;
		const DLorentzDeflections *lorentz_def;//< Correction to FDC cathodes due to Lorentz force

		class hit_info_t{
			public:
				// Inputs
				DReferenceTrajectory *rt;
				const DCoordinateSystem *wire;
				double tdrift;
				
				// Outputs
				double doca;
				double dist;
				double tof;
				double u;
				double u_lorentz;
				int LRfit;
				bool LRis_correct;
				TVector3 pos_doca;
				TVector3 mom_doca;
				TVector3 pos_wire;

				void FindLR(vector<const DMCTrackHit*> &mctrackhits, const DLorentzDeflections *lorentz_def=NULL);
		};
		
		
	private:
		jerror_t init(void);	///< Invoked via DEventProcessor virtual method
		jerror_t brun(JEventLoop *loop, int runnumber);
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
		jerror_t erun(void);					///< Invoked via DEventProcessor virtual method
		jerror_t fini(void);					///< Invoked via DEventProcessor virtual method
			
		pthread_mutex_t mutex;
		
		int NLRbad, NLRgood, NLRfit_unknown;
		int Nevents;
		DReferenceTrajectory *rt;
		
		int Nwarnings;
};

#endif // _DEventProcessor_candidate_tree_

