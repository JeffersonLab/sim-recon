// $Id$
//
//    File: DEventProcessor_fdc_covariance_hists.h
// Created: Mon Apr 20 10:18:30 EDT 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.6.0 Darwin Kernel Version 9.6.0)
//

#ifndef _DEventProcessor_fdc_covariance_hists_
#define _DEventProcessor_fdc_covariance_hists_

#include <pthread.h>
#include <map>
using std::map;

#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile2D.h>

#include <JANA/JFactory.h>
#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>

#include <PID/DKinematicData.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <TRACKING/DMCTrackHit.h>
#include <TRACKING/DTrack.h>
#include <TRACKING/DMCThrown.h>
#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCHit.h>
#include <FDC/DFDCWire.h>
#include <HDGEOMETRY/DLorentzDeflections.h>
#include <HDGEOMETRY/DMagneticFieldMap.h>

#include "dchit.h"
class DEventProcessor_fdc_covariance_hists:public JEventProcessor{

	public:
		DEventProcessor_fdc_covariance_hists();
		~DEventProcessor_fdc_covariance_hists();

		TTree *fdchits;
		dchit fdchit;
		dchit *fdchit_ptr;
		//TH2D *fdc_cov_numerator;
		//TH2D *fdc_cov_denominator;
		TProfile2D *fdc_cov;
		TProfile2D *fdc_cath_cov;
		
		const DMagneticFieldMap *bfield;
		const DLorentzDeflections *lorentz_def;//< Correction to FDC cathodes due to Lorentz force
		DReferenceTrajectory *rt;

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
				DVector3 pos_doca;
				DVector3 mom_doca;
				DVector3 pos_wire;

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
		double Z_fdc1;
};

#endif // _DEventProcessor_fdc_covariance_hists_

