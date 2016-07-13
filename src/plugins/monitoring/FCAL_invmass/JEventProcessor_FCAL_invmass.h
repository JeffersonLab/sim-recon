// $Id$
//
//    File: JEventProcessor_FCAL_invmass.h
// Created: Tue May 31 09:44:35 EDT 2016
// Creator: adesh (on Linux ifarm1101 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _DEventProcessor_FCAL_Shower_
#define _DEventProcessor_FCAL_Shower_

#include <JANA/JEventProcessor.h>
#include <JANA/JApplication.h>

#include <ANALYSIS/DEventWriterROOT.h>
#include <HDDM/DEventWriterREST.h>
#include <ANALYSIS/DHistogramActions.h>
#include "ANALYSIS/DAnalysisUtilities.h"
#include "TRACKING/DTrackFinder.h"

#include "DLorentzVector.h"
#include "TMatrixD.h"


using namespace jana;
using namespace std;

class JEventProcessor_FCAL_invmass : public jana::JEventProcessor
{
	public:
		JEventProcessor_FCAL_invmass(){};
		~JEventProcessor_FCAL_invmass(){};
		const char* className(void){return "JEventProcessor_FCAL_invmass";}
		//DVector3 Calc_CrudeVertex(const deque< const DKinematicData* > & locParticles) const;
		
	       	
	private:
		//const DAnalysisUtilities* dAnalysisUtilities;
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop* locEventLoop, int32_t locRunNumber);	///< Called every time a new run number is detected.
		jerror_t evnt(jana::JEventLoop* locEventLoop, uint64_t locEventNumber);	///< Called every event.
		jerror_t erun(void);						///< Called every time run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed
		//jerror_t fillHists();
//double m_x;

//vector<vector<string> > ParseTSV(const char* s);
	 int XYtoAbsNum(int my_x, int my_y);
	 pair<int,int> AbsNumtoXY(int channel);

        DFCALGeometry *m_fcalgeom;
        DFCALGeometry* mygeom;




	// bool read_gains;

        double z_diff;

       
        TMatrixD m_nhits;
       
		
};

#endif 
