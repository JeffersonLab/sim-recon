// $Id$
//
//    File: DEventProcessor_BCAL_gainmatrix.h
// Created: Fri Oct 10 16:41:18 EDT 2014
// Creator: wmcginle (on Linux ifarm1101 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#ifndef _DEventProcessor_BCAL_gainmatrix_
#define _DEventProcessor_BCAL_gainmatrix_

#include <JANA/JEventProcessor.h>
#include <JANA/JApplication.h>

#include <ANALYSIS/DEventWriterROOT.h>
#include <HDDM/DEventWriterREST.h>
#include <ANALYSIS/DHistogramActions.h>
#include "ANALYSIS/DAnalysisUtilities.h"
//#include "TRACKING/DTrackFinder.h"

#include "DLorentzVector.h"
#include "TMatrixD.h"


using namespace jana;
using namespace std;

class DEventProcessor_BCAL_gainmatrix : public jana::JEventProcessor
{
	public:
		DEventProcessor_BCAL_gainmatrix(){};
		~DEventProcessor_BCAL_gainmatrix(){};
		const char* className(void){return "DEventProcessor_BCAL_gainmatrix";}
		TTree *BCAL_Neutrals;
	       	uint32_t eventnum  ;
		Float_t E1  ;
		Float_t E1_raw  ;
		Float_t E2  ;
		Float_t E2_raw  ;
		Float_t t1  ;
		Float_t t2  ;
		Float_t z1  ;
		Float_t z2  ;
		Float_t x1  ;
		Float_t x2  ;
		Float_t y1  ;
		Float_t y2  ;
		Float_t psi  ;
		Float_t vertexz ;
		Float_t vertexZ  ;
		Float_t vertexX  ;
		Float_t vertexY  ;
		uint32_t Run_Number  ;
		Int_t num_tracks ;
		Int_t num_showers ;
		Float_t inv_mass  ;
		Float_t inv_mass_raw  ;
		uint32_t logical1;
		uint32_t logical2;
		vector<double> point_energy;
		vector<double> point1_energy_calib;
		vector<double> point2_energy_calib;
		vector<double> point1_channel;
		vector<double> point2_channel;
		vector<pair < double , int > > frac_en;

	private:
		const DAnalysisUtilities* dAnalysisUtilities;
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop* locEventLoop, int32_t locRunNumber);	///< Called every time a new run number is detected.
		jerror_t evnt(jana::JEventLoop* locEventLoop, uint64_t locEventNumber);	///< Called every event.
		jerror_t erun(void);						///< Called every time run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		double m_massbias;

		TMatrixD m_mC;
		TMatrixD m_mD;
		TMatrixD m_mL;
		TMatrix m_nhits;	
		TH2F* h2D_mC;
		TH1F* h1D_mL;
		TH1F* h1D_mD;
		TH1F* h1D_massbias;
		TH1F* h1D_nhits;

		const DEventWriterROOT* dEventWriterROOT;
		const DEventWriterREST* dEventWriterREST;
};

#endif // _DEventProcessor_BCAL_gainmatrix_

