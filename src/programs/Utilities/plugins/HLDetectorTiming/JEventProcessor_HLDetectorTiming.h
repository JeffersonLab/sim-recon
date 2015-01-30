// $Id$
//
//    File: JEventProcessor_HLDetectorTiming.h
// Created: Mon Jan 12 14:37:56 EST 2015
// Creator: mstaib (on Linux egbert 2.6.32-431.20.3.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_HLDetectorTiming_
#define _JEventProcessor_HLDetectorTiming_

#include <DAQ/DEPICSvalue.h>
#include <JANA/JEventProcessor.h>
#include <BCAL/DBCALHit.h>
#include <CDC/DCDCHit.h>
#include <FCAL/DFCALHit.h>
#include <FDC/DFDCHit.h>
#include <TOF/DTOFHit.h>
#include <START_COUNTER/DSCHit.h>
#include <TAGGER/DTAGHHit.h>
#include <TAGGER/DTAGMHit.h>
#include <PID/DParticleID.h>

#include "TFitResult.h"
#include "TF1.h"
#include "TH1D.h"
#include "TObjArray.h"
#include "TMath.h"

//#include "HistogramTools.h"

//class JEventProcessor_HLDetectorTiming:public jana::JEventProcessor, public HistogramTools{
class JEventProcessor_HLDetectorTiming:public jana::JEventProcessor{
    public:
		JEventProcessor_HLDetectorTiming();
		~JEventProcessor_HLDetectorTiming();
		const char* className(void){return "JEventProcessor_HLDetectorTiming";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

        //HistogramTools *histoTools;
        const DParticleID* dParticleID;

        void DoRoughTiming();
        void DoTDCADCAlign();
        void DoTrackBased();

        int GetCCDBIndexTOF(const DTOFHit *);
        int GetCCDBIndexBCAL(const DBCALHit *); ///< Not implimented

        double BEAM_CURRENT;
        int DO_ROUGH_TIMING, DO_TDC_ADC_ALIGN, DO_TRACK_BASED, DO_VERIFY, REQUIRE_BEAM;

        // The final setp requires some shifts relative to the previous values, need to store them

        double fcal_t_base, bcal_t_base, tof_t_base_fadc, tof_t_base_tdc;
        double tagm_fadc_time_offsets[103], tagm_tdc_time_offsets[103];
        double tagh_fadc_time_offsets[275], tagh_tdc_time_offsets[275];
        vector<double> sc_tdc_time_offsets;
        vector<double> tof_tdc_time_offsets;
};

#endif // _JEventProcessor_HLDetectorTiming_

