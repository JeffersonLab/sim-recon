// $Id$
//
//    File: JEventProcessor_CDC_amp.h
// Created: Tue Sep  6 10:13:02 EDT 2016
// Creator: njarvis (on Linux egbert 2.6.32-642.3.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_CDC_amp_
#define _JEventProcessor_CDC_amp_


#include <JANA/JEventProcessor.h>


#include "TRACKING/DTrackFitter.h"
#include "TRACKING/DTrackTimeBased.h"

#include "TRACKING/DMCThrown.h"

#include "CDC/DCDCTrackHit.h"
#include "CDC/DCDCHit.h"
#include "CDC/DCDCDigiHit.h"
#include "DAQ/Df125CDCPulse.h"

#include "DAQ/Df125Config.h"
#include "TRIGGER/DTrigger.h"

#include <stdint.h>
#include <vector>
#include <TMath.h>

#include <TDirectory.h>

#include <TH2.h>
#include <TH1.h>

using namespace std;
using namespace jana;





class JEventProcessor_CDC_amp:public jana::JEventProcessor{
	public:
		JEventProcessor_CDC_amp();
		~JEventProcessor_CDC_amp();
		const char* className(void){return "JEventProcessor_CDC_amp";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

                int32_t run_number;

                // default scaling factors will be overridden by Df125Config if present

                uint16_t ASCALE = 1;   //amplitude  this was 8 for runs before 40,000
                uint16_t PSCALE = 1;   //ped

                // histograms

                TH1I *asum = NULL; 
                TH2I *an = NULL; 

                TH1I *atsum = NULL; 
                TH2I *atn = NULL; 

                TH1I *attsum = NULL; 
                TH2I *attn = NULL; 

                TH2D *atheta = NULL;

                TH1D *qsum = NULL; 
                TH2D *qn = NULL; 

                TH1D *qtsum = NULL; 
                TH2D *qtn = NULL; 

                TH1D *qttsum = NULL; 
                TH2D *qttn = NULL; 

                TH2D *qtheta = NULL;

};

#endif // _JEventProcessor_CDC_amp_

