// $Id$
//
//    File: JEventProcessor_BCAL_attenlength_gainratio.h
// Created: Mon Aug 10 10:17:48 EDT 2015
// Creator: dalton (on Linux gluon02.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_BCAL_attenlength_gainratio_
#define _JEventProcessor_BCAL_attenlength_gainratio_

#include <JANA/JEventProcessor.h>
#include "TH2.h"

class JEventProcessor_BCAL_attenlength_gainratio:public jana::JEventProcessor{
	public:
		JEventProcessor_BCAL_attenlength_gainratio();
		~JEventProcessor_BCAL_attenlength_gainratio();
		const char* className(void){return "JEventProcessor_BCAL_attenlength_gainratio";}

		static const int nummodule=48;
		static const int numlayer=4;
		static const int numsector=4;

		// Summary histograms
		TH1I *hist_attenlength;
		TH1I *hist_gainratio;
		TH1I *hist_attenlength_err;
		TH1I *hist_gainratio_err;
		TH1I *hist_attenlength_relerr;
		TH1I *hist_gainratio_relerr;
		TH2F *hist2D_peakattenlength;
		TH2F *hist2D_peakgainratio;
		TH2F *hist2D_intattenlength;
		TH2F *hist2D_intgainratio;

		// Channel by channel histograms
		TH2I *logpeakratiovsZ_all;
		TH2I *logintratiovsZ_all;
		TH2I *logpeakratiovsZ[nummodule][numlayer][numsector];
		TH2I *logintratiovsZ[nummodule][numlayer][numsector];
		TH2I *EvsZ[nummodule][numlayer][numsector];

		// Debug histograms to help understand data
		TH2I *EvsZ_all;
		TH2I *EvsZ_layer[4];
		TH2F *hist2D_aveZ;

		double z_target_center;

	private:
		uint32_t VERBOSE;

		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_BCAL_attenlength_gainratio_

