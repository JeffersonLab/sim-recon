// $Id$
//
//    File: JEventProcessor_occupancy_online.h
// Created: Tue Apr 12 09:43:54 EDT 2016
// Creator: zihlmann (on Linux gluon47.jlab.org 2.6.32-358.23.2.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_occupancy_online_
#define _JEventProcessor_occupancy_online_

#include <JANA/JEventProcessor.h>



#include <TDirectory.h>
#include <TH2.h>
#include <TH1.h>

using namespace jana;
using namespace std;

class JEventProcessor_occupancy_online:public jana::JEventProcessor{
	public:
		JEventProcessor_occupancy_online();
		~JEventProcessor_occupancy_online();
		const char* className(void){return "JEventProcessor_occupancy_online";}


		//------------------------ BCAL -----------------------
		TH2I *bcal_adc_occ;
		TH2I *bcal_tdc_occ;

		//------------------------ FCAL -----------------------
		TH1I *fcal_num_events;
		TH2I* fcal_occ;

		//------------------------ FDC ------------------------

		//------------------------ PS/PSC ---------------------

		//------------------------ RF -------------------------

		//------------------------ ST -------------------------
		TH1I *st_adc_occ;
		TH1I *st_tdc_occ;

		//------------------------ TAGH -----------------------

		//------------------------ TAGM -----------------------

		//------------------------ TPOL -----------------------

		//------------------------ CDC ------------------------
		TH1I *cdc_num_events;
		TH2F *cdc_occ_ring[28];

		//------------------------ TOF ------------------------
		TH1I *tdcOccS;
		TH1I *tdcOccN;
		TH1I *tdcOccU;
		TH1I *tdcOccD;

		TH1I *adcOccS;
		TH1I *adcOccN;
		TH1I *adcOccU;
		TH1I *adcOccD;


	private:
		jerror_t init(void);
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);
		jerror_t erun(void);
		jerror_t fini(void);
};

#endif // _JEventProcessor_occupancy_online_

