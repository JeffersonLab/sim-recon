// $Id$
//
//    File: JEventProcessor_lowlevel_online.h
// Created: Tue Apr 12 09:43:54 EDT 2016
// Creator: zihlmann (on Linux gluon47.jlab.org 2.6.32-358.23.2.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_lowlevel_online_
#define _JEventProcessor_lowlevel_online_

#include <JANA/JEventProcessor.h>

#include <GlueX.h>
#include <PAIR_SPECTROMETER/DPSGeometry.h>

#include <TDirectory.h>
#include <TH2.h>
#include <TH1.h>

#include <vector>

using namespace jana;
using namespace std;

class JEventProcessor_lowlevel_online:public jana::JEventProcessor{
	public:
		JEventProcessor_lowlevel_online();
		~JEventProcessor_lowlevel_online();
		const char* className(void){return "JEventProcessor_lowlevel_online";}

        bool INDIVIDUAL_CHANNEL_DATA;
        vector<int> Nstraws_integrated;

		//------------------------ BCAL -----------------------
        TH1I *bcal_num_events;

        TH1I *bcal_adc_integral;
        TH1I *bcal_adc_integral_pedsub;
        TH1I *bcal_adc_peak;
        TH1I *bcal_adc_peak_pedsub;
        TH1I *bcal_adc_time;
        TH1I *bcal_adc_pedestal;
        TH1I *bcal_adc_quality;
        TH2I *bcal_adc_integral_chan;
        TH2I *bcal_adc_integral_pedsub_chan;
        TH2I *bcal_adc_peak_chan;
        TH2I *bcal_adc_peak_pedsub_chan;
        TH2I *bcal_adc_time_chan;
        TH2I *bcal_adc_pedestal_chan;
        TH2I *bcal_adc_quality_chan;

		//------------------------ CDC ------------------------
		TH1I *cdc_num_events;
		TH2F *cdc_occ_ring[28];

		//------------------------ FCAL -----------------------
		TH1I *fcal_num_events;
		TH2F* fcal_occ;

		//------------------------ FDC ------------------------
        TH1I *fdc_num_events;
        TH2F *fdc_cathode_occ;
        TH2F *fdc_wire_occ;

		//------------------------ PS/PSC ---------------------
		TH1I *ps_num_events;
		TH1I *psc_adc_left_occ;
		TH1I *psc_adc_right_occ;
		TH1I *psc_tdc_left_occ;
		TH1I *psc_tdc_right_occ;
		TH1I *ps_left_occ;
		TH1I *ps_right_occ;

		//------------------------ RF -------------------------
		TH1I *rf_num_events;
		TH1D* rf_occ; //TH1D ON PURPOSE!
		map<DetectorSystem_t, double> dRFBinValueMap;

		//------------------------ ST -------------------------
		TH1I *st_num_events;
		TH1I *st_adc_occ;
		TH1I *st_tdc_occ;

		//------------------------ TAGH -----------------------
		TH1I *tag_num_events;
		TH1I *tagh_adc_occ;
		TH1I *tagh_tdc_occ;

		//------------------------ TAGM -----------------------
		TH1I *tagm_adc_occ;
		TH1I *tagm_tdc_occ;

		//------------------------ TPOL -----------------------
		TH1I *tpol_occ;

		//------------------------ TOF ------------------------
		TH1I *tof_num_events;
		TH1I *tof_tdc_S_occ;
		TH1I *tof_tdc_N_occ;
		TH1I *tof_tdc_U_occ;
		TH1I *tof_tdc_D_occ;

		TH1I *tof_adc_S_occ;
		TH1I *tof_adc_N_occ;
		TH1I *tof_adc_U_occ;
		TH1I *tof_adc_D_occ;


	private:
		jerror_t init(void);
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);
		jerror_t erun(void);
		jerror_t fini(void);
};

#endif // _JEventProcessor_lowlevel_online_

