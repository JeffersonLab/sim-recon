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
        bool CHECK_EMULATED_DATA;
        bool ANALYZE_F125_DATA;
        bool ANALYZE_F250_DATA;
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

        TH1I *bcal_adc_emudelta_integral;
        TH1I *bcal_adc_emudelta_peak;
        TH1I *bcal_adc_emudelta_pedestal;
        TH1I *bcal_adc_emudelta_coarsetime;
        TH1I *bcal_adc_emudelta_finetime;


		//------------------------ CDC ------------------------
        TH1I *cdc_num_events;
        TH1I *cdc_adc_integral;
        TH1I *cdc_adc_integral_pedsub;
        TH1I *cdc_adc_peak;
        TH1I *cdc_adc_peak_pedsub;
        TH1I *cdc_adc_time;
        TH1I *cdc_adc_pedestal;
        TH1I *cdc_adc_quality;
        TH2I *cdc_adc_integral_chan;
        TH2I *cdc_adc_integral_pedsub_chan;
        TH2I *cdc_adc_peak_chan;
        TH2I *cdc_adc_peak_pedsub_chan;
        TH2I *cdc_adc_time_chan;
        TH2I *cdc_adc_pedestal_chan;
        TH2I *cdc_adc_quality_chan;


		//------------------------ FCAL -----------------------
        TH1I *fcal_num_events;
        TH1I *fcal_adc_integral;
        TH1I *fcal_adc_integral_pedsub;
        TH1I *fcal_adc_peak;
        TH1I *fcal_adc_peak_pedsub;
        TH1I *fcal_adc_time;
        TH1I *fcal_adc_pedestal;
        TH1I *fcal_adc_quality;
        TH2I *fcal_adc_integral_chan;
        TH2I *fcal_adc_integral_pedsub_chan;
        TH2I *fcal_adc_peak_chan;
        TH2I *fcal_adc_peak_pedsub_chan;
        TH2I *fcal_adc_time_chan;
        TH2I *fcal_adc_pedestal_chan;
        TH2I *fcal_adc_quality_chan;

        TH1I *fcal_adc_emudelta_integral;
        TH1I *fcal_adc_emudelta_peak;
        TH1I *fcal_adc_emudelta_pedestal;
        TH1I *fcal_adc_emudelta_coarsetime;
        TH1I *fcal_adc_emudelta_finetime;


		//------------------------ FDC ------------------------
        TH1I *fdc_num_events;
        TH1I *fdc_adc_integral;
        TH1I *fdc_adc_integral_pedsub;
        TH1I *fdc_adc_peak;
        TH1I *fdc_adc_peak_pedsub;
        TH1I *fdc_adc_time;
        TH1I *fdc_adc_pedestal;
        TH1I *fdc_adc_quality;
        TH2I *fdc_adc_integral_chan;
        TH2I *fdc_adc_integral_pedsub_chan;
        TH2I *fdc_adc_peak_chan;
        TH2I *fdc_adc_peak_pedsub_chan;
        TH2I *fdc_adc_time_chan;
        TH2I *fdc_adc_pedestal_chan;
        TH2I *fdc_adc_quality_chan;


		//------------------------ PSC ---------------------
        TH1I *psc_adc_integral;
        TH1I *psc_adc_integral_pedsub;
        TH1I *psc_adc_peak;
        TH1I *psc_adc_peak_pedsub;
        TH1I *psc_adc_time;
        TH1I *psc_adc_pedestal;
        TH1I *psc_adc_quality;
        TH2I *psc_adc_integral_chan;
        TH2I *psc_adc_integral_pedsub_chan;
        TH2I *psc_adc_peak_chan;
        TH2I *psc_adc_peak_pedsub_chan;
        TH2I *psc_adc_time_chan;
        TH2I *psc_adc_pedestal_chan;
        TH2I *psc_adc_quality_chan;

        TH1I *psc_adc_emudelta_integral;
        TH1I *psc_adc_emudelta_peak;
        TH1I *psc_adc_emudelta_pedestal;
        TH1I *psc_adc_emudelta_coarsetime;
        TH1I *psc_adc_emudelta_finetime;


		//------------------------ PS ---------------------
        TH1I *ps_num_events;

        TH1I *ps_adc_integral;
        TH1I *ps_adc_integral_pedsub;
        TH1I *ps_adc_peak;
        TH1I *ps_adc_peak_pedsub;
        TH1I *ps_adc_time;
        TH1I *ps_adc_pedestal;
        TH1I *ps_adc_quality;
        TH2I *ps_adc_integral_chan;
        TH2I *ps_adc_integral_pedsub_chan;
        TH2I *ps_adc_peak_chan;
        TH2I *ps_adc_peak_pedsub_chan;
        TH2I *ps_adc_time_chan;
        TH2I *ps_adc_pedestal_chan;
        TH2I *ps_adc_quality_chan;

        TH1I *ps_adc_emudelta_integral;
        TH1I *ps_adc_emudelta_peak;
        TH1I *ps_adc_emudelta_pedestal;
        TH1I *ps_adc_emudelta_coarsetime;
        TH1I *ps_adc_emudelta_finetime;


		//------------------------ ST -------------------------
        TH1I *st_num_events;
        TH1I *st_adc_integral;
        TH1I *st_adc_integral_pedsub;
        TH1I *st_adc_peak;
        TH1I *st_adc_peak_pedsub;
        TH1I *st_adc_time;
        TH1I *st_adc_pedestal;
        TH1I *st_adc_quality;
        TH2I *st_adc_integral_chan;
        TH2I *st_adc_integral_pedsub_chan;
        TH2I *st_adc_peak_chan;
        TH2I *st_adc_peak_pedsub_chan;
        TH2I *st_adc_time_chan;
        TH2I *st_adc_pedestal_chan;
        TH2I *st_adc_quality_chan;

        TH1I *st_adc_emudelta_integral;
        TH1I *st_adc_emudelta_peak;
        TH1I *st_adc_emudelta_pedestal;
        TH1I *st_adc_emudelta_coarsetime;
        TH1I *st_adc_emudelta_finetime;


		//------------------------ TAGH -----------------------
        TH1I *tag_num_events;

        TH1I *tagh_adc_integral;
        TH1I *tagh_adc_integral_pedsub;
        TH1I *tagh_adc_peak;
        TH1I *tagh_adc_peak_pedsub;
        TH1I *tagh_adc_time;
        TH1I *tagh_adc_pedestal;
        TH1I *tagh_adc_quality;
        TH2I *tagh_adc_integral_chan;
        TH2I *tagh_adc_integral_pedsub_chan;
        TH2I *tagh_adc_peak_chan;
        TH2I *tagh_adc_peak_pedsub_chan;
        TH2I *tagh_adc_time_chan;
        TH2I *tagh_adc_pedestal_chan;
        TH2I *tagh_adc_quality_chan;

        TH1I *tagh_adc_emudelta_integral;
        TH1I *tagh_adc_emudelta_peak;
        TH1I *tagh_adc_emudelta_pedestal;
        TH1I *tagh_adc_emudelta_coarsetime;
        TH1I *tagh_adc_emudelta_finetime;


		//------------------------ TAGM -----------------------
        TH1I *tagm_adc_integral;
        TH1I *tagm_adc_integral_pedsub;
        TH1I *tagm_adc_peak;
        TH1I *tagm_adc_peak_pedsub;
        TH1I *tagm_adc_time;
        TH1I *tagm_adc_pedestal;
        TH1I *tagm_adc_quality;
        TH2I *tagm_adc_integral_chan;
        TH2I *tagm_adc_integral_pedsub_chan;
        TH2I *tagm_adc_peak_chan;
        TH2I *tagm_adc_peak_pedsub_chan;
        TH2I *tagm_adc_time_chan;
        TH2I *tagm_adc_pedestal_chan;
        TH2I *tagm_adc_quality_chan;

        TH1I *tagm_adc_emudelta_integral;
        TH1I *tagm_adc_emudelta_peak;
        TH1I *tagm_adc_emudelta_pedestal;
        TH1I *tagm_adc_emudelta_coarsetime;
        TH1I *tagm_adc_emudelta_finetime;


		//------------------------ TOF ------------------------
        TH1I *tof_num_events;
        TH1I *tof_adc_integral;
        TH1I *tof_adc_integral_pedsub;
        TH1I *tof_adc_peak;
        TH1I *tof_adc_peak_pedsub;
        TH1I *tof_adc_time;
        TH1I *tof_adc_pedestal;
        TH1I *tof_adc_quality;
        TH2I *tof_adc_integral_chan;
        TH2I *tof_adc_integral_pedsub_chan;
        TH2I *tof_adc_peak_chan;
        TH2I *tof_adc_peak_pedsub_chan;
        TH2I *tof_adc_time_chan;
        TH2I *tof_adc_pedestal_chan;
        TH2I *tof_adc_quality_chan;

        TH1I *tof_adc_emudelta_integral;
        TH1I *tof_adc_emudelta_peak;
        TH1I *tof_adc_emudelta_pedestal;
        TH1I *tof_adc_emudelta_coarsetime;
        TH1I *tof_adc_emudelta_finetime;



	private:
		jerror_t init(void);
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);
		jerror_t erun(void);
		jerror_t fini(void);
};

#endif // _JEventProcessor_lowlevel_online_

