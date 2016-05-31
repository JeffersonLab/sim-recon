// Classes to store configuration information for mcsmear

#ifndef _MCSMEAR_CONFIG_H_
#define _MCSMEAR_CONFIG_H_

#include "units.h"

#include <DANA/DApplication.h>
#include <JANA/JEventLoop.h>

using namespace jana;

// Overall configuration parameters
class mcsmear_config_t 
{
  public:
	mcsmear_config_t() {
		// default values
		DROP_TRUTH_HITS = false;
		ADD_NOISE      = false;
		SMEAR_HITS     = true;
		//SMEAR_BCAL     = true;
		IGNORE_SEEDS   = false;
		
		TRIGGER_LOOKBACK_TIME = -100; // ns
	}
	~mcsmear_config_t() {}

	bool ADD_NOISE;
	bool DROP_TRUTH_HITS;
	bool SMEAR_HITS;
	//bool SMEAR_BCAL;
	//bool FDC_ELOSS_OFF;
	bool IGNORE_SEEDS;
	double TRIGGER_LOOKBACK_TIME;
};


// Detector-specific classes
class bcal_config_t 
{
  public:
	bcal_config_t(JEventLoop *loop) {
		// default values
 		BCAL_DARKRATE_GHZ         = 0.;// 0.0176 (from calibDB BCAL/bcal_parms) for 4x4 array
 		BCAL_SIGMA_SIG_RELATIVE   = 0.;// 0.105  (from calibDB BCAL/bcal_parms)
 		BCAL_SIGMA_PED_RELATIVE   = 0.;// 0.139  (from calibDB BCAL/bcal_parms)
 		BCAL_SIPM_GAIN_VARIATION  = 0.;// 0.04   (from calibDB BCAL/bcal_parms)
 		BCAL_XTALK_FRACT          = 0.;// 0.157  (from calibDB BCAL/bcal_parms)
 		BCAL_INTWINDOW_NS         = 0.;// 100    (from calibDB BCAL/bcal_parms)
 		BCAL_DEVICEPDE            = 0.;// 0.21   (from calibDB BCAL/bcal_parms)
 		BCAL_SAMPLING_FRACT       = 0.;// 0.095  (from calibDB BCAL/bcal_parms)
 		BCAL_PHOTONSPERSIDEPERMEV_INFIBER = 0.0;// 75 (from calibDB BCAL/bcal_parms)
 		BCAL_AVG_DARK_DIGI_VALS_PER_EVENT = 0.0; // 240 used to set thresholds
 		BCAL_SAMPLINGCOEFA        = 0.0; // 0.042 (from calibDB BCAL/bcal_parms)
 		BCAL_SAMPLINGCOEFB        = 0.0; // 0.013 (from calibDB BCAL/bcal_parms)
 		BCAL_TIMEDIFFCOEFA        = 0.0; // 0.07 * sqrt( 2 ) (from calibDB BCAL/bcal_parms)
 		BCAL_TIMEDIFFCOEFB        = 0.0; // 0.00 * sqrt( 2 ) (from calibDB BCAL/bcal_parms)
 		BCAL_TWO_HIT_RESOL        = 0.0; // 50. (from calibDB BCAL/bcal_parms)
 		BCAL_mevPerPE             = 0.31; // Energy corresponding to one pixel firing in MeV - FIX - in CCDB??

		// FIX - Pull from geometry?
 		BCAL_NUM_MODULES = 48;
 		BCAL_NUM_LAYERS = 4;
 		BCAL_NUM_SECTORS = 4;

 		BCAL_BASE_TIME_OFFSET     = 0; // -100.0 (from calibDB BCAL/base_time_offset)
 		BCAL_TDC_BASE_TIME_OFFSET = 0; // -100.0 (from calibDB BCAL/base_time_offset)

		// FIX - put in CCDB?
 		BCAL_ADC_THRESHOLD_MEV    = 2.2;  // MeV (To be updated/improved)
 		BCAL_FADC_TIME_RESOLUTION = 0.3;  // ns (To be updated/improved)
 		BCAL_TDC_TIME_RESOLUTION  = 0.3;  // ns (To be updated/improved)
 		BCAL_MEV_PER_ADC_COUNT    = 0.029;  // MeV per integrated ADC count (based on Spring 2015 calibrations)
 		BCAL_NS_PER_ADC_COUNT     = 0.0;  // 0.0625 ns per ADC count (from calibDB BCAL/digi_scales)
 		BCAL_NS_PER_TDC_COUNT     = 0.0;  // 0.0559 ns per TDC count (from calibDB BCAL/digi_scales)

		// BCAL flags
 		NO_T_SMEAR = false;
 		NO_DARK_PULSES = false;
 		NO_SAMPLING_FLUCTUATIONS = false;
 		NO_SAMPLING_FLOOR_TERM = false;
 		NO_POISSON_STATISTICS = false;

		// Load parameters from CCDB
     	cout << "get BCAL/bcal_parms parameters from CCDB..." << endl;
     	map<string, double> bcalparms;
     	if(loop->GetCalib("BCAL/bcal_parms", bcalparms)) {
     		jerr << "Problem loading BCAL/bcal_parms from CCDB!" << endl;
     	} else {
     		BCAL_DARKRATE_GHZ         = bcalparms["BCAL_DARKRATE_GHZ"];
     		BCAL_SIGMA_SIG_RELATIVE   = bcalparms["BCAL_SIGMA_SIG_RELATIVE"];
     		BCAL_SIGMA_PED_RELATIVE   = bcalparms["BCAL_SIGMA_PED_RELATIVE"];
     		BCAL_SIPM_GAIN_VARIATION  = bcalparms["BCAL_SIPM_GAIN_VARIATION"];
     		BCAL_XTALK_FRACT          = bcalparms["BCAL_XTALK_FRACT"];
     		BCAL_INTWINDOW_NS         = bcalparms["BCAL_INTWINDOW_NS"];
     		BCAL_DEVICEPDE         	  = bcalparms["BCAL_DEVICEPDE"];
     		BCAL_SAMPLING_FRACT       = bcalparms["BCAL_SAMPLING_FRACT"];
     		BCAL_AVG_DARK_DIGI_VALS_PER_EVENT   = bcalparms["BCAL_AVG_DARK_DIGI_VALS_PER_EVENT"];
     		BCAL_PHOTONSPERSIDEPERMEV_INFIBER   = bcalparms["BCAL_PHOTONSPERSIDEPERMEV_INFIBER"];
     		BCAL_SAMPLINGCOEFA 		  = bcalparms["BCAL_SAMPLINGCOEFA"];
     		BCAL_SAMPLINGCOEFB 		  = bcalparms["BCAL_SAMPLINGCOEFB"];
     		BCAL_TIMEDIFFCOEFA 		  = bcalparms["BCAL_TIMEDIFFCOEFA"];
     		BCAL_TIMEDIFFCOEFB 		  = bcalparms["BCAL_TIMEDIFFCOEFB"];
     		BCAL_TWO_HIT_RESOL 		  = bcalparms["BCAL_TWO_HIT_RESOL"];
		}
		
    	cout << "Get BCAL/attenuation_parameters from CCDB..." <<endl;
     	vector< vector<double> > in_atten_parameters;
     	if(loop->GetCalib("BCAL/attenuation_parameters", in_atten_parameters)) {
     		jerr << "Problem loading BCAL/bcal_parms from CCDB!" << endl;
		} else {
     		attenuation_parameters.clear();

     		for (unsigned int i = 0; i < in_atten_parameters.size(); i++) {
     			attenuation_parameters.push_back( in_atten_parameters.at(i) );
     		}
		}
		
     	cout << "Get BCAL/effective_velocities parameters from CCDB..." << endl;
     	vector <double> in_effective_velocities;
     	if(loop->GetCalib("BCAL/effective_velocities", in_effective_velocities)) {
     		jerr << "Problem loading BCAL/effective_velocities from CCDB!" << endl;
     	} else {
     		for (unsigned int i = 0; i < in_effective_velocities.size(); i++) {
       			effective_velocities.push_back( in_effective_velocities.at(i) );
     		}
     	}
   
     	cout << "Get BCAL/digi_scales parameters from CCDB..." << endl;
     	map<string, double> bcaldigiscales;
     	if(loop->GetCalib("BCAL/digi_scales", bcaldigiscales)) {
     		jerr << "Problem loading BCAL/digi_scales from CCDB!" << endl;
     	} else {
     		BCAL_NS_PER_ADC_COUNT = bcaldigiscales["BCAL_ADC_TSCALE"];
     		BCAL_NS_PER_TDC_COUNT = bcaldigiscales["BCAL_TDC_SCALE"];
   		}

     	cout << "Get BCAL/base_time_offset parameters from CCDB..." << endl;
     	map<string, double> bcaltimeoffsets;
     	if(loop->GetCalib("BCAL/base_time_offset", bcaltimeoffsets)) {
     		jerr << "Problem loading BCAL/base_time_offset from CCDB!" << endl;
 		} else {
     		BCAL_BASE_TIME_OFFSET = bcaltimeoffsets["BCAL_BASE_TIME_OFFSET"];
    		BCAL_TDC_BASE_TIME_OFFSET = bcaltimeoffsets["BCAL_TDC_BASE_TIME_OFFSET"];
   		}

	}
	~bcal_config_t() {}

	double BCAL_DARKRATE_GHZ;
	double BCAL_SIGMA_SIG_RELATIVE;
	double BCAL_SIGMA_PED_RELATIVE;
	double BCAL_SIPM_GAIN_VARIATION;
	double BCAL_XTALK_FRACT;
	double BCAL_INTWINDOW_NS;
	double BCAL_DEVICEPDE;
	double BCAL_SAMPLING_FRACT;
	double BCAL_PHOTONSPERSIDEPERMEV_INFIBER;
	double BCAL_AVG_DARK_DIGI_VALS_PER_EVENT;
	double BCAL_SAMPLINGCOEFA;
	double BCAL_SAMPLINGCOEFB;
	double BCAL_TIMEDIFFCOEFA;
	double BCAL_TIMEDIFFCOEFB;
	double BCAL_TWO_HIT_RESOL;
	double BCAL_mevPerPE;

	int BCAL_NUM_MODULES;
	int BCAL_NUM_LAYERS;
	int BCAL_NUM_SECTORS;

	double BCAL_BASE_TIME_OFFSET;
	double BCAL_TDC_BASE_TIME_OFFSET;

	double BCAL_ADC_THRESHOLD_MEV;
	double BCAL_FADC_TIME_RESOLUTION;
	double BCAL_TDC_TIME_RESOLUTION;
	double BCAL_MEV_PER_ADC_COUNT;
	double BCAL_NS_PER_ADC_COUNT;
	double BCAL_NS_PER_TDC_COUNT;

	// BCAL flags
	bool NO_T_SMEAR;
	bool NO_DARK_PULSES;
	bool NO_SAMPLING_FLUCTUATIONS;
	bool NO_SAMPLING_FLOOR_TERM;
	bool NO_POISSON_STATISTICS;

	vector<vector<double> > attenuation_parameters; // Avg. of 525 (from calibDB BCAL/attenuation_parameters)
	vector<double> effective_velocities; // 16.75 (from calibDB BCAL/effective_velocities)

};

class fcal_config_t 
{
  public:
	fcal_config_t(JEventLoop *loop) {
		// default values
		FCAL_PHOT_STAT_COEF   = 0.0; //0.035;
		FCAL_BLOCK_THRESHOLD  = 0.0; //20.0*k_MeV;
		// FCAL_TSIGMA           = 0.0; // 200 ps
		FCAL_TSIGMA           = 0.2; // 200 ps - FIX THIS!!

		// Get values from CCDB
		cout << "Get FCAL/fcal_parms parameters from CCDB..." << endl;
     	map<string, double> fcalparms;
     	if(loop->GetCalib("FCAL/fcal_parms", fcalparms)) { 
     		jerr << "Problem loading FCAL/fcal_parms from CCDB!" << endl;
     	} else {
       		FCAL_PHOT_STAT_COEF   = fcalparms["FCAL_PHOT_STAT_COEF"]; 
       		FCAL_BLOCK_THRESHOLD  = fcalparms["FCAL_BLOCK_THRESHOLD"];
		}
		
	}
	~fcal_config_t();

	double FCAL_PHOT_STAT_COEF;
	double FCAL_BLOCK_THRESHOLD;
	double FCAL_TSIGMA;
};

class cdc_config_t 
{
  public:
	cdc_config_t(JEventLoop *loop) {
		// default values
		CDC_TDRIFT_SIGMA      = 0.0;
 		CDC_TIME_WINDOW       = 0.0;
 		CDC_PEDESTAL_SIGMA    = 0.0;
 		CDC_THRESHOLD_FACTOR  = 0.0;
 		
 		// load data from CCDB
 		jerr << "get CDC/cdc_parms parameters from CCDB..." << endl;
     	map<string, double> cdcparms;
     	if(loop->GetCalib("CDC/cdc_parms", cdcparms)) {
     		jerr << "Problem loading CDC/cdc_parms from CCDB!" << endl;
     	} else {
       		CDC_TDRIFT_SIGMA   = cdcparms["CDC_TDRIFT_SIGMA"]; 
       		CDC_TIME_WINDOW    = cdcparms["CDC_TIME_WINDOW"];
       		CDC_PEDESTAL_SIGMA = cdcparms["CDC_PEDESTAL_SIGMA"]; 
       		CDC_THRESHOLD_FACTOR = cdcparms["CDC_THRESHOLD_FACTOR"];
		}

	}
	~cdc_config_t();

	double CDC_TDRIFT_SIGMA;
	double CDC_TIME_WINDOW;
	double CDC_PEDESTAL_SIGMA;
	double CDC_THRESHOLD_FACTOR; // number of pedestal sigmas for determining sparcification threshold

};

class fdc_config_t 
{
  public:
	fdc_config_t(JEventLoop *loop) {
		// default values
		FDC_TDRIFT_SIGMA      = 0.0;
 		FDC_CATHODE_SIGMA     = 0.0;
 		FDC_PED_NOISE         = 0.0;
 		FDC_THRESHOLD_FACTOR  = 0.0;
 		FDC_TIME_WINDOW       = 0.0;
 		FDC_THRESH_KEV        = 0.0;

   		// Calculate ped noise level based on position resolution
   		//   FDC_PED_NOISE = -0.004594 + 0.008711*FDC_CATHODE_SIGMA +
   		//                    0.000010*FDC_CATHODE_SIGMA*FDC_CATHODE_SIGMA; //pC
   		FDC_PED_NOISE = -0.0938 + 0.0485*FDC_CATHODE_SIGMA;

		// load data from CCDB
		cout << "Get FDC/fdc_parms parameters from CCDB..." << endl;
     	map<string, double> fdcparms;
     	if(loop->GetCalib("FDC/fdc_parms", fdcparms)) {
     		jerr << "Problem loading FDC/fdc_parms from CCDB!" << endl;
     	} else {
       		FDC_TDRIFT_SIGMA      = fdcparms["FDC_TDRIFT_SIGMA"];
       		FDC_CATHODE_SIGMA     = fdcparms["FDC_CATHODE_SIGMA"];
       		FDC_THRESHOLD_FACTOR  = fdcparms["FDC_THRESHOLD_FACTOR"];
     		//FDC_PED_NOISE         = fdcparms["FDC_PED_NOISE"];  // ???
       		FDC_TIME_WINDOW       = fdcparms["FDC_TIME_WINDOW"];
       		//FDC_HIT_DROP_FRACTION = fdcparms["FDC_HIT_DROP_FRACTION"];   // ???
       		FDC_THRESH_KEV 		  = fdcparms["FDC_THRESH_KEV"]; 
		}
	}
	~fdc_config_t();

	double FDC_TDRIFT_SIGMA;
	double FDC_CATHODE_SIGMA;
	double FDC_PED_NOISE;         // in pC calculated in above
	double FDC_THRESHOLD_FACTOR;  // number of pedestal sigmas for determining sparcification threshold
	//double FDC_HIT_DROP_FRACTION = 0.0; // 1000.0E-9
	double FDC_TIME_WINDOW;
	double FDC_THRESH_KEV;        // fdc anode discriminator threshold

};

class ccal_config_t 
{
  public:
	ccal_config_t(JEventLoop *loop) {
		// default values
		// (This is just a rough estimate 11/30/2010 DL)
		CCAL_PHOT_STAT_COEF = 0.035/2.0;
		CCAL_BLOCK_THRESHOLD = 20.0*k_MeV;
		CCAL_SIGMA = 200.0e-3;
	}
	~ccal_config_t();

	// Time smearing factor
	double CCAL_SIGMA;

	// Photon-statistics factor for smearing hit energy for CompCal
	double CCAL_PHOT_STAT_COEF;

	// Single block energy threshold (applied after smearing)
	double CCAL_BLOCK_THRESHOLD;

};

class ftof_config_t 
{
  public:
	ftof_config_t(JEventLoop *loop) {
		// default values
 		TOF_SIGMA = 100.*k_psec;
 		TOF_PHOTONS_PERMEV = 400.;
 		TOF_BAR_THRESHOLD    = 0.0;

		// Load data from CCDB
     	cout<<"Get TOF/tof_parms parameters from CCDB..."<<endl;
     	map<string, double> tofparms;
     	if(loop->GetCalib("TOF/tof_parms", tofparms)) {
     		jerr << "Problem loading TOF/tof_parms from CCDB!" << endl;
     		return;
     	}
     	
     	TOF_SIGMA =  tofparms["TOF_SIGMA"];
     	TOF_PHOTONS_PERMEV =  tofparms["TOF_PHOTONS_PERMEV"];
	
	}
	~ftof_config_t();

	double TOF_SIGMA;
	double TOF_PHOTONS_PERMEV;
	double TOF_BAR_THRESHOLD;

};

class sc_config_t 
{
  public:
	sc_config_t(JEventLoop *loop) {
		// default values
		START_SIGMA           = 0.0; // 300ps
		START_PHOTONS_PERMEV  = 0.0; // used to be 8000 should be more like 200
		START_PADDLE_THRESHOLD  = 0.0;

		// Load data from CCDB
     	cout << "Get START_COUNTER/start_parms parameters from CCDB..." << endl;
     	map<string, double> startparms;
     	if(loop->GetCalib("START_COUNTER/start_parms", startparms)) {
			jerr << "Problem loading START_COUNTER/start_parms from CCDB!" << endl;
		} else {
     		START_SIGMA = startparms["START_SIGMA"] ;
     		START_PHOTONS_PERMEV = startparms["START_PHOTONS_PERMEV"];
		}
		
	}
	~sc_config_t();

	double START_SIGMA;
	double START_PHOTONS_PERMEV;
	double START_PADDLE_THRESHOLD;

};

class tagm_config_t 
{
  public:
	tagm_config_t(JEventLoop *loop) {
		// default values
		TAGM_TSIGMA = 0.200;        // ns
		TAGM_FADC_TSIGMA = 0.350;   // ns
		TAGM_NPIX_PER_GEV = 1.e5;
	}
	~tagm_config_t();

	double TAGM_TSIGMA;
	double TAGM_FADC_TSIGMA;
	double TAGM_NPIX_PER_GEV;

};

class tagh_config_t 
{
  public:
	tagh_config_t(JEventLoop *loop) {
		// default values
		TAGH_TSIGMA = 0.350;        // ns
		TAGH_FADC_TSIGMA = 0.450;   // ns
		TAGH_NPE_PER_GEV = 5.e5;

	}
	~tagh_config_t();

	double TAGH_TSIGMA;
	double TAGH_FADC_TSIGMA;
	double TAGH_NPE_PER_GEV;

};

class ps_config_t 
{
  public:
	ps_config_t(JEventLoop *loop) {
		// default values
		PS_SIGMA = 0.200; // ns
		PS_NPIX_PER_GEV = 1.e5;
		PS_THRESHOLD          = 0.0;
	}
	~ps_config_t();

	double PS_SIGMA;
	double PS_NPIX_PER_GEV;
	double PS_THRESHOLD;

};

class psc_config_t 
{
  public:
	psc_config_t(JEventLoop *loop){
		// default values
		PSC_SIGMA = 0.200; //ns
		PSC_PHOTONS_PERMEV = 5.e5;
		PSC_THRESHOLD         = 0.0;
	}
	~psc_config_t();

	double PSC_SIGMA;
	double PSC_PHOTONS_PERMEV;
	double PSC_THRESHOLD;

};

class fmwpc_config_t 
{
  public:
	fmwpc_config_t(JEventLoop *loop) {
		// default values
		FMWPC_TSIGMA = 10.0;  // ns
 		FMWPC_ASIGMA = 0.5E-6;
 		FMWPC_THRESHOLD = 0.0;
	}
	~fmwpc_config_t();

	// FMWPC resolutions and threshold
	double FMWPC_TSIGMA;
	double FMWPC_ASIGMA;
	double FMWPC_THRESHOLD;

};



#endif  // _MCSMEAR_CONFIG_H_
