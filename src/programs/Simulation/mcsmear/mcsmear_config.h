// Classes to store configuration information for mcsmear

#ifndef _MCSMEAR_CONFIG_H_
#define _MCSMEAR_CONFIG_H_

// Overall configuration parameters
class mcsmear_config_t 
{
  public:
	mcsmear_config_t() {
		// default values
		DROP_TRUTH_HITS = false;
		//ADD_NOISE      = false;
		SMEAR_HITS     = true;
		SMEAR_BCAL     = true;
		IGNORE_SEEDS   = false;
		
		TRIGGER_LOOKBACK_TIME = -100; // ns
	}
	~mcsmear_config_t() {}

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
	bcal_config_t() {
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

		vector<vector<double> > attenuation_parameters; // Avg. of 525 (from calibDB BCAL/attenuation_parameters)
		vector<double> effective_velocities; // 16.75 (from calibDB BCAL/effective_velocities)

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
	fcal_config_t() {
		// default values
		FCAL_PHOT_STAT_COEF   = 0.0; //0.035;
		FCAL_BLOCK_THRESHOLD  = 0.0; //20.0*k_MeV;
		// FCAL_TSIGMA           = 0.0; // 200 ps
		FCAL_TSIGMA           = 0.2; // 200 ps - FIX THIS!!
	}
	~fcal_config_t();

	double FCAL_PHOT_STAT_COEF;
	double FCAL_BLOCK_THRESHOLD;
	double FCAL_TSIGMA;
};

class cdc_config_t 
{
  public:
	cdc_config_t() {
		// default values
		CDC_TDRIFT_SIGMA      = 0.0;
 		CDC_TIME_WINDOW       = 0.0;
 		CDC_PEDESTAL_SIGMA    = 0.0;
 		CDC_THRESHOLD_FACTOR  = 0.0;
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
	fdc_config_t() {
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
	ccal_config_t() {
		// default values
		// (This is just a rough estimate 11/30/2010 DL)
		CCAL_PHOT_STAT_COEF = 0.035/2.0;
		CCAL_BLOCK_THRESHOLD = 20.0*k_MeV;
		CCAL_SIGMA = 200.0e-3;
	}
	~ccal_config_t();

	// Photon-statistics factor for smearing hit energy for CompCal
	double CCAL_PHOT_STAT_COEF;

	// Single block energy threshold (applied after smearing)
	double CCAL_BLOCK_THRESHOLD;

};

class ftof_config_t 
{
  public:
	ftof_config_t() {
		// default values
 		TOF_SIGMA = 100.*k_psec;
 		TOF_PHOTONS_PERMEV = 400.;
 		TOF_BAR_THRESHOLD    = 0.0;

		// Load from CCDB
     	cout<<"get TOF/tof_parms parameters from CCDB..."<<endl;
     	map<string, double> tofparms;
     	jcalib->Get("TOF/tof_parms", tofparms);
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
	sc_config_t() {
		// default values
		START_SIGMA           = 0.0; // 300ps
		START_PHOTONS_PERMEV  = 0.0; // used to be 8000 should be more like 200
		START_PADDLE_THRESHOLD  = 0.0;

	}
	~sc_config_t();

	double START_SIGMA;
	double START_PHOTONS_PERMEV;
	double START_PADDLE_THRESHOLD;

};

class tagm_config_t 
{
  public:
	tagm_config_t() {
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
	tagh_config_t() {
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
	ps_config_t() {
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
	psc_config_t(){
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
	fmwpc_config_t() {
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
