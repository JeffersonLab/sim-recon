#include <vector>
#include "TH1F.h"
#include "TH2F.h"

extern vector<vector<float> >fdc_smear_parms; 

extern vector<unsigned int> NCDC_STRAWS;
extern vector<double> CDC_RING_RADIUS;

extern vector<double> FDC_LAYER_Z;
extern double FDC_RATE_COEFFICIENT;

// Do we or do we not add noise hits
extern bool ADD_NOISE;

// Do we or do we not smear real hits
extern bool SMEAR_HITS;

// Use or ignore random number seeds found in HDDM file
extern bool IGNORE_SEEDS;

// The error on the drift time in the CDC. The drift times
// for the actual CDC hits coming from the input file
// are smeared by a gaussian with this sigma.
extern double CDC_TDRIFT_SIGMA;

// The time window for which CDC hits are accumulated.
// This is used to determine the number of background
// hits in the CDC for a given event.
extern double CDC_TIME_WINDOW;

// The error in the energy deposition measurement in the CDC due to pedestal noise
extern double CDC_PEDESTAL_SIGMA;

// Number of sigmas above threshold for sparcification of CDC data
extern double CDC_THRESHOLD_FACTOR;
 
// If the following flag is true, then include the drift-distance
// dependency on the error in the FDC position. Otherwise, use a
// flat distribution given by the FDC_TDRIFT_SIGMA below.
extern bool FDC_USE_PARAMETERIZED_SIGMA;

// The error on the drift time in the FDC. The drift times
// for the actual FDC hits coming from the input file
// are smeared by a gaussian with this sigma.
extern double FDC_TDRIFT_SIGMA;

// The error in the distance along the wire as measured by
// the cathodes. This should NOT include the Lorentz
// effect which is already included in hdgeant. It
// should include any fluctuations due to ion trail density
// etc.
extern double FDC_CATHODE_SIGMA;

// The FDC pedestal noise is used to smear the cathode ADC
// values such that the position along the wire has the resolution
// specified by FDC_CATHODE_SIGMA.
extern double FDC_PED_NOISE; //pC (calculated from FDC_CATHODE_SIGMA in SmearFDC)

// Number of sigmas above threshold for sparcification of FDC data
extern double FDC_THRESHOLD_FACTOR;
 
// If energy loss was turned off in the FDC then the pedestal
// noise will not be scaled properly to give the nominal 200 micron
// resolution along the wires. This flag is used to indicated
// the magic scale factor should be applied to FDC_PED_NOISE
// when it is calculated below to get the correct resolution.
extern bool FDC_ELOSS_OFF;

// Time window for acceptance of FDC hits
extern double FDC_TIME_WINDOW;

// Fraction of FDC hits to randomly drop (0=drop nothing 1=drop everything)
extern double FDC_HIT_DROP_FRACTION;

// FDC discriminator threshold in keV
extern double FDC_THRESH_KEV;

// Single FTOF bar energy threshold (applied after smearing)
extern double FTOF_BAR_THRESHOLD;

// Single STC paddle energy threshold (applied after smearing)
extern double STC_PADDLE_THRESHOLD;

// Single pair spectrometer energy threshold (applied after smearing)
extern double PSC_THRESHOLD;

// PS counter resolutions and energy reponses
extern double PS_SIGMA; // ns
extern double PSC_SIGMA; //ns
extern double PS_NPIX_PER_GEV;
extern double PSC_PHOTONS_PERMEV;

// Tagging counter time resolutions (s)
extern double TAGM_TSIGMA;
extern double TAGH_TSIGMA;

// Tagging counter fADC responses
extern double TAGM_FADC_TSIGMA;
extern double TAGH_FADC_TSIGMA;
extern double TAGM_NPIX_PER_GEV;
extern double TAGH_NPE_PER_GEV;

// Photon-statistics factor for smearing hit energy (from Criss's MC)
// (copied from DFCALMCResponse_factory.cc 7/2/2009 DL)
extern double FCAL_PHOT_STAT_COEF;

// Single block energy threshold (applied after smearing)
extern double FCAL_BLOCK_THRESHOLD;

// Forward TOF resolution
extern double TOF_SIGMA;
extern double TOF_PHOTONS_PERMEV;

// Start counter resolution
extern double START_SIGMA ;
extern double START_PHOTONS_PERMEV;

// Pair spectrometer resolution
extern double PS_SIGMA;
extern double PS_PHOTONS_PERMEV;


extern bool SMEAR_BCAL;

// Beginning of readout window
extern double TRIGGER_LOOKBACK_TIME;

// Mutex used to control accessing the ROOT global memory
extern pthread_mutex_t root_mutex;

// Flag used specifically for BCAL
extern bool SMEAR_BCAL;

// The following are all false by default, but can be
// set to true via command line parameters. Setting
// one of these to true will turn OFF the feature.
extern bool NO_E_SMEAR;
extern bool NO_T_SMEAR;
extern bool NO_DARK_PULSES;
extern bool NO_SAMPLING_FLUCTUATIONS;
extern bool NO_SAMPLING_FLOOR_TERM;
extern bool NO_POISSON_STATISTICS;
extern bool NO_THRESHOLD_CUT;

extern double BCAL_FADC_TIME_RESOLUTION; // ns
extern double BCAL_TDC_TIME_RESOLUTION; // ns
extern double BCAL_ADC_THRESHOLD_MEV; // MeV

// setup response parameters
extern double BCAL_SAMPLINGCOEFA;               // 0.042 (from calibDB BCAL/bcal_parms)
extern double BCAL_SAMPLINGCOEFB;               // 0.013 (from calibDB BCAL/bcal_parms)
extern double BCAL_TWO_HIT_RESOL;               // 50. (from calibDB BCAL/bcal_parms)
extern double BCAL_mevPerPE;                    // (defined below)
extern double BCAL_NS_PER_ADC_COUNT;            // 0.0625 (defined in mcsmear.cc)
extern double BCAL_NS_PER_TDC_COUNT;            // 0.0559 (defined in mcsmear.cc)
extern double BCAL_MEV_PER_ADC_COUNT;           // 0.029 (defined in mcsmear.cc)

extern vector<vector<double> > attenuation_parameters; // Avg. of 525 (from calibDB BCAL/attenuation_parameters)
extern vector<double> effective_velocities;     // 16.75 (from calibDB BCAL/effective_velocities)

extern int BCAL_NUM_MODULES;
extern int BCAL_NUM_LAYERS;
extern int BCAL_NUM_SECTORS;

extern double BCAL_BASE_TIME_OFFSET;            // -100.0 (from calibDB BCAL/base_time_offset)
extern double BCAL_TDC_BASE_TIME_OFFSET;        // -100.0 (from calibDB BCAL/base_time_offset)

extern double BCAL_XTALK_FRACT;

// The following are not currently in use
extern double BCAL_DARKRATE_GHZ;                // 0.0176 (from calibDB BCAL/bcal_parms) for 4x4 array
extern double BCAL_DEVICEPDE;                   // 0.21   (from calibDB BCAL/bcal_parms)
extern double BCAL_SAMPLING_FRACT;              // 0.095  (from calibDB BCAL/bcal_parms)
extern double BCAL_PHOTONSPERSIDEPERMEV_INFIBER;// 75 (from calibDB BCAL/bcal_parms)
extern double BCAL_SIGMA_SIG_RELATIVE;          // 0.105  (from calibDB BCAL/bcal_parms)
extern double BCAL_SIGMA_PED_RELATIVE;          // 0.139  (from calibDB BCAL/bcal_parms)
extern double BCAL_SIPM_GAIN_VARIATION;         // 0.04   (from calibDB BCAL/bcal_parms)
extern double BCAL_INTWINDOW_NS;                // 100    (from calibDB BCAL/bcal_parms)
extern double BCAL_AVG_DARK_DIGI_VALS_PER_EVENT;// 240 used to set thresholds
extern double BCAL_TIMEDIFFCOEFA;               // 0.07 * sqrt( 2 ) (from calibDB BCAL/bcal_parms)
extern double BCAL_TIMEDIFFCOEFB;               // 0.00 * sqrt( 2 ) (from calibDB BCAL/bcal_parms)
// end "not-in-use"
