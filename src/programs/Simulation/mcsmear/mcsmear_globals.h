#include <vector>
#include "TH1F"
#include "TH2F"

extern vector<vector<float> >fdc_smear_parms; 
extern TH2F *fdc_drift_time_smear_hist;
extern TH2F *fdc_drift_dist_smear_hist;
extern TH2F *fdc_drift_time;
extern TH1F *fdc_cathode_charge;
extern TH2F *cdc_drift_time;
extern TH1F *cdc_charge;
extern TH1F *fdc_anode_mult;
extern TH2F *cdc_drift_smear;

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

// FMWPC resolutions and threshold
extern double FMWPC_TSIGMA;
extern double FMWPC_ASIGMA;
extern double FMWPC_THRESHOLD;

extern bool SMEAR_BCAL;

// Beginning of readout window
extern double TRIGGER_LOOKBACK_TIME;

extern bool DROP_TRUTH_HITS;


