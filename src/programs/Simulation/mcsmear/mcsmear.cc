// $Id: mcsmear.cc 19023 2015-07-14 20:23:27Z beattite $
//
// Created June 22, 2005  David Lawrence


// The following flag can be used to switch from the classic mode where
// the event loop is implemented in main() to a JANA based event-loop.
// NOTE:  for consistency, one should set the value of #define JANA_ENABLED
// in smear.cc to the same value!
#define USE_JANA 1

#include <iostream>
#include <iomanip>
using namespace std;

#include <TF1.h>
#include <TFile.h>
#include <TH2.h>
#include <TH1.h>

#include <signal.h>
#include <time.h>

#if USE_JANA
#include <DANA/DApplication.h>
#include "MyProcessor.h"
#include "JFactoryGenerator_ThreadCancelHandler.h"

#include <CDC/DCDCWire.h>
//#include <FDC/DFDCWire.h>
#include <HDGEOMETRY/DGeometry.h>

#endif

#include "units.h"
#include "HDDM/hddm_s.hpp"

void Smear(hddm_s::HDDM *record);
void ParseCommandLineArguments(int narg, char* argv[]);
void Usage(void);

extern void SetSeeds(const char *vals);

#if ! USE_JANA
void ctrlCHandleMCSmear(int x);
#endif

char *INFILENAME = NULL;
char *OUTFILENAME = NULL;
int QUIT = 0;

bool ADD_NOISE      = false;
bool SMEAR_HITS     = true;
bool SMEAR_BCAL     = true;
bool FDC_ELOSS_OFF  = false;
bool IGNORE_SEEDS   = false;

// setup response parameters
double BCAL_DARKRATE_GHZ         = 0.;// 0.0176 (from calibDB BCAL/bcal_parms) for 4x4 array
double BCAL_SIGMA_SIG_RELATIVE   = 0.;// 0.105  (from calibDB BCAL/bcal_parms)
double BCAL_SIGMA_PED_RELATIVE   = 0.;// 0.139  (from calibDB BCAL/bcal_parms)
double BCAL_SIPM_GAIN_VARIATION  = 0.;// 0.04   (from calibDB BCAL/bcal_parms)
double BCAL_XTALK_FRACT          = 0.;// 0.157  (from calibDB BCAL/bcal_parms)
double BCAL_INTWINDOW_NS         = 0.;// 100    (from calibDB BCAL/bcal_parms)
double BCAL_DEVICEPDE            = 0.;// 0.21   (from calibDB BCAL/bcal_parms)
double BCAL_SAMPLING_FRACT       = 0.;// 0.095  (from calibDB BCAL/bcal_parms)
double BCAL_PHOTONSPERSIDEPERMEV_INFIBER = 0.0;// 75 (from calibDB BCAL/bcal_parms)
double BCAL_AVG_DARK_DIGI_VALS_PER_EVENT = 0.0; // 240 used to set thresholds
double BCAL_SAMPLINGCOEFA        = 0.0; // 0.042 (from calibDB BCAL/bcal_parms)
double BCAL_SAMPLINGCOEFB        = 0.0; // 0.013 (from calibDB BCAL/bcal_parms)
double BCAL_TIMEDIFFCOEFA        = 0.0; // 0.07 * sqrt( 2 ) (from calibDB BCAL/bcal_parms)
double BCAL_TIMEDIFFCOEFB        = 0.0; // 0.00 * sqrt( 2 ) (from calibDB BCAL/bcal_parms)
double BCAL_TWO_HIT_RESOL        = 0.0; // 50. (from calibDB BCAL/bcal_parms)
double BCAL_mevPerPE             = 0.31; // Energy corresponding to one pixel firing in MeV

int BCAL_NUM_MODULES = 48;
int BCAL_NUM_LAYERS = 4;
int BCAL_NUM_SECTORS = 4;

double BCAL_BASE_TIME_OFFSET     = 0; // -100.0 (from calibDB BCAL/base_time_offset)
double BCAL_TDC_BASE_TIME_OFFSET = 0; // -100.0 (from calibDB BCAL/base_time_offset)

vector<vector<double> > attenuation_parameters; // Avg. of 525 (from calibDB BCAL/attenuation_parameters)
vector<double> effective_velocities; // 16.75 (from calibDB BCAL/effective_velocities)

double BCAL_ADC_THRESHOLD_MEV    = 2.2;  // MeV (To be updated/improved)
double BCAL_FADC_TIME_RESOLUTION = 0.3;  // ns (To be updated/improved)
double BCAL_TDC_TIME_RESOLUTION  = 0.3;  // ns (To be updated/improved)
double BCAL_MEV_PER_ADC_COUNT    = 0.029;  // MeV per integrated ADC count (based on Spring 2015 calibrations)
double BCAL_NS_PER_ADC_COUNT     = 0.0;  // 0.0625 ns per ADC count (from calibDB BCAL/digi_scales)
double BCAL_NS_PER_TDC_COUNT     = 0.0;  // 0.0559 ns per TDC count (from calibDB BCAL/digi_scales)

// BCAL flags
bool NO_T_SMEAR = false;
bool NO_DARK_PULSES = false;
bool NO_SAMPLING_FLUCTUATIONS = false;
bool NO_SAMPLING_FLOOR_TERM = false;
bool NO_POISSON_STATISTICS = false;

double FTOF_BAR_THRESHOLD    = 0.0;
double STC_PADDLE_THRESHOLD  = 0.0;
double PSC_THRESHOLD         = 0.0;
double PS_THRESHOLD          = 0.0;

double TAGM_TSIGMA = 0.200;        // ns
double TAGH_TSIGMA = 0.350;        // ns
double TAGM_FADC_TSIGMA = 0.350;   // ns
double TAGH_FADC_TSIGMA = 0.450;   // ns
double TAGM_NPIX_PER_GEV = 1.e5;
double TAGH_NPE_PER_GEV = 5.e5;

// The following are place holders to be filled in with realistic numbers later
double PS_SIGMA = 0.200; // ns
double PSC_SIGMA = 0.200; //ns
double PS_NPIX_PER_GEV = 1.e5;
double PSC_PHOTONS_PERMEV = 5.e5;

double FCAL_PHOT_STAT_COEF   = 0.0; //0.035;
double FCAL_BLOCK_THRESHOLD  = 0.0; //20.0*k_MeV;

double CDC_TDRIFT_SIGMA      = 0.0; // 150.0/55.0*1E-9 seconds
double CDC_TIME_WINDOW       = 0.0; // 1000.0E-9 seconds
double CDC_PEDESTAL_SIGMA    = 0.0; // in fC
double CDC_THRESHOLD_FACTOR  = 0.0; // number of pedestal sigmas for determining sparcification threshold

double FDC_TDRIFT_SIGMA      = 0.0; // 200.0/55.0*1.0E-9 seconds
double FDC_CATHODE_SIGMA     = 0.0; // 150.0 microns
double FDC_PED_NOISE         = 0.0; // in pC calculated in SmearFDC
double FDC_THRESHOLD_FACTOR  = 0.0; // number of pedestal sigmas for determining sparcification threshold
double FDC_HIT_DROP_FRACTION = 0.0; // 1000.0E-9
double FDC_TIME_WINDOW       = 0.0; // 1000.0E-9 in seconds
double FDC_THRESH_KEV        = 0.0; // fdc anode discriminator threshold

double START_SIGMA           = 0.0; // 300ps
double START_PHOTONS_PERMEV  = 0.0; // used to be 8000 should be more like 200

// TOF parameters will be read from data base later
double TOF_SIGMA = 100.*k_psec;
double TOF_PHOTONS_PERMEV = 400.;

double FMWPC_TSIGMA = 10.0;  // ns
double FMWPC_ASIGMA = 0.5E-6;
double FMWPC_THRESHOLD = 0.0;

double TRIGGER_LOOKBACK_TIME = -100; // ns

bool DROP_TRUTH_HITS=false;

// CDC geometry parameters (for noise)
vector<unsigned int> NCDC_STRAWS;
vector<double> CDC_RING_RADIUS;

// FDC geometry and rate parameters (for noise)
vector<double> FDC_LAYER_Z;
double FDC_RATE_COEFFICIENT;

#include <JANA/JCalibrationFile.h>
using namespace jana;
static JCalibration *jcalib=NULL;

// histogram
pthread_mutex_t root_mutex = PTHREAD_MUTEX_INITIALIZER;
TH2F *fdc_drift_time_smear_hist;
TH2F *fdc_drift_dist_smear_hist;
TH2F *fdc_drift_time;
TH1F *fdc_cathode_charge;
TH2F *cdc_drift_time;
TH1F *cdc_charge;
TH1F *fdc_anode_mult;
TH2F *cdc_drift_smear;

//-----------
// main
//-----------
int main(int narg,char* argv[])
{
#if ! USE_JANA
   // Set up to catch SIGINTs for graceful exits
   signal(SIGINT,ctrlCHandleMCSmear);
#endif
   ParseCommandLineArguments(narg, argv);

   // Create DApplication object and use it to create JCalibration object
   DApplication dapp(narg, argv);
   dapp.AddFactoryGenerator(new JFactoryGenerator_ThreadCancelHandler());
   jcalib = dapp.GetJCalibration(1);

   // Make sure jcalib is set
   if(!jcalib){
     _DBG_<<"ERROR - jcalib not set!"<<endl;
     _DBG_<<"ERROR - Exiting ..."<<endl;
     return 0;
   }
   
   // get the TOF parameters
   {
     cout<<"get TOF/tof_parms parameters from calibDB"<<endl;
     map<string, double> tofparms;
     jcalib->Get("TOF/tof_parms", tofparms);
     TOF_SIGMA =  tofparms["TOF_SIGMA"];
     TOF_PHOTONS_PERMEV =  tofparms["TOF_PHOTONS_PERMEV"];
   }

   // get the BCAL parameters
   {
     cout<<"get BCAL/bcal_parms parameters from calibDB"<<endl;
     map<string, double> bcalparms;
     jcalib->Get("BCAL/bcal_parms", bcalparms);
     BCAL_DARKRATE_GHZ         =  bcalparms["BCAL_DARKRATE_GHZ"];
     BCAL_SIGMA_SIG_RELATIVE   = bcalparms["BCAL_SIGMA_SIG_RELATIVE"];
     BCAL_SIGMA_PED_RELATIVE   = bcalparms["BCAL_SIGMA_PED_RELATIVE"];
     BCAL_SIPM_GAIN_VARIATION   = bcalparms["BCAL_SIPM_GAIN_VARIATION"];
     BCAL_XTALK_FRACT         = bcalparms["BCAL_XTALK_FRACT"];
     BCAL_INTWINDOW_NS         = bcalparms["BCAL_INTWINDOW_NS"];
     BCAL_DEVICEPDE         = bcalparms["BCAL_DEVICEPDE"];
     BCAL_SAMPLING_FRACT      = bcalparms["BCAL_SAMPLING_FRACT"];
     BCAL_AVG_DARK_DIGI_VALS_PER_EVENT      = bcalparms["BCAL_AVG_DARK_DIGI_VALS_PER_EVENT"];
     BCAL_PHOTONSPERSIDEPERMEV_INFIBER = bcalparms["BCAL_PHOTONSPERSIDEPERMEV_INFIBER"];
     BCAL_SAMPLINGCOEFA = bcalparms["BCAL_SAMPLINGCOEFA"];
     BCAL_SAMPLINGCOEFB = bcalparms["BCAL_SAMPLINGCOEFB"];
     BCAL_TIMEDIFFCOEFA = bcalparms["BCAL_TIMEDIFFCOEFA"];
     BCAL_TIMEDIFFCOEFB = bcalparms["BCAL_TIMEDIFFCOEFB"];
     BCAL_TWO_HIT_RESOL = bcalparms["BCAL_TWO_HIT_RESOL"];
   }

   {
     cout<<"get BCAL/attenuation_parameters from calibDB"<<endl;
     vector< vector<double> > in_atten_parameters;
     jcalib->Get("BCAL/attenuation_parameters", in_atten_parameters);
     attenuation_parameters.clear();
     int channel = 0;
     for (int module=1; module<=BCAL_NUM_MODULES; module++) {
   	  for (int layer=1; layer<=BCAL_NUM_LAYERS; layer++) {
   		for (int sector=1; sector<=BCAL_NUM_SECTORS; sector++) {
			//int cell_id = GetCalibIndex(module,layer,sector);

			vector<double> new_params(3,0.);
			//attenuation_parameters[cell_id][0] = in_atten_parameters[channel][0];
			//attenuation_parameters[cell_id][1] = in_atten_parameters[channel][1];
			//attenuation_parameters[cell_id][2] = in_atten_parameters[channel][2];
			// hack to workaround odd CCDB behavior
			//attenuation_parameters[cell_id][0] = in_atten_parameters[channel][1];
			//attenuation_parameters[cell_id][1] = in_atten_parameters[channel][2];
			//attenuation_parameters[cell_id][2] = in_atten_parameters[channel][0];
			 
			new_params[0] = in_atten_parameters[channel][0];
			new_params[1] = in_atten_parameters[channel][1];
			new_params[2] = in_atten_parameters[channel][2];
			attenuation_parameters.push_back( new_params );

			channel++;
		}
	  }
     }
   }

   {
     cout<<"get BCAL/effective_velocities parameters from calibDB"<<endl;
     vector <double> effective_velocities_temp;
     jcalib->Get("BCAL/effective_velocities", effective_velocities_temp);
     for (unsigned int i = 0; i < effective_velocities_temp.size(); i++){
       effective_velocities.push_back(effective_velocities_temp.at(i));
     }
   }

   {
     cout<<"get BCAL/digi_scales parameters from calibDB"<<endl;
     map<string, double> bcaldigiscales;
     jcalib->Get("BCAL/digi_scales", bcaldigiscales);
     BCAL_NS_PER_ADC_COUNT = bcaldigiscales["BCAL_ADC_TSCALE"];
     BCAL_NS_PER_TDC_COUNT = bcaldigiscales["BCAL_TDC_SCALE"];
   }

   {
     cout<<"get BCAL/base_time_offset parameters from calibDB"<<endl;
     map<string, double> bcaltimeoffsets;
     jcalib->Get("BCAL/base_time_offset", bcaltimeoffsets);
     BCAL_BASE_TIME_OFFSET = bcaltimeoffsets["BCAL_BASE_TIME_OFFSET"];
     BCAL_TDC_BASE_TIME_OFFSET = bcaltimeoffsets["BCAL_TDC_BASE_TIME_OFFSET"];
   }

   {
     cout<<"get FCAL/fcal_parms parameters from calibDB"<<endl;
     map<string, double> fcalparms;
     jcalib->Get("FCAL/fcal_parms", fcalparms);
     if (FCAL_PHOT_STAT_COEF == 0.0)
       FCAL_PHOT_STAT_COEF   = fcalparms["FCAL_PHOT_STAT_COEF"]; 
     if (FCAL_BLOCK_THRESHOLD == 0.0)
       FCAL_BLOCK_THRESHOLD  = fcalparms["FCAL_BLOCK_THRESHOLD"];
   }
   {
     cout<<"get CDC/cdc_parms parameters from calibDB"<<endl;
     map<string, double> cdcparms;
     jcalib->Get("CDC/cdc_parms", cdcparms);
     if (CDC_TDRIFT_SIGMA == 0.0)
       CDC_TDRIFT_SIGMA   = cdcparms["CDC_TDRIFT_SIGMA"]; 
     if (CDC_TIME_WINDOW == 0.0)
       CDC_TIME_WINDOW    = cdcparms["CDC_TIME_WINDOW"];
     if (CDC_PEDESTAL_SIGMA == 0.0)
       CDC_PEDESTAL_SIGMA = cdcparms["CDC_PEDESTAL_SIGMA"]; 
     if (CDC_THRESHOLD_FACTOR == 0.0)
       CDC_THRESHOLD_FACTOR = cdcparms["CDC_THRESHOLD_FACTOR"];
   }

   {
     cout<<"get FDC/fdc_parms parameters from calibDB"<<endl;
     map<string, double> fdcparms;
     jcalib->Get("FDC/fdc_parms", fdcparms);

     if (FDC_TDRIFT_SIGMA == 0.0)
       FDC_TDRIFT_SIGMA      = fdcparms["FDC_TDRIFT_SIGMA"];
     if (FDC_CATHODE_SIGMA ==0.0)
       FDC_CATHODE_SIGMA     = fdcparms["FDC_CATHODE_SIGMA"];
     if (FDC_THRESHOLD_FACTOR == 0.0)
       FDC_THRESHOLD_FACTOR = fdcparms["FDC_THRESHOLD_FACTOR"];
     FDC_PED_NOISE         = fdcparms["FDC_PED_NOISE"];

     if (FDC_TIME_WINDOW == 0.0)
       FDC_TIME_WINDOW       = fdcparms["FDC_TIME_WINDOW"];

     if (FDC_HIT_DROP_FRACTION == 0.0)
       FDC_HIT_DROP_FRACTION = fdcparms["FDC_HIT_DROP_FRACTION"];  
     if (FDC_THRESH_KEV == 0.0)
       FDC_THRESH_KEV = fdcparms["FDC_THRESH_KEV"]; 
   }

   {
     cout<<"get START_COUNTER/start_parms parameters from calibDB"<<endl;
     map<string, double> startparms;
     jcalib->Get("START_COUNTER/start_parms", startparms);

     START_SIGMA = startparms["START_SIGMA"] ;
     START_PHOTONS_PERMEV = startparms["START_PHOTONS_PERMEV"];

   }

   // hist file
   TFile *hfile = new TFile("smear.root","RECREATE","smearing histograms");
   fdc_drift_time_smear_hist=new TH2F("fdc_drift_time_smear_hist","Drift time smearing for FDC",
                  300,0.0,0.6,400,-200,200);
   fdc_drift_dist_smear_hist=new TH2F("fdc_drift_dist_smear_hist","Drift distance smearing for FDC",
                  100,0.0,0.6,400,-0.5,0.5);
   double tmax=TRIGGER_LOOKBACK_TIME+FDC_TIME_WINDOW;
   int num_time_bins=int(FDC_TIME_WINDOW);
   fdc_drift_time=new TH2F("fdc_drift_time","FDC drift distance vs. time",num_time_bins,TRIGGER_LOOKBACK_TIME,tmax,100,0,1.);
   
   fdc_anode_mult = new TH1F("fdc_anode_mult","wire hit multiplicity",20,-0.5,19.5);
   fdc_cathode_charge = new TH1F("fdc_cathode_charge","charge on strips",1000,0,1000);

   tmax=TRIGGER_LOOKBACK_TIME+CDC_TIME_WINDOW;
   num_time_bins=int(CDC_TIME_WINDOW);
   cdc_drift_time = new TH2F("cdc_drift_time","CDC drift distance vs time",num_time_bins,TRIGGER_LOOKBACK_TIME,tmax,80,0.,0.8);

   cdc_drift_smear = new TH2F("cdc_drift_smear","CDC drift smearing",
               100,0.0,800.0,100,-0.1,0.1);
   
   cdc_charge  = new TH1F("cdc_charge","Measured charge in straw",1000,-10e3,40e3);



#if ! USE_JANA
   cout << " input file: " << INFILENAME << endl;
   cout << " output file: " << OUTFILENAME << endl;
   
   // Open Input file
   ifstream ifs(INFILENAME);
   if (!ifs.is_open()) {
      cout << " Error opening input file \"" << INFILENAME << "\"!" << endl;
      exit(-1);
   }
   hddm_s::istream fin(ifs);
   
   // Output file
   ofstream ofs(OUTFILENAME);
   if (!ofs.is_open()){
      cout << " Error opening output file \"" << OUTFILENAME << "\"!" << endl;
      exit(-1);
   }
   hddm_s::HDDM fout(ofs);
   
   // Loop over events in input file
   hddm_s::HDDM *record;
   int NEvents = 0;
   time_t last_time = time(NULL);
   while (ifs->good()) {
      fin >> *record;
      NEvents++;
      time_t now = time(NULL);
      if(now != last_time){
         cout << "  " << NEvents << " events processed      \r";
         cout.flush();
         last_time = now;
      }
      
      // Smear values
      Smear(record);
      
      // Write event to output file
      *fout << *record;
      
      if (QUIT) break;
   }
   cout << endl;
   
   // close input and output files
   ifs.close();
   ofs.close();

   cout << " " << NEvents << " events read" << endl;

#else

   DGeometry *dgeom=dapp.GetDGeometry(1);
   
   // Get number of cdc wires per ring and the radii of each ring
   vector<vector<DCDCWire *> >cdcwires;
   dgeom->GetCDCWires(cdcwires);
   for (unsigned int i=0;i<cdcwires.size();i++) {
      NCDC_STRAWS.push_back(cdcwires[i].size());
      CDC_RING_RADIUS.push_back(cdcwires[i][0]->origin.Perp());
   }  
   // Get the FDC z positions for the wire planes
   dgeom->GetFDCZ(FDC_LAYER_Z);

   // Coefficient used to calculate FDCsingle wire rate. We calculate
   // it once here just to save calculating it for every wire in every event
   FDC_RATE_COEFFICIENT = exp(-log(4.0)/23.0)/2.0/log(24.0)*FDC_TIME_WINDOW/1000.0E-9;
   
   // Something is a little off in my calculation above so I scale it down via
   // an emprical factor:
   FDC_RATE_COEFFICIENT *= 0.353;

   MyProcessor myproc;
   
   dapp.Run(&myproc);

#endif
   
   hfile->Write();
   hfile->Close();

   return 0;
}

//-----------
// ParseCommandLineArguments
//-----------
void ParseCommandLineArguments(int narg, char* argv[])
{
   bool warn_obsolete = false;

  for (int i=1; i<narg; i++) {
    char *ptr = argv[i];
    
    if (ptr[0] == '-') {
      switch(ptr[1]) {
      case 'h': Usage();                                   break;
      case 'o': OUTFILENAME = strdup(&ptr[2]);             break;
      case 'n': warn_obsolete=true;                        break;
      case 'N': ADD_NOISE=true;                            break;
      case 's': SMEAR_HITS=false;                          break;
      case 'i': IGNORE_SEEDS=true;                         break;
      case 'r': SetSeeds(&ptr[2]);                         break;
      case 'u': CDC_TDRIFT_SIGMA=atof(&ptr[2])*1.0E-9;     break;
      case 't': CDC_TIME_WINDOW=atof(&ptr[2])*1.0E-9;      break;
      case 'U': FDC_TDRIFT_SIGMA=atof(&ptr[2])*1.0E-9;     break;
      case 'C': FDC_CATHODE_SIGMA=atof(&ptr[2])*1.0E-6;    break;
      case 'T': FDC_TIME_WINDOW=atof(&ptr[2])*1.0E-9;      break;
      case 'e': FDC_ELOSS_OFF = true;                      break;
      case 'E': CDC_PEDESTAL_SIGMA = atof(&ptr[2])*k_keV;  break;
      case 'd': DROP_TRUTH_HITS=true;                      break;
      case 'p': FCAL_PHOT_STAT_COEF = atof(&ptr[2]);       break;
      case 'b': FCAL_BLOCK_THRESHOLD = atof(&ptr[2])*k_MeV; break;
      case 'B': SMEAR_BCAL = false;                        break;
      case 'V': BCAL_ADC_THRESHOLD_MEV = atof(&ptr[2]);    break;
      case 'X': BCAL_FADC_TIME_RESOLUTION = atof(&ptr[2]); break;
      case 'G': NO_T_SMEAR = true;                         break;
      case 'H': NO_DARK_PULSES = true;                     break;
      case 'K': NO_SAMPLING_FLUCTUATIONS = true;           break;
      case 'L': NO_SAMPLING_FLOOR_TERM = true;             break;
      case 'M': NO_POISSON_STATISTICS = true;              break;
      case 'f': TOF_SIGMA= atof(&ptr[2])*k_psec;           break;
      case 'S': START_SIGMA= atof(&ptr[2])*k_psec;         break;
      }
    }
    else {
      INFILENAME = argv[i];
    }
  }

   if (!INFILENAME){
      cout << endl << "You must enter a filename!" << endl << endl;
      Usage();
   }
   
   if (warn_obsolete) {
      cout << endl;
      cout << "WARNING: Use of the \"-n\" option is obsolete. Random noise"
           << endl;
      cout << "         hits are disabled by default now. To turn them back"
           << endl;
      cout << "         on use the \"-N\" option." << endl;
      cout << endl;
   }
   
   // Generate output filename based on input filename
   if (OUTFILENAME == NULL) {
      char *ptr, *path_stripped;
      path_stripped = ptr = strdup(INFILENAME);
      while((ptr = strstr(ptr, "/")))path_stripped = ++ptr;
      ptr = strstr(path_stripped, ".hddm");
      if(ptr)*ptr=0;
      char str[256];
      sprintf(str, "%s_%ssmeared.hddm", path_stripped, ADD_NOISE ? "n":"");
      OUTFILENAME = strdup(str);
   }
   
   cout << "BCAL values will " <<  (SMEAR_BCAL ? "":"not")  << " be smeared"

        << endl;
   cout << "BCAL values will " <<  (SMEAR_BCAL ? "":"not")  << " be added"
        << endl;
}


//-----------
// Usage
//-----------
void Usage(void)
{
   cout << endl << "Usage:" << endl;
   cout << "     mcsmear [options] file.hddm" << endl;
   cout << endl;
   cout << " Read the given, Geant-produced HDDM file as input and smear" << endl;
   cout << "the truth values for \"hit\" data before writing out to a" << endl;
   cout << "separate file. The truth values for the thrown particles are" << endl;
   cout << "not changed. Noise hits can also be added using the -n option." << endl;
   cout << "Note that all smearing is done using Gaussians, with the " << endl;
   cout << "sigmas configurable with the options below." << endl;
   cout << endl;
   cout << "  options:" << endl;
   cout << "    -ofname  Write output to a file named \"fname\" (default auto-generate name)" << endl;
   cout << "    -N       Add random background hits to CDC and FDC (default is not to add)" << endl;
   cout << "    -s       Don't smear real hits (see -B for BCAL, default is to smear)" << endl;
   cout << "    -i       Ignore random number seeds found in input HDDM file" << endl;
   cout << "    -r\"s1 s2 s3\" Set initial random number seeds" << endl;
   cout << "    -u#      Sigma CDC anode drift time in ns (def:" << CDC_TDRIFT_SIGMA*1.0E9 << "ns)" << endl;
   cout << "             (NOTE: this is only used if -y is also specified!)" << endl;
   cout << "    -y       Do NOT apply drift distance dependence error to" << endl;
   cout << "             CDC (default is to apply)" << endl;
   cout << "    -Y       Apply constant sigma smearing for FDC drift time. "  << endl;
   cout << "             Default is to use a drift-distance dependent parameterization."  << endl;
   cout << "    -t#      CDC time window for background hits in ns (def:" << CDC_TIME_WINDOW*1.0E9 << "ns)" << endl;
   cout << "    -U#      Sigma FDC anode drift time in ns (def:" << FDC_TDRIFT_SIGMA*1.0E9 << "ns)" << endl;
   cout << "    -C#      Sigma FDC cathode strips in microns (def:" << FDC_TDRIFT_SIGMA << "ns)" << endl;
   cout << "    -T#      FDC time window for background hits in ns (def:" << FDC_TIME_WINDOW*1.0E9 << "ns)" << endl;
   cout << "    -e       hdgeant was run with LOSS=0 so scale the FDC cathode" << endl;
   cout << "             pedestal noise (def:false)" << endl;
   cout << "    -d       Drop truth hits (default: keep truth hits)" << endl;
   cout << "    -p#      FCAL photo-statistics smearing factor in GeV^3/2 (def:" << FCAL_PHOT_STAT_COEF << ")" << endl;
   cout << "    -b#      FCAL single block threshold in MeV (def:" << FCAL_BLOCK_THRESHOLD/k_MeV << ")" << endl;
   cout << "    -B       Don't process BCAL hits at all (def. process)" << endl;
   cout << "    -Vthresh BCAL ADC threshold (def. " << BCAL_ADC_THRESHOLD_MEV << " MeV)" << endl;
   cout << "    -Xsigma  BCAL fADC time resolution (def. " << BCAL_FADC_TIME_RESOLUTION << " ns)" << endl;
   cout << "    -G       Don't smear BCAL times (def. smear)" << endl;
   cout << "    -H       Don't add BCAL dark hits (def. add)" << endl;
   cout << "    -K       Don't apply BCAL sampling fluctuations (def. apply)" << endl;
   cout << "    -L       Don't apply BCAL sampling floor term (def. apply)" << endl;
   cout << "    -M       Don't apply BCAL Poisson statistics (def. apply)" << endl;
   cout << "    -f#      TOF sigma in psec (def: " <<  TOF_SIGMA/k_psec << ")" << endl;
   cout << "    -h       Print this usage statement." << endl;
   cout << endl;
   cout << " Example:" << endl;
   cout << endl;
   cout << "     mcsmear -u3.5 -t500 hdgeant.hddm" << endl;
   cout << endl;
   cout << " This will produce a file named hdgeant_nsmeared.hddm that" << endl;
   cout << " includes the hit information from the input file hdgeant.hddm" << endl;
   cout << " but with the FDC and CDC hits smeared out. The CDC hits will" << endl;
   cout << " have their drift times smeared via a gaussian with a 3.5ns width" << endl;
   cout << " while the FDC will be smeared using the default values." << endl;
   cout << " In addition, background hits will be added, the exact number of" << endl;
   cout << " of which are determined by the time windows specified for the" << endl;
   cout << " CDC and FDC. In this examplem the CDC time window was explicitly" << endl;
   cout << " set to 500 ns." << endl;
   cout << endl;

   exit(0);
}

#if ! USE_JANA
//-----------------------------------------------------------------
// ctrlCHandleMCSmear
//-----------------------------------------------------------------
void ctrlCHandleMCSmear(int x)
{
   QUIT++;
   cerr << endl << "SIGINT received (" << QUIT << ")....." << endl;
}
#endif
