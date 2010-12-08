// $Id$
//
// Created June 22, 2005  David Lawrence

#include <iostream>
#include <iomanip>
using namespace std;

#include <TF1.h>
#include <TFile.h>
#include <TH2.h>
#include <TH1.h>

#include <signal.h>
#include <time.h>

#include "units.h"
#include "HDDM/hddm_s.h"

void Smear(s_HDDM_t *hddm_s);
void ParseCommandLineArguments(int narg, char* argv[]);
void Usage(void);
void ctrlCHandle(int x);

char *INFILENAME = NULL;
char *OUTFILENAME = NULL;
int QUIT = 0;

extern bool ADD_NOISE;
extern bool SMEAR_HITS;
extern bool CDC_USE_PARAMETERIZED_SIGMA;
extern bool FDC_USE_PARAMETERIZED_SIGMA;
extern double CDC_TDRIFT_SIGMA;
extern double CDC_TIME_WINDOW;
extern double CDC_PEDESTAL_SIGMA;
extern double FDC_TDRIFT_SIGMA;
extern double FDC_CATHODE_SIGMA;
extern double FDC_PED_NOISE;
extern bool FDC_ELOSS_OFF;
extern double FDC_TIME_WINDOW;
extern double FDC_HIT_DROP_FRACTION;
extern double FCAL_PHOT_STAT_COEF;
extern double FCAL_BLOCK_THRESHOLD;
extern double TOF_SIGMA;
extern double START_SIGMA;

vector<vector<float> >fdc_smear_parms;
TF1 *fdc_smear_function;

#include <JANA/JCalibrationFile.h>
using namespace jana;
static JCalibration *jcalib=NULL;

// histogram
TH2F *fdc_drift_time_smear_hist;
TH2F *fdc_drift_dist_smear_hist;
TH2F *fdc_drift_time;

//-----------
// main
//-----------
int main(int narg,char* argv[])
{
	// Set up to catch SIGINTs for graceful exits
	signal(SIGINT,ctrlCHandle);

	ParseCommandLineArguments(narg, argv);
	
	// hist file
	TFile *hfile = new TFile("smear.root","RECREATE","smearing histograms");
	fdc_drift_time_smear_hist=new TH2F("fdc_drift_time_smear_hist","Drift time smearing for FDC",
					   300,0.0,0.6,400,-200,200);
	fdc_drift_dist_smear_hist=new TH2F("fdc_drift_dist_smear_hist","Drift distance smearing for FDC",
					   100,0.0,0.6,400,-0.5,0.5);
	fdc_drift_time=new TH2F("fdc_drift_time","FDC drift distance vs. time",100,-20,380,100,0,1.);
	

	// Create a JCalibration object using the JANA_CALIB_URL environment variable
	// Right now, we hardwire this to use JCalibrationFile.
	const char *url = getenv("JANA_CALIB_URL");
	if(!url){
		_DBG_<<"JANA_CALIB_URL environment not set."<<endl;
		exit(-1);
	}
	jcalib = new JCalibrationFile(url, 1, "");
	// Make sure jcalib is set
	if(!jcalib){
	  _DBG_<<"ERROR - jcalib not set!"<<endl;
	  _DBG_<<"ERROR - Exiting ..."<<endl;
	  return 0;
	}
	
	if (FDC_USE_PARAMETERIZED_SIGMA==true){ // Get the drift time smearing parameters from the calib DB
	  fdc_smear_function=new TF1("f1","gaus(0)+gaus(3)+gaus(6)",-0.5,0.5);
	  
	  vector< map<string, float> > tvals;
	  jcalib->Get("FDC/drift_smear_parms", tvals);
	  // Notify user
	  cout<<"Read "<<tvals.size()<<" values from FDC/drift_smear_parms in calibDB"
	      <<endl;
	  cout << "Columns:  " ;
	  map<string,float>::iterator iter;
	  for(iter=tvals[0].begin(); iter!=tvals[0].end(); iter++)cout<<iter->first<<" ";
	  cout<<endl;

	  for(unsigned int i=0; i<tvals.size(); i++){
	    map<string, float> &row = tvals[i];
	    vector<float>dummy;
	    dummy.push_back(row["h0"]);
	    dummy.push_back(row["m0"]);
	    dummy.push_back(row["s0"]);
	    dummy.push_back(row["h1"]);
	    dummy.push_back(row["m1"]);
	    dummy.push_back(row["s1"]);
	    dummy.push_back(row["h2"]);
	    dummy.push_back(row["m2"]);
	    dummy.push_back(row["s2"]);
	    fdc_smear_parms.push_back(dummy);
	    dummy.clear();
	  }
	}


	cout<<" input file: "<<INFILENAME<<endl;
	cout<<" output file: "<<OUTFILENAME<<endl;
	
	// Open Input file
	s_iostream_t *fin = open_s_HDDM(INFILENAME);
	if(!fin){
		cout<<" Error opening input file \""<<INFILENAME<<"\"!"<<endl;
		exit(-1);
	}
	
	// Output file
	s_iostream_t *fout = init_s_HDDM(OUTFILENAME);
	if(!fout){
		cout<<" Error opening output file \""<<OUTFILENAME<<"\"!"<<endl;
		exit(-1);
	}
	
	// Loop over events in input file
	s_HDDM_t *hddm_s;
	int NEvents = 0;
	time_t last_time = time(NULL);
	while((hddm_s = read_s_HDDM(fin))){
		NEvents++;
		time_t now = time(NULL);
		if(now != last_time){
			cout<<"  "<<NEvents<<" events processed      \r";cout.flush();
			last_time = now;
		}
		
		// Smear values
		Smear(hddm_s);
		
		// Write event to output file
		flush_s_HDDM(hddm_s, fout);
		
		if(QUIT)break;
	}
	
	// close input and output files
	close_s_HDDM(fin);
	close_s_HDDM(fout);

	cout<<" "<<NEvents<<" events read"<<endl;
	
	hfile->Write();

	return 0;
}

//-----------
// ParseCommandLineArguments
//-----------
void ParseCommandLineArguments(int narg, char* argv[])
{
	bool warn_obsolete = false;

  for(int i=1; i<narg; i++){
    char *ptr = argv[i];
    
    if(ptr[0] == '-'){
      switch(ptr[1]){
      case 'h': Usage();													break;
      case 'n': warn_obsolete=true;										break;
      case 'N': ADD_NOISE=true;											break;
      case 's': SMEAR_HITS=false;										break;
      case 'u': CDC_TDRIFT_SIGMA=atof(&ptr[2])*1.0E-9;			break;
      case 'y': CDC_USE_PARAMETERIZED_SIGMA=false;					break; 
      case 'Y': FDC_USE_PARAMETERIZED_SIGMA=false;					break;
      case 't': CDC_TIME_WINDOW=atof(&ptr[2])*1.0E-9;				break;
      case 'U': FDC_TDRIFT_SIGMA=atof(&ptr[2])*1.0E-9;			break;
      case 'C': FDC_CATHODE_SIGMA=atof(&ptr[2])*1.0E-6;			break;
      case 'T': FDC_TIME_WINDOW=atof(&ptr[2])*1.0E-9;				break;
      case 'e': FDC_ELOSS_OFF = true;	break;
      case 'E': CDC_PEDESTAL_SIGMA = atof(&ptr[2])*k_keV; break;
      case 'd': FDC_HIT_DROP_FRACTION=atof(&ptr[2]);				break;
      case 'p': FCAL_PHOT_STAT_COEF = atof(&ptr[2]);				break;
      case 'b': FCAL_BLOCK_THRESHOLD = atof(&ptr[2])*k_MeV;		break;
      case 'f': TOF_SIGMA= atof(&ptr[2])*k_psec; break;
      case 'S': START_SIGMA= atof(&ptr[2])*k_psec; break;
      }
    }else{
      INFILENAME = argv[i];
    }
  }

	if(!INFILENAME){
		cout<<endl<<"You must enter a filename!"<<endl<<endl;
		Usage();
	}
	
	if(warn_obsolete){
		cout<<endl;
		cout<<"WARNING: Use of the \"-n\" option is obsolete. Random noise"<<endl;
		cout<<"         hits are disabled by default now. To turn them back"<<endl;
		cout<<"         on use the \"-N\" option."<<endl;
		cout<<endl;
	}
	
	// Generate output filename based on input filename
	char *ptr, *path_stripped;
	path_stripped = ptr = strdup(INFILENAME);
	while((ptr = strstr(ptr, "/")))path_stripped = ++ptr;
	ptr = strstr(path_stripped, ".hddm");
	if(ptr)*ptr=0;
	char str[256];
	sprintf(str, "%s_%ssmeared.hddm", path_stripped, ADD_NOISE ? "n":"");
	OUTFILENAME = strdup(str);
}


//-----------
// Usage
//-----------
void Usage(void)
{
	cout<<endl<<"Usage:"<<endl;
	cout<<"     mcsmear [options] file.hddm"<<endl;
	cout<<endl;
	cout<<" Read the given, Geant-produced HDDM file as input and smear"<<endl;
	cout<<"the truth values for \"hit\" data before writing out to a"<<endl;
	cout<<"separate file. The truth values for the thrown particles are"<<endl;
	cout<<"not changed. Noise hits can also be added using the -n option."<<endl;
	cout<<"Note that all smearing is done using Gaussians, with the "<<endl;
	cout<<"sigmas configurable with the options below."<<endl;
	cout<<endl;
	cout<<"  options:"<<endl;
	cout<<"    -N       Add random background hits to CDC and FDC (default is not to add)"<<endl;
	cout<<"    -s       Don't smear real hits (default is to smear)"<<endl;
	cout<<"    -u#      Sigma CDC anode drift time in ns (def:"<<CDC_TDRIFT_SIGMA*1.0E9<<"ns)"<<endl;
	cout<<"             (NOTE: this is only used if -y is also specified!)"<<endl;
	cout<<"    -y       Do NOT apply drift distance dependence error to"<<endl;
	cout<<"             CDC (default is to apply)"<<endl;
	cout<<"    -Y       Apply constant sigma smearing for FDC drift time. " <<endl;
	cout<<"             Default is to use a drift-distance dependent parameterization." <<endl;
	cout<<"    -t#      CDC time window for background hits in ns (def:"<<CDC_TIME_WINDOW*1.0E9<<"ns)"<<endl;
	cout<<"    -U#      Sigma FDC anode drift time in ns (def:"<<FDC_TDRIFT_SIGMA*1.0E9<<"ns)"<<endl;
	cout<<"    -C#      Sigma FDC cathode strips in microns (def:"<<FDC_TDRIFT_SIGMA<<"ns)"<<endl;
	cout<<"    -t#      FDC time window for background hits in ns (def:"<<FDC_TIME_WINDOW*1.0E9<<"ns)"<<endl;
	cout<<"    -e       hdgeant was run with LOSS=0 so scale the FDC cathode"<<endl;
	cout<<"             pedestal noise (def:false)"<<endl;
	cout<<"    -d#      Randomly drop this fraction of FDC hits (0=drop none  1=drop all)"<<endl;
	cout<<"             default is to drop none."<<endl;
	cout<<"    -p#      FCAL photo-statistics smearing factor in GeV^3/2 (def:"<<FCAL_PHOT_STAT_COEF<<")"<<endl;
	cout<<"    -b#      FCAL single block threshold in MeV (def:"<<FCAL_BLOCK_THRESHOLD/k_MeV<<")"<<endl;
	cout<<"    -f#      TOF sigma in psec (def: "<< TOF_SIGMA/k_psec<<")"<<endl;
	cout<<"    -h       Print this usage statement."<<endl;
	cout<<endl;
	cout<<" Example:"<<endl;
	cout<<endl;
	cout<<"     mcsmear -u3.5 -t500 hdgeant.hddm"<<endl;
	cout<<endl;
	cout<<" This will produce a file named hdgeant_nsmeared.hddm that"<<endl;
	cout<<" includes the hit information from the input file hdgeant.hddm"<<endl;
	cout<<" but with the FDC and CDC hits smeared out. The CDC hits will"<<endl;
	cout<<" have their drift times smeared via a gaussian with a 3.5ns width"<<endl;
	cout<<" while the FDC will be smeared using the default values."<<endl;
	cout<<" In addition, background hits will be added, the exact number of"<<endl;
	cout<<" of which are determined by the time windows specified for the"<<endl;
	cout<<" CDC and FDC. In this examplem the CDC time window was explicitly"<<endl;
	cout<<" set to 500 ns."<<endl;
	cout<<endl;

	exit(0);
}

//-----------------------------------------------------------------
// ctrlCHandle
//-----------------------------------------------------------------
void ctrlCHandle(int x)
{
	QUIT++;
	cerr<<endl<<"SIGINT received ("<<QUIT<<")....."<<endl;
}
