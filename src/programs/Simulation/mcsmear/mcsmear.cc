// $Id$
//
// Created June 22, 2005  David Lawrence

#include <iostream>
#include <iomanip>
using namespace std;

#include <signal.h>
#include <time.h>

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
extern double CDC_TDRIFT_SIGMA;
extern double CDC_TIME_WINDOW;
extern double FDC_TDRIFT_SIGMA;
extern double FDC_CATHODE_SIGMA;
extern double FDC_PED_NOISE;
extern double FDC_TIME_WINDOW;

//-----------
// main
//-----------
int main(int narg,char* argv[])
{
	// Set up to catch SIGINTs for graceful exits
	signal(SIGINT,ctrlCHandle);

	ParseCommandLineArguments(narg, argv);
	
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

	return 0;
}

//-----------
// ParseCommandLineArguments
//-----------
void ParseCommandLineArguments(int narg, char* argv[])
{

	for(int i=1; i<narg; i++){
		char *ptr = argv[i];
		
		if(ptr[0] == '-'){
			switch(ptr[1]){
				case 'h': Usage();													break;
				case 'n': ADD_NOISE=false;											break;
				case 's': SMEAR_HITS=false;										break;
				case 'u': CDC_TDRIFT_SIGMA=atof(&ptr[2])*1.0E-9;			break;
				case 't': CDC_TIME_WINDOW=atof(&ptr[2])*1.0E-9;				break;
				case 'U': FDC_TDRIFT_SIGMA=atof(&ptr[2])*1.0E-9;			break;
				case 'C': FDC_CATHODE_SIGMA=atof(&ptr[2])*1.0E-6;			break;
				case 'T': FDC_TIME_WINDOW=atof(&ptr[2])*1.0E-9;				break;
			}
		}else{
			INFILENAME = argv[i];
		}
	}

	if(!INFILENAME){
		cout<<endl<<"You must enter a filename!"<<endl<<endl;
		Usage();
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
	cout<<"    -n       Don't add background hits to CDC and FDC (default is to add)"<<endl;
	cout<<"    -s       Don't smear real hits (default is to smear)"<<endl;
	cout<<"    -u       Sigma CDC anode drift time in ns (def:"<<CDC_TDRIFT_SIGMA*1.0E9<<"ns)"<<endl;
	cout<<"    -t       CDC time window for background hits in ns (def:"<<CDC_TIME_WINDOW*1.0E9<<"ns)"<<endl;
	cout<<"    -U       Sigma FDC anode drift time in ns (def:"<<FDC_TDRIFT_SIGMA*1.0E9<<"ns)"<<endl;
	cout<<"    -C       Sigma FDC cathode strips in microns (def:"<<FDC_TDRIFT_SIGMA<<"ns)"<<endl;
	cout<<"    -t       FDC time window for background hits in ns (def:"<<FDC_TIME_WINDOW*1.0E9<<"ns)"<<endl;
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
