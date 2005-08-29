// $Id$
//
// Created June 22, 2005  David Lawrence

#include <iostream>
#include <iomanip>
using namespace std;

#include <signal.h>
#include <time.h>

#include "hddm_s.h"

void Smear(s_HDDM_t *hddm_s);
void ParseCommandLineArguments(int narg, char* argv[]);
void Usage(void);
void ctrlCHandle(int x);

char *INFILENAME = NULL;
char *OUTFILENAME = NULL;
int QUIT = 0;

extern bool ADD_NOISE;
extern float CDC_R;
extern float CDC_Z_SIGMA;
extern float CDC_AVG_NOISE_HITS;
extern float FDC_R;
extern float FDC_AVG_NOISE_HITS;


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
				case 'n': ADD_NOISE=true;											break;
				case 'r': CDC_R=atof(&ptr[2]);									break;
				case 'z': CDC_Z_SIGMA=atof(&ptr[2]);							break;
				case 'j': CDC_AVG_NOISE_HITS=3240.0*atof(&ptr[2])/100.0;	break;
				case 'R': FDC_R=atof(&ptr[2]);									break;
				case 'J': CDC_AVG_NOISE_HITS=2856.0*atof(&ptr[2])/100.0;	break;
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
	sprintf(str, "%s_smeared.hddm", path_stripped);
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
	cout<<"    -n       Add noise hits to CDC and FDC (def:"<<ADD_NOISE<<")"<<endl;
	cout<<"    -r       Sigma CDC r-direction (def:"<<CDC_R<<"cm)"<<endl;
	cout<<"    -z       Sigma CDC z-direction (def:"<<CDC_Z_SIGMA<<"cm)"<<endl;
	cout<<"    -j       CDC occupancy r-direction (def:"<<100.0*CDC_AVG_NOISE_HITS/3240.0<<"%)"<<endl;
	cout<<"    -R       Sigma FDC r-direction (def:"<<FDC_R<<"cm)"<<endl;
	cout<<"    -J       FDC occupancy r-direction (def:"<<100.0*FDC_AVG_NOISE_HITS/2856.0<<"%)"<<endl;
	cout<<"    -h       Print this usage statement."<<endl;
	cout<<endl;
	cout<<" Example:"<<endl;
	cout<<endl;
	cout<<"     mcsmear -n -r1.2 -J5.5"<<endl;
	cout<<endl;
	cout<<" This will turn on CDC and FDC noise hits, set the r-direction"<<endl;
	cout<<" sigma for CDC to 1.2cm, and set the noise hit occupancy level"<<endl;
	cout<<" to 5.5%. Note that the occupancy values apply only to the"<<endl;
	cout<<" noise hits so the actual occupancy will be higher due to the"<<endl;
	cout<<" presence of the real hits."<<endl;
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
