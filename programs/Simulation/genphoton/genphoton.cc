

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <vector>
using namespace std;

#include <TVector3.h>
#include <TLorentzVector.h>

unsigned int MAX_EVENTS=10000;
double P_MIN=0.100;
double P_MAX=6.000;
double PHI_MIN = 0.0;
double PHI_MAX = 2.0*M_PI;
double THETA_MIN = 0.0;
double THETA_MAX = M_PI;
bool IS_POSITIVE = true;

int RUN_NUMBER=100;
string OUTPUT_FILENAME="genphoton.ascii";

#define GAMMA_TYPE 1

#define _DBG_ cout<<__FILE__<<":"<<__LINE__<<" "
#define _DBG__ _DBG_<<endl

class photonX{
	public:
		double px,py,pz,E; // photonX
};


double SampleSin2Theta(void);
void ParseCommandLineArguments(int narg, char* argv[]);
void Usage(void);


//----------------------------
// main
//----------------------------
int main(int narg, char* argv[])
{

	// Parse the command line
	ParseCommandLineArguments(narg, argv);

	// Open file for output
	ofstream of(OUTPUT_FILENAME.c_str());
	if(!of.is_open()){
		cout<<"Unable to open \""<<OUTPUT_FILENAME<<"\" for writing!"<<endl;
		exit(-1);
	}

	// Seed random number generator
	srandom(time(NULL));	

	// Loop over events
	unsigned int nevents;
	for(nevents=1; nevents<=MAX_EVENTS; nevents++){
	
		vector<photonX> photonXs;
		photonX p;
	
		// Randomly sample the energy and angles of the pion
		double mom = (double)random()/(double)RAND_MAX*(P_MAX-P_MIN) + P_MIN;
		double phi = (double)random()/(double)RAND_MAX*(PHI_MAX-PHI_MIN) + PHI_MIN;
		double theta = (double)random()/(double)RAND_MAX*(THETA_MAX-THETA_MIN) + THETA_MIN;

		p.E = mom;	
		p.px = mom*sin(theta)*cos(phi);
		p.py = mom*sin(theta)*sin(phi);
		p.pz = mom*cos(theta);
		photonXs.push_back(p);
			
		// Write event to file
		unsigned int type = GAMMA_TYPE;
		of<<RUN_NUMBER<<" "<<nevents<<" "<<photonXs.size()<<endl;
		for(unsigned int j=0; j<photonXs.size(); j++){
			photonX &p = photonXs[j];
			unsigned int index = j+1;
			
			of<<index<<" "<<type<<" "<<0<<endl;
			of<<"   "<<0<<" "<<p.px<<" "<<p.py<<" "<<p.pz<<" "<<p.E<<endl;			
		}
		
		// Update screen
		if(nevents%1000 == 0){
			cout<<"  "<<nevents<<" events generated"<<endl;
		}

	}
	
	of.close();

	cout<<endl<<""<<nevents-1<<" events written to "<<OUTPUT_FILENAME<<endl;
	
	return 0;
}


//----------------------------
// ParseCommandLineArguments
//----------------------------
void ParseCommandLineArguments(int narg, char* argv[])
{
	if(narg==1)Usage();

	for(int i=1; i<narg; i++){
		string arg = argv[i];
		if(arg=="-h" || arg=="--help"){
			Usage();
		}else if(arg=="-M"){
			if(i==narg-1){cout<<"-M requires an argument!"<<endl; Usage();}
			MAX_EVENTS = atoi(argv[++i]);
		}else if(arg=="-Pmin"){
			if(i==narg-1){cout<<"-Pmin requires an argument!"<<endl; Usage();}
			P_MIN = atof(argv[++i]);
		}else if(arg=="-Pmax"){
			if(i==narg-1){cout<<"-Pmax requires an argument!"<<endl; Usage();}
			P_MAX = atof(argv[++i]);
		}else if(arg=="-Phimin"){
			if(i==narg-1){cout<<"-Phimin requires an argument!"<<endl; Usage();}
			PHI_MIN = M_PI/180.0*atof(argv[++i]);
		}else if(arg=="-Phimax"){
			if(i==narg-1){cout<<"-Phimax requires an argument!"<<endl; Usage();}
			PHI_MAX = M_PI/180.0*atof(argv[++i]);
		}else if(arg=="-Thetamin"){
			if(i==narg-1){cout<<"-Thetamin requires an argument!"<<endl; Usage();}
			THETA_MIN = M_PI/180.0*atof(argv[++i]);
		}else if(arg=="-Thetamax"){
			if(i==narg-1){cout<<"-Thetamax requires an argument!"<<endl; Usage();}
			THETA_MAX = M_PI/180.0*atof(argv[++i]);
		}else if(arg=="-o"){
			if(i==narg-1){cout<<"-o requires an argument!"<<endl; Usage();}
			OUTPUT_FILENAME = argv[++i];
		}else if(arg=="-n"){
			IS_POSITIVE = false;
		}
	}
	
	cout<<"---- genphoton will use the following settings: ----"<<endl;
	cout<<"MAX_EVENTS     = "<<MAX_EVENTS<<endl;
	cout<<"P_MIN          = "<<P_MIN<<endl;
	cout<<"P_MAX          = "<<P_MAX<<endl;
	cout<<"PHI_MIN        = "<<PHI_MIN<<endl;
	cout<<"PHI_MAX        = "<<PHI_MAX<<endl;
	cout<<"THETA_MIN      = "<<THETA_MIN<<endl;
	cout<<"THETA_MAX      = "<<THETA_MAX<<endl;
	cout<<"OUTPUT_FILENAME = \""<<OUTPUT_FILENAME<<"\""<<endl;
	cout<<"-------------------------------------------------"<<endl;
	
}

//----------------------------
// Usage
//----------------------------
void Usage(void)
{
	cout<<endl;
	cout<<"Usage:"<<endl;
	cout<<"      genphoton -M numEvents [options]"<<endl;
	cout<<endl;
	cout<<"  options:"<<endl;
	cout<<"    -h                print this help message"<<endl;
	cout<<"    -M numEvents      set the number of events to produced (required)"<<endl;
	cout<<"    -Pmin p           set the lower photon momentum in GeV/c"<<endl;
	cout<<"    -Pmax p           set the upper photon momentum in GeV/c"<<endl;
	cout<<"    -Phimin phi       set the lower photon phi angle in degrees"<<endl;
	cout<<"    -Phimax phi       set the upper photon phi angle in degrees"<<endl;
	cout<<"    -Thetamin theta   set the lower photon theta angle in degrees"<<endl;
	cout<<"    -Thetamax theta   set the upper photon theta angle in degrees"<<endl;
	cout<<"    -o filename       set the output filename"<<endl;
	cout<<endl;
	cout<<"This program is essentially a single photon particle gun."<<endl;
	cout<<"It can produce single photon events that can be converted using"<<endl;
	cout<<"genr8_2_hddm into a format suitable as input to hdgeant."<<endl;
	cout<<"\"Why go to this trouble\" you ask, \"when hdgeant already has a"<<endl;
	cout<<"built-in particle gun?\" The reason is because I want to overlay"<<endl;
	cout<<"EM background events on these and hdgeant will only do that if the"<<endl;
	cout<<"event is read from an external source."<<endl;

	cout<<endl;
	exit(0);
}


