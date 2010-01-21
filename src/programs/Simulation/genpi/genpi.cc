

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

double PI_CHARGED_MASS = 0.139568;
unsigned int MAX_EVENTS=10000;
double P_MIN=0.100;
double P_MAX=6.000;
double PHI_MIN = 0.0;
double PHI_MAX = 2.0*M_PI;
double THETA_MIN = 0.0;
double THETA_MAX = M_PI;
bool IS_POSITIVE = true;

int RUN_NUMBER=100;
string OUTPUT_FILENAME="genpi.ascii";

#define GAMMA_TYPE 1
#define PI_PLUS_TYPE 8
#define PI_MINUS_TYPE 9

#define _DBG_ cout<<__FILE__<<":"<<__LINE__<<" "
#define _DBG__ _DBG_<<endl

class piX{
	public:
		double px,py,pz,E; // piX
};


extern "C"{
double rtsafe(void (*funcd)(double, double *, double *), double x1, double x2, double xacc);
double rtnewt(void (*funcd)(double, double *, double *), double x1, double x2, double xacc);
double zbrent(double (*func)(double), double x1, double x2, double tol);
void funcd(double, double *f, double *df);
double func(double);
}

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
	
		vector<piX> piXs;
		piX p;
	
		// Randomly sample the energy and angles of the pion
		double mom = (double)random()/(double)RAND_MAX*(P_MAX-P_MIN) + P_MIN;
		double phi = (double)random()/(double)RAND_MAX*(PHI_MAX-PHI_MIN) + PHI_MIN;
		double theta = (double)random()/(double)RAND_MAX*(THETA_MAX-THETA_MIN) + THETA_MIN;

		p.E = sqrt(mom*mom + PI_CHARGED_MASS*PI_CHARGED_MASS);	
		p.px = mom*sin(theta)*cos(phi);
		p.py = mom*sin(theta)*sin(phi);
		p.pz = mom*cos(theta);
		piXs.push_back(p);
			
		// Write event to file
		unsigned int type = PI_PLUS_TYPE;
		of<<RUN_NUMBER<<" "<<nevents<<" "<<piXs.size()<<endl;
		for(unsigned int j=0; j<piXs.size(); j++){
			piX &p = piXs[j];
			unsigned int index = j+1;
			
			of<<index<<" "<<type<<" "<<0<<endl;
			of<<"   "<<0<<" "<<p.px<<" "<<p.py<<" "<<p.pz<<" "<<p.E<<endl;
			
			// alternate bewtween pi+ and pi-
			type = IS_POSITIVE ? PI_PLUS_TYPE:PI_MINUS_TYPE;
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
	
	cout<<"---- genpiX will use the following settings: ----"<<endl;
	cout<<"MAX_EVENTS     = "<<MAX_EVENTS<<endl;
	cout<<"P_MIN          = "<<P_MIN<<endl;
	cout<<"P_MAX          = "<<P_MAX<<endl;
	cout<<"PHI_MIN        = "<<PHI_MIN<<endl;
	cout<<"PHI_MAX        = "<<PHI_MAX<<endl;
	cout<<"THETA_MIN      = "<<THETA_MIN<<endl;
	cout<<"THETA_MAX      = "<<THETA_MAX<<endl;
	cout<<"PI_CHARGE      = "<<(IS_POSITIVE ? "+1":"-1")<<endl;
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
	cout<<"      genpiX -M numEvents [options]"<<endl;
	cout<<endl;
	cout<<"  options:"<<endl;
	cout<<"    -h                print this help message"<<endl;
	cout<<"    -M numEvents      set the number of events to produced (required)"<<endl;
	cout<<"    -N numpiXs        set the number of piXs to generate per event."<<endl;
	cout<<"    -Pmin p           set the lower pion momentum in GeV/c"<<endl;
	cout<<"    -Pmax p           set the upper pion momentum in GeV/c"<<endl;
	cout<<"    -Phimin phi       set the lower pion phi angle in degrees"<<endl;
	cout<<"    -Phimax phi       set the upper pion phi angle in degrees"<<endl;
	cout<<"    -Thetamin theta   set the lower pion theta angle in degrees"<<endl;
	cout<<"    -Thetamax theta   set the upper pion theta angle in degrees"<<endl;
	cout<<"    -n                set the particle type to a pi-"<<endl;
	cout<<"    -o filename       set the output filename"<<endl;
	cout<<endl;
	cout<<"This program is essentially a single pion particle gun."<<endl;
	cout<<"It can produce single pion events that can be converted using"<<endl;
	cout<<"genr8_2_hddm into a format suitable as input to hdgeant."<<endl;
	cout<<"\"Why go to this trouble\" you ask, \"when hdgeant already has a"<<endl;
	cout<<"built-in particle gun?\" The reason is because I want to overlay"<<endl;
	cout<<"EM background events on these and hdgeant will only do that if the"<<endl;
	cout<<"event is read from an external source."<<endl;

	cout<<endl;
	exit(0);
}


