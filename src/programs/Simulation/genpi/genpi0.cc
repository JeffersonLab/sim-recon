

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <vector>
using namespace std;

double PI_ZERO_MASS = 0.13497;
unsigned int MAX_EVENTS=10000;
int NUM_TO_GEN=1;
double E_BEAM_MIN=PI_ZERO_MASS;
double E_BEAM_MAX=1.0;
int RUN_NUMBER=100;
string OUTPUT_FILENAME="genpi0.ascii";

double FORCE_THETA = -1000.0;
double FORCE_PHI = -1000.0;
double THETA_PHOTON_MIN = 0.0; // minimum angle for photons in radians
double THETA_PHOTON_MAX = M_PI; // minimum angle for photons in radians

#define GAMMA_TYPE 1
#define PI_TYPE 7

class pi0{
	public:
		double px,py,pz,E; // pi0
		double px1, py1, pz1, E1; // decay photon 1
		double px2, py2, pz2, E2; // decay photon 2
};


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
	
		// Determine inicident photon energy. (This is just
		// evenly sampled from the specified range.)
		double Etot = (double)random()/(double)RAND_MAX*(E_BEAM_MAX-E_BEAM_MIN) + E_BEAM_MIN;
	
		// Generate pi0s
		vector<pi0> pi0s;
		for(int i=0; i<NUM_TO_GEN; i++){
		
			// Beam energy is semi-randomly distributed over all pions.
			// If this is the last pion it just gets the energy left
			// in Etot. Otherwise, it gets between 1/4 and 3/4 of the
			// remaining energy. (This needs something more sophisticated,
			// but for now it will do).
			double E_pi0 = 0.0;
			if(i==NUM_TO_GEN-1){
				E_pi0 = Etot;
			}else{
				E_pi0 = Etot*(0.25+0.5*(double)random()/(double)RAND_MAX);
			}
		
			// 4-vector of pi0. Pick a random direction.
			double phi_pi0 = 2.0*M_PI*((double)random()/(double)RAND_MAX);
			double theta_pi0 = M_PI*((double)random()/(double)RAND_MAX);
			double p_pi0 = sqrt(E_pi0*E_pi0 - PI_ZERO_MASS*PI_ZERO_MASS);
			
			// If user specified angles, use them instead
			if(FORCE_THETA>-1000.0)theta_pi0=FORCE_THETA;
			if(FORCE_PHI>-1000.0)phi_pi0=FORCE_PHI;
			pi0 p;
			p.E = E_pi0;
			p.px = p_pi0*sin(theta_pi0)*cos(phi_pi0);
			p.py = p_pi0*sin(theta_pi0)*sin(phi_pi0);
			p.pz = p_pi0*cos(theta_pi0);

			// 4-vectors of 2 decay photons in rest frame of pi0
			// where they are isotropic and back to back.
			double phi = 2.0*M_PI*((double)random()/(double)RAND_MAX);
			double theta = M_PI*((double)random()/(double)RAND_MAX);
			
			p.E1 = PI_ZERO_MASS/2.0;
			p.px1 = p.E1*sin(theta)*cos(phi);
			p.py1 = p.E1*sin(theta)*sin(phi);
			p.pz1 = p.E1*cos(theta);
			
			p.E2 = p.E1;
			p.px2 = -p.px1;
			p.py2 = -p.py1;
			p.pz2 = -p.pz1;
			
			// Boost photons into lab frame using 4-vector of pi0
			// (http://rd11.web.cern.ch/RD11/rkb/PH14pp/node105.html)
			p.E1 = (p.E1*p.E + p.px1*p.px + p.py1*p.py + p.pz1*p.pz)/PI_ZERO_MASS;
			p.px1 += p.px*(PI_ZERO_MASS/2.0 + p.E1)/(p.E + PI_ZERO_MASS);
			p.py1 += p.py*(PI_ZERO_MASS/2.0 + p.E1)/(p.E + PI_ZERO_MASS);
			p.pz1 += p.pz*(PI_ZERO_MASS/2.0 + p.E1)/(p.E + PI_ZERO_MASS);

			p.E2 = (p.E2*p.E + p.px2*p.px + p.py2*p.py + p.pz2*p.pz)/PI_ZERO_MASS;
			p.px2 += p.px*(PI_ZERO_MASS/2.0 + p.E2)/(p.E + PI_ZERO_MASS);
			p.py2 += p.py*(PI_ZERO_MASS/2.0 + p.E2)/(p.E + PI_ZERO_MASS);
			p.pz2 += p.pz*(PI_ZERO_MASS/2.0 + p.E2)/(p.E + PI_ZERO_MASS);
			
			// Check that photons are within limits set by user
			double gtheta1 = acos(p.pz1/p.E1);
			double gtheta2 = acos(p.pz2/p.E2);
			if(gtheta1<THETA_PHOTON_MIN || gtheta1>THETA_PHOTON_MAX
				|| gtheta2<THETA_PHOTON_MIN || gtheta2>THETA_PHOTON_MAX){
				i--;
				continue;
			}
			
			pi0s.push_back(p);
			Etot -= E_pi0; // subtract this from beam energy available
		}
		
		// Write event to file
		of<<RUN_NUMBER<<" "<<nevents<<" "<<pi0s.size()*2<<endl;
		for(unsigned int j=0; j<pi0s.size(); j++){
			pi0 &p = pi0s[j];
			unsigned int index = 2*j+1;
			
			of<<index<<" "<<GAMMA_TYPE<<" "<<0<<endl;
			of<<"   "<<0<<" "<<p.px1<<" "<<p.py1<<" "<<p.pz1<<" "<<p.E1<<endl;
			index++;
			of<<index<<" "<<GAMMA_TYPE<<" "<<0<<endl;
			of<<"   "<<0<<" "<<p.px2<<" "<<p.py2<<" "<<p.pz2<<" "<<p.E2<<endl;
			
#if 0			
			double m2 = pow(p.E1+p.E2, 2.0)
						- pow(p.px1+p.px2, 2.0)
						- pow(p.py1+p.py2, 2.0)
						- pow(p.pz1+p.pz2, 2.0);
			cout<<"m="<<sqrt(m2)<<endl;
#endif

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
		}else if(arg=="-N"){
			if(i==narg-1){cout<<"-N requires an argument!"<<endl; Usage();}
			NUM_TO_GEN = atoi(argv[++i]);
		}else if(arg=="-Emin"){
			if(i==narg-1){cout<<"-Emin requires an argument!"<<endl; Usage();}
			E_BEAM_MIN = atof(argv[++i]);
		}else if(arg=="-Emax"){
			if(i==narg-1){cout<<"-Emax requires an argument!"<<endl; Usage();}
			E_BEAM_MAX = atof(argv[++i]);
		}else if(arg=="-o"){
			if(i==narg-1){cout<<"-o requires an argument!"<<endl; Usage();}
			OUTPUT_FILENAME = argv[++i];
		}else if(arg=="-theta"){
			if(i==narg-1){cout<<"-theta requires an argument!"<<endl; Usage();}
			FORCE_THETA = atof(argv[++i]);
		}else if(arg=="-phi"){
			if(i==narg-1){cout<<"-phi requires an argument!"<<endl; Usage();}
			FORCE_PHI = atof(argv[++i]);
		}else if(arg=="-m"){
			if(i==narg-1){cout<<"-m requires an argument!"<<endl; Usage();}
			PI_ZERO_MASS = atof(argv[++i]);
		}else if(arg=="-gtheta_min"){
			if(i==narg-1){cout<<"-gtheta_min requires an argument!"<<endl; Usage();}
			THETA_PHOTON_MIN = atof(argv[++i])*M_PI/180.0;
		}else if(arg=="-gtheta_max"){
			if(i==narg-1){cout<<"-gtheta_max requires an argument!"<<endl; Usage();}
			THETA_PHOTON_MAX = atof(argv[++i])*M_PI/180.0;
		}
	}
	
	cout<<"---- genpi0 will use the following settings: ----"<<endl;
	cout<<"MAX_EVENTS      = "<<MAX_EVENTS<<endl;
	cout<<"NUM_TO_GEN      = "<<NUM_TO_GEN<<endl;
	cout<<"E_BEAM_MIN      = "<<E_BEAM_MIN<<endl;
	cout<<"E_BEAM_MAX      = "<<E_BEAM_MAX<<endl;
	cout<<"THETA           = ";
		if(FORCE_THETA>-1000.0)cout<<FORCE_THETA; else cout<<"random";
		cout<<endl;
	cout<<"PHI             = ";
		if(FORCE_PHI>-1000.0)cout<<FORCE_PHI; else cout<<"random";
		cout<<endl;
	cout<<"THETA_PHOTON_MIN= "<<THETA_PHOTON_MIN*180.0/M_PI<<" degrees"<<endl;
	cout<<"THETA_PHOTON_MAX= "<<THETA_PHOTON_MAX*180.0/M_PI<<" degrees"<<endl;
	cout<<"PI_ZERO_MASS    = "<<PI_ZERO_MASS<<endl;
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
	cout<<"      genpi0 -M numEvents [options]"<<endl;
	cout<<endl;
	cout<<"  options:"<<endl;
	cout<<"    -h                  print this help message"<<endl;
	cout<<"    -M     numEvents    set the number of events to produced (required)"<<endl;
	cout<<"    -N     numPi0s      set the number of pi0s to generate per event."<<endl;
	cout<<"    -Emin  E            set the lower beam energy limit in GeV"<<endl;
	cout<<"    -Emax  E            set the upper beam energy limit in GeV"<<endl;
	cout<<"    -theta angle        set theta angle in rad (default: random)"<<endl;
	cout<<"    -phi   angle        set phi angle in rad (default: random)"<<endl;
	cout<<"    -gtheta_min angle   set min. theta angle for photons in deg. (default: 0)"<<endl;
	cout<<"    -gtheta_max angle   set max. theta angle for photons in deg. (default: 180)"<<endl;
	cout<<"    -o     filename     set the output filename"<<endl;
	cout<<"    -m     mass         set the rest mass of the pi0(GeV) (e.g. make it an eta!)"<<endl;
	cout<<endl;
	cout<<"This program will produce events with one or more pi0s and write"<<endl;
	cout<<"out the resulting decay photons in an ASCII file of the same"<<endl;
	cout<<"format as produced by genr8. There is pretty much no real physics"<<endl;
	cout<<"here other than the pi0 decays isotropically in its rest frame."<<endl;
	cout<<"The total beam energy is evenly sampled from the given range and"<<endl;
	cout<<"The pi0s divide that up in a semi-random way so that each should"<<endl;
	cout<<"get a reasonable amount of kinetic energy. The total momentum is"<<endl;
	cout<<"NOT conserved however as each pion is given a random angle"<<endl;
	cout<<"distributed isotropically."<<endl;
	cout<<endl;
	cout<<"This is intended for testing the calorimeter (BCAL and FCAL)"<<endl;
	cout<<"reconstruction code."<<endl;

	cout<<endl;
	exit(0);
}


