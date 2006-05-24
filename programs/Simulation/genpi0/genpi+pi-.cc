

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <vector>
using namespace std;

double PI_CHARGED_MASS = 0.139568;
uint MAX_EVENTS=10000;
int NUM_TO_GEN=1;
double E_BEAM_MIN=PI_CHARGED_MASS;
double E_BEAM_MAX=1.0;
int RUN_NUMBER=100;
string OUTPUT_FILENAME="genpiX.ascii";

#define GAMMA_TYPE 1
#define PI_PLUS_TYPE 8
#define PI_MINUS_TYPE 9

class piX{
	public:
		double px,py,pz,E; // piX
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
	uint nevents;
	for(nevents=1; nevents<=MAX_EVENTS; nevents++){
	
		// Determine inicident photon energy. (This is just
		// evenly sampled from the specified range.)
		double Etot = (double)random()/(double)RAND_MAX*(E_BEAM_MAX-E_BEAM_MIN) + E_BEAM_MIN;
	
		// Generate piXs
		vector<piX> piXs;
		for(int i=0; i<NUM_TO_GEN; i++){
		
			// Beam energy is semi-randomly distributed over all pions.
			// If this is the last pion it just gets the energy left
			// in Etot. Otherwise, it gets between 1/4 and 3/4 of the
			// remaining energy. (This needs something more sophisticated,
			// but for now it will do).
			double E_piX = 0.0;
			if(i==NUM_TO_GEN-1){
				E_piX = Etot;
			}else{
				E_piX = Etot*(0.25+0.5*(double)random()/(double)RAND_MAX);
			}
			Etot -= E_piX; // subtract this from beam energy available
		
			// 4-vector of piX. Pick a random direction.
			double phi_piX = 2.0*M_PI*((double)random()/(double)RAND_MAX);
			double theta_piX = M_PI*((double)random()/(double)RAND_MAX);
			double p_piX = sqrt(E_piX*E_piX - PI_CHARGED_MASS*PI_CHARGED_MASS);
			piX p;
			p.E = E_piX;
			p.px = p_piX*cos(theta_piX)*cos(phi_piX);
			p.py = p_piX*cos(theta_piX)*sin(phi_piX);
			p.pz = p_piX*sin(theta_piX);
		
			piXs.push_back(p);
		}
		
		// Write event to file
		of<<RUN_NUMBER<<" "<<nevents<<" "<<piXs.size()*2<<endl;
		for(uint j=0; j<piXs.size(); j++){
			piX &p = piXs[j];
			uint index = 2*j+1;
			
			of<<index<<" "<<PI_PLUS_TYPE<<" "<<0<<endl;
			of<<"   "<<0<<" "<<p.px<<" "<<p.py<<" "<<p.pz<<" "<<p.E<<endl;
			index++;
			of<<index<<" "<<PI_MINUS_TYPE<<" "<<0<<endl;
			of<<"   "<<0<<" "<<p.px<<" "<<p.py<<" "<<p.pz<<" "<<p.E<<endl;
			
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
		}else if(arg=="-m"){
			if(i==narg-1){cout<<"-m requires an argument!"<<endl; Usage();}
			PI_CHARGED_MASS = atof(argv[++i]);
		}
	}
	
	cout<<"---- genpiX will use the following settings: ----"<<endl;
	cout<<"MAX_EVENTS      = "<<MAX_EVENTS<<endl;
	cout<<"NUM_TO_GEN      = "<<NUM_TO_GEN<<endl;
	cout<<"E_BEAM_MIN      = "<<E_BEAM_MIN<<endl;
	cout<<"E_BEAM_MAX      = "<<E_BEAM_MAX<<endl;
	cout<<"PI_CHARGED_MASS    = "<<PI_CHARGED_MASS<<endl;
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
	cout<<"    -h            print this help message"<<endl;
	cout<<"    -M numEvents  set the number of events to produced (required)"<<endl;
	cout<<"    -N numpiXs    set the number of piXs to generate per event."<<endl;
	cout<<"    -Emin E       set the lower beam energy limit in GeV"<<endl;
	cout<<"    -Emax E       set the upper beam energy limit in GeV"<<endl;
	cout<<"    -o filename   set the output filename"<<endl;
	cout<<"    -m mass       set the rest mass of the pi(GeV) (e.g. make it an eta!)"<<endl;
	cout<<endl;
	cout<<"This program will produce events with one or more piXs and write"<<endl;
	cout<<"out the resulting decay photons in an ASCII file of the same"<<endl;
	cout<<"format as produced by genr8. There is pretty much no real physics"<<endl;
	cout<<"here other than the piX decays isotropically in its rest frame."<<endl;
	cout<<"The total beam energy is evenly sampled from the given range and"<<endl;
	cout<<"The piXs divide that up in a semi-random way so that each should"<<endl;
	cout<<"get a reasonable amount of kinetic energy. The total momentum is"<<endl;
	cout<<"NOT conserved however as each pion is given a random angle"<<endl;
	cout<<"distributed isotropically."<<endl;
	cout<<endl;
	cout<<"This is intended for testing the calorimeter (BCAL and FCAL)"<<endl;
	cout<<"reconstruction code."<<endl;

	cout<<endl;
	exit(0);
}


