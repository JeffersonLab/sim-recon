

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <vector>
using namespace std;

double MUON_CHARGED_MASS = 0.10566;
unsigned int MAX_EVENTS=10000;
int NUM_TO_GEN=2;
double E_BEAM_MIN=4.0*MUON_CHARGED_MASS;
double E_BEAM_MAX=1.0;
int RUN_NUMBER=100;
string OUTPUT_FILENAME="genmuX.ascii";

#define GAMMA_TYPE 1
#define MUON_PLUS_TYPE 5
#define MUON_MINUS_TYPE 6

class muX{
	public:
		double px,py,pz,E; // muX
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

	double pmax = sqrt(E_BEAM_MAX*E_BEAM_MAX - MUON_CHARGED_MASS*MUON_CHARGED_MASS);
	double pmin = sqrt(E_BEAM_MIN*E_BEAM_MIN - MUON_CHARGED_MASS*MUON_CHARGED_MASS);

	// Loop over events
	unsigned int nevents;
	for(nevents=1; nevents<=MAX_EVENTS; nevents++){
	
		// Determine inicident photon energy. (This is just
		// evenly sampled from the specified range.)
		//double Etot = (double)random()/(double)RAND_MAX*(E_BEAM_MAX-E_BEAM_MIN) + E_BEAM_MIN;
	
		// Generate muXs
		vector<muX> muXs;
		for(int i=0; i<NUM_TO_GEN; i++){

#if 0		
			// Beam energy is semi-randomly distributed over all pions.
			// If this is the last pion it just gets the energy left
			// in Etot. Otherwise, it gets between 1/4 and 3/4 of the
			// remaining energy. (This needs something more sophisticated,
			// but for now it will do).
			double E_muX = 0.0;
			if(i==NUM_TO_GEN-1){
				E_muX = Etot;
			}else{
				E_muX = Etot*(0.25+0.5*(double)random()/(double)RAND_MAX);
			}
			Etot -= E_muX; // subtract this from beam energy available
			double p_muX = sqrt(E_muX*E_muX - MUON_CHARGED_MASS*MUON_CHARGED_MASS);
#else
			// This is changed from the above. For most studies, we actually
			// want an iso-momentum distribution of pions.
			double p_muX = (double)random()/(double)RAND_MAX*(pmax-pmin) + pmin;
			double E_muX = sqrt(MUON_CHARGED_MASS*MUON_CHARGED_MASS + p_muX*p_muX);
#endif		

			// 4-vector of first pion. Pick a random direction.
			double phi_muX = 2.0*M_PI*((double)random()/(double)RAND_MAX);
			double theta_muX = M_PI*((double)random()/(double)RAND_MAX);
			
			muX p;
			p.E = E_muX;
			p.px = p_muX*sin(theta_muX)*cos(phi_muX);
			p.py = p_muX*sin(theta_muX)*sin(phi_muX);
			p.pz = p_muX*cos(theta_muX);
			
			// Move angles to keep tracks apart
			//phi_muX += 2.0*M_PI/(double)NUM_TO_GEN;
			//if(phi_muX>2.0*M_PI)phi_muX-=2.0*M_PI;
			//theta_muX += 2.0*(M_MUON_2 - theta_muX);
		
			muXs.push_back(p);
		}
		
		// Write event to file
		unsigned int type = MUON_PLUS_TYPE;
		of<<RUN_NUMBER<<" "<<nevents<<" "<<muXs.size()<<endl;
		for(unsigned int j=0; j<muXs.size(); j++){
			muX &p = muXs[j];
			unsigned int index = j+1;
			
			of<<index<<" "<<type<<" "<<0<<endl;
			of<<"   "<<0<<" "<<p.px<<" "<<p.py<<" "<<p.pz<<" "<<p.E<<endl;
			
			// alternate bewtween pi+ and pi-
			type = type==MUON_PLUS_TYPE ? MUON_MINUS_TYPE:MUON_PLUS_TYPE;
			
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
			MUON_CHARGED_MASS = atof(argv[++i]);
		}
	}
	
	cout<<"---- genmuX will use the following settings: ----"<<endl;
	cout<<"MAX_EVENTS      = "<<MAX_EVENTS<<endl;
	cout<<"NUM_TO_GEN      = "<<NUM_TO_GEN<<endl;
	cout<<"E_BEAM_MIN      = "<<E_BEAM_MIN<<endl;
	cout<<"E_BEAM_MAX      = "<<E_BEAM_MAX<<endl;
	cout<<"MUON_CHARGED_MASS    = "<<MUON_CHARGED_MASS<<endl;
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
	cout<<"      genmuX -M numEvents [options]"<<endl;
	cout<<endl;
	cout<<"  options:"<<endl;
	cout<<"    -h            print this help message"<<endl;
	cout<<"    -M numEvents  set the number of events to produced (required)"<<endl;
	cout<<"    -N nummuXs    set the number of muXs to generate per event."<<endl;
	cout<<"    -Emin E       set the lower beam energy limit in GeV"<<endl;
	cout<<"    -Emax E       set the upper beam energy limit in GeV"<<endl;
	cout<<"    -o filename   set the output filename"<<endl;
	cout<<"    -m mass       set the rest mass of the pi(GeV) (e.g. make it an eta!)"<<endl;
	cout<<endl;
	cout<<"This program will produce events with one or more muXs and write"<<endl;
	cout<<"out the resulting decay photons in an ASCII file of the same"<<endl;
	cout<<"format as produced by genr8. There is pretty much no real physics"<<endl;
	cout<<"here other than the muX decays isotropically in its rest frame."<<endl;
	cout<<"The total beam energy is evenly sampled from the given range and"<<endl;
	cout<<"The muXs divide that up in a semi-random way so that each should"<<endl;
	cout<<"get a reasonable amount of kinetic energy. The total momentum is"<<endl;
	cout<<"NOT conserved however as each pion is given a random angle"<<endl;
	cout<<"distributed isotropically."<<endl;
	cout<<endl;
	cout<<"This is intended for testing the calorimeter (BCAL and FCAL)"<<endl;
	cout<<"reconstruction code."<<endl;

	cout<<endl;
	exit(0);
}


