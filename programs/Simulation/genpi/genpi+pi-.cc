

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
double PROTON_MASS = 0.938;
double MESON_MASS = 0.770;
unsigned int MAX_EVENTS=10000;
int NUM_TO_GEN=2;
double E_BEAM_MIN=4.0*PI_CHARGED_MASS;
double E_BEAM_MAX=1.0;
int RUN_NUMBER=100;
string OUTPUT_FILENAME="genpiX.ascii";
double gINTEGRAL_FRACTION;

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
	
		// Determine incident photon energy. (This is just
		// evenly sampled from the specified range.)
		double Eg = (double)random()/(double)RAND_MAX*(E_BEAM_MAX-E_BEAM_MIN) + E_BEAM_MIN;


#if 0
		// tmin occurs when the target proton goes out perpendicular to
		// the beam. This corresponds to the smallest theta angle the
		// rho can go to. Calculate th rho angle and momentum here.
		double rho_mom = 2.0*pow(Eg, 2.0);
		rho_mom -= pow(PROTON_MASS,2.0) - pow(MESON_MASS, 2.0);
		rho_mom /= 2.0*Eg;
		rho_mom = sqrt(pow(rho_mom, 2.0) - pow(MESON_MASS, 2.0));
#endif

		// Calculate rho theta and phi angles
		double rho_theta = 0.0;
		double rho_phi = 2.0*M_PI*((double)random()/(double)RAND_MAX);
#if 1
		// Calculate rho momentum based on angle
		double M_N = 938.0;
		double L = Eg/(Eg + M_N);
		double K = (Eg*M_N + pow(MESON_MASS,2.0)/2.0)/(Eg + M_N);
		double A = pow(L*cos(rho_theta),2.0) - 1.0;
		double B = 2.0*K*L*cos(rho_theta);
		double C = pow(K,2.0) - pow(MESON_MASS,2.0);
		double rho_mom = ((-B) - sqrt(B*B - 4.0*A*C))/(2.0*A);
#endif
		TVector3 rho_p;
		rho_p.SetMagThetaPhi(rho_mom, rho_theta, rho_phi);
		
		// Find decay angle of pions in rest frame of rho
		double pi1_theta = SampleSin2Theta();
		double pi1_phi = 2.0*M_PI*((double)random()/(double)RAND_MAX);
		double pi2_theta = M_PI - pi1_theta;
		double pi2_phi = pi1_phi + M_PI;
		if(pi2_phi>2.0*M_PI)pi2_phi-=2.0*M_PI;
		double pi_E = MESON_MASS/2.0;
		double pi_mom = sqrt(pow(pi_E,2.0) - pow(PI_CHARGED_MASS,2.0));
		TLorentzVector pi1( pi_mom*sin(pi1_theta)*cos(pi1_phi)
								, pi_mom*sin(pi1_theta)*sin(pi1_phi)
								, pi_mom*cos(pi1_theta)
								, pi_E);
		TLorentzVector pi2( pi_mom*sin(pi2_theta)*cos(pi2_phi)
								, pi_mom*sin(pi2_theta)*sin(pi2_phi)
								, pi_mom*cos(pi2_theta)
								, pi_E);

		// Boost the pions into the lab frame
		TVector3 beta = (1.0/sqrt(rho_p.Mag2() + pow(MESON_MASS,2.0)))*rho_p;
		pi1.Boost(beta);
		pi2.Boost(beta);

		// Generate piXs
		vector<piX> piXs;
		
		piX p;
		p.E = pi1.E();
		p.px = pi1.Px();
		p.py = pi1.Py();
		p.pz = pi1.Pz();
		piXs.push_back(p);

		p.E = pi2.E();
		p.px = pi2.Px();
		p.py = pi2.Py();
		p.pz = pi2.Pz();
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
			type = type==PI_PLUS_TYPE ? PI_MINUS_TYPE:PI_PLUS_TYPE;
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
// SampleSin2Theta
//----------------------------
double SampleSin2Theta(void)
{
	// Randomly sample a value of theta between 0 and pi
	// such than the distributions goes like sin^2.
	
	// Get random number evenly sampled from 0 to 1
	gINTEGRAL_FRACTION = (double)random()/(double)RAND_MAX;
	
	// The integral fraction of sin^2 as a function
	// of theta over the range zero pi is given by
	//
	// f = (theta - sin(2*theta))/pi
	//
	// Unfortunately, this is a trancendental equation
	// so we must find the root numerically.
	// We use the numerical recipes routine rtsafe
	// and define the function and derivative in
	// funcd.
	return zbrent(func, 0.0, M_PI, 1.0E-5);
}

//----------------------------
// funcd
//----------------------------
void funcd(double theta, double *f, double *df)
{
	*f = gINTEGRAL_FRACTION - (theta - sin(2.0*theta))/M_PI;
	*df = -(1.0 - 2.0*cos(2.0*theta))/M_PI;
}

//----------------------------
// func
//----------------------------
double func(double theta)
{
	return gINTEGRAL_FRACTION - (theta - 0.5*sin(2.0*theta))/M_PI;
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


