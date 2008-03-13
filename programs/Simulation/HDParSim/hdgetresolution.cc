
// This is s simple example program that shows how to use a DTrackResolution
// class to apply detector derived resolutions to charged tracks.
// There are detailed comments that should step you through the program
// with enough information to understand it and allow you to modify it to
// your own needs.


#include <iostream>
using namespace std;

#include <particleType.h>

#include "DTrackingResolutionGEANT.h"

void Usage(void);
void ParseCommandLineArguments(int narg, char *argv[]);

double PTOT=-1.0;
double THETA=-1.0;
double PTOTDELTA=0.0;
double THETADELTA=0.0;
int NPVALS=1;
int NTHETAVALS=1;
bool OUTPUT_C_PLUS_PLUS=false;

//------------
// main
//------------
int main(int narg, char *argv[])
{
	ParseCommandLineArguments(narg, argv);

	// Create a DTrackingResolution object. Specifically,
	// we create a DTrackingResolutionGEANT object so as
	// to use the resolutions derived from the GEANT-based
	// simulations.
	DTrackingResolution *res = new DTrackingResolutionGEANT();

	if(OUTPUT_C_PLUS_PLUS)cout<<"/*"<<endl;
	cout<<"# ptot is in GeV/c and theta is in degrees"<<endl;
	cout<<"# Momentum resolution is fractional"<<endl;
	cout<<"# Angular resolutions in milliradians"<<endl; 
	cout<<"# ptot\ttheta\tdpt/pt  \tdtheta\tdphi"<<endl;
	if(OUTPUT_C_PLUS_PLUS){
		cout<<"*/"<<endl;
		cout<<endl;
		cout<<"int NPVALS = "<<NPVALS<<";"<<endl;
		cout<<"int NTHETAVALS = "<<NTHETAVALS<<";"<<endl;
		cout<<"int Nresvals = "<<NPVALS*NTHETAVALS<<"; // (NPVALS*NTHETAVALS)"<<endl;
		cout<<"double PTOT = "<<PTOT<<";"<<endl;
		cout<<"double THETA = "<<THETA<<";"<<endl;
		cout<<"double PTOTDELTA = "<<PTOTDELTA<<";"<<endl;
		cout<<"double THETADELTA = "<<THETADELTA<<";"<<endl;
		cout<<endl;
		cout<<"class resvals_t{"<<endl;
		cout<<"	public:"<<endl;
		cout<<"	double p, theta, sigma_pt_over_pt, sigma_theta, sigma_phi;"<<endl;
		cout<<"};"<<endl;
		cout<<endl;
		cout<<"resvals_t resvals["<<NPVALS*NTHETAVALS<<"] = {"<<endl;
	}

	// Loop over p vals
	for(int i=0; i<NPVALS; i++){
		double p = PTOT;
		if(i>0)p += (double)i*PTOTDELTA/((double)NPVALS-1.0);
		
		// Loop over theta vals
		for(int j=0; j<NTHETAVALS; j++){
			double theta = THETA;
			if(j>0)theta += (double)j*THETADELTA/((double)NTHETAVALS-1.0);
		
			TVector3 mom;
			mom.SetMagThetaPhi(p, theta*M_PI/180.0, 0.0);
			double pt_res, theta_res, phi_res;
			res->GetResolution(8, mom, pt_res, theta_res, phi_res);
			
			if(OUTPUT_C_PLUS_PLUS){
				cout<<"{"<<p<<", "<<theta<<", "<<pt_res<<", "<<theta_res<<", "<<phi_res<<"}";
				if(((j+1)*(i+1))!=(NTHETAVALS*NPVALS))cout<<",";
				cout<<endl;
			}else{
				cout<<p<<"\t"<<theta<<"\t"<<pt_res<<"\t"<<theta_res<<"\t"<<phi_res<<endl;
			}
		}
	}
	
	if(OUTPUT_C_PLUS_PLUS){
		cout<<"};"<<endl<<endl;
	}

	return 0;
}

//------------
// ParseCommandLineArguments
//------------
void ParseCommandLineArguments(int narg, char *argv[])
{
	if(narg<3)Usage();

	for(int i=1; i<narg; i++){
		string arg(argv[i]);

		if(arg=="-ptotdelta"){
			if(i<narg-1){
				PTOTDELTA = atof(argv[++i]);
			}else{
				cerr<<"-ptotdelta requires an argument!"<<endl;
				Usage();
			}
		}else if(arg=="-thetadelta"){
			if(i<narg-1){
				THETADELTA = atof(argv[++i]);
			}else{
				cerr<<"-thetadelta requires an argument!"<<endl;
				Usage();
			}
		}else if(arg=="-npvals"){
			if(i<narg-1){
				NPVALS = atoi(argv[++i]);
			}else{
				cerr<<"-npvals requires an argument!"<<endl;
				Usage();
			}
		}else if(arg=="-nthetavals"){
			if(i<narg-1){
				NTHETAVALS = atoi(argv[++i]);
			}else{
				cerr<<"-nthetavals requires an argument!"<<endl;
				Usage();
			}
		}else if(arg=="-cpp"){
			OUTPUT_C_PLUS_PLUS=true;
		}else if(arg[0]!='-'){
			if(PTOT<0.0){
				PTOT = atof(arg.c_str());
			}else if(THETA<0.0){
				THETA = atof(arg.c_str());
			}else{
				cerr<<"Too many (non-switch) arguments!";
				Usage();
			} //if PTOT
		}else{
			cerr<<endl<<"Unkown argument \""<<arg<<"\"!"<<endl;
			Usage();
		} // if arg
	} // for
	
	if(THETA<0.0){
		cerr<<"No value given for theta!"<<endl;
		Usage();
	}

	if(PTOT<0.0){
		cerr<<"No value given for ptot!"<<endl;
		Usage();
	}
}

//------------
// Usage
//------------
void Usage(void)
{
	cout<<endl;
	cout<<"Usage:"<<endl;
	cout<<"  hdgetresolution [options] ptot theta"<<endl;
	cout<<endl;
	cout<<" options:"<<endl;
	cout<<"    -h             Print this help message"<<endl;
	cout<<"    -ptotdelta X   ptotdelta (in GeV/c)"<<endl;
	cout<<"    -thetadelta X  thetadelta (degrees)"<<endl;
	cout<<"    -npvals X      Number of moemtum values to print"<<endl;
	cout<<"    -nthetavals X  Number of theta values to print"<<endl;
	cout<<"    -cpp           Print output in C++ compilable format"<<endl; 
	cout<<endl;
	cout<<"   Print the charged particle resolutions to the sceen for a"<<endl;
	cout<<"specified total momentum and polar angle. If only the ptot"<<endl;
	cout<<"and theta arguments are given, then the resolutions for those"<<endl;
	cout<<"are printed. The -ptotdelta and -thetadelta switches can be"<<endl;
	cout<<"used to specify a range of values where the range starts at"<<endl;
	cout<<"ptot (theta) and ends at ptot+ptotdelta (theta+thetadelta)."<<endl;
	cout<<"The -npvals and -nthetavals switches allow one to set the"<<endl;
	cout<<"total number of values printed to the screen. The values will"<<endl;
	cout<<"be evenly spaced such that both end points are printed."<<endl;
	cout<<"Therefore, if the -ptotdelta argument is given, one also need"<<endl;
	cout<<"to supply a -npvals argument that is at least 2 or more."<<endl;
	cout<<"Otherwise, only the values at ptot will be printed. Similarly"<<endl;
	cout<<"for theta."<<endl;
	cout<<endl;
	cout<<"Note that for this to work with the GEANT derived resolutions"<<endl;
	cout<<"one needs to have the \"hd_res_charged.root\" file in the working"<<endl;
	cout<<"directory when this is run."<<endl;
	cout<<endl;
	
	exit(0);
}

