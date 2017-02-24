// $Id$
//
//    File: hdbeam_current.cc
// Created: Tue Feb 21 20:51:51 EST 2017
// Creator: davidl (on Linux gluon48.jlab.org 2.6.32-431.20.3.el6.x86_64 x86_64)
//

#include <iostream>
#include <iomanip>
#include <string>
using namespace std;

#include <stdlib.h>

#include <DANA/DApplication.h>
#include <JANA/JParameterManager.h>

#include <DAQ/DBeamCurrent_factory.h>


bool QUIET = false;
bool PRINT_FIDUCIAL = false;
bool PRINT_TOTAL = false;
bool PRINT_TRIP_TIME = false;
bool PRINT_NTRIPS = false;
double BEAM_TRIP_MIN_T = 3.0; // seconds
double BEAM_ON_MIN_nA = 10.0; // nA
double T_MIN = 0.0;
double T_MAX = 0.0;
int32_t RUN=0;
		

void Usage(void);
void ParseCommandLineArgs(int narg, char* argv[]);

//-----------
// main
//-----------
int main(int narg, char *argv[])
{

	// Parse command line arguments and then create a DApplication
	ParseCommandLineArgs(narg, argv);
	DApplication *dapp = new DApplication(0, NULL);
	//dapp->Init();
	
	// Create a JEventLoop so we'll have a DBeamCurrent_factory
	JEventLoop *loop = new JEventLoop(dapp);
	loop->GetJEvent().SetRunNumber(RUN);

	// Get factory
	DBeamCurrent_factory *fac = (DBeamCurrent_factory*)loop->GetFactory("DBeamCurrent");
	if(!fac){
		cerr << "Unable to find DBeamCurrent_factory!" << endl;
		return -1;
	}
	
	// Read in trip map from CCDB
	gPARMS->SetParameter("BEAM_ON_MIN_nA",  BEAM_ON_MIN_nA);
	gPARMS->SetParameter("BEAM_TRIP_MIN_T", BEAM_TRIP_MIN_T);
	fac->init();
	fac->brun(loop, RUN);
	
	// Get some values from factory
	double t_total = fac->IntegratedTime();
	if(T_MAX==0.0) T_MAX = t_total;
	double t_fiducial = fac->IntegratedFiducialTime(T_MIN, T_MAX);
	uint32_t Nboundaries = fac->boundaries.size();
	uint32_t Ntrips = fac->trip.size();
	
	// Get avg. beam current as well as fraction of time beam was on
	double last_t     = 0.0;
	double last_Ibeam = 0.0;
	double Ion_sum    = 0.0;
	double Ion_t      = 0.0;
	double Ioff_t     = 0.0;
	for(auto &b : fac->boundaries){
		double delta_t = b.t - last_t;
		if(last_Ibeam >= BEAM_ON_MIN_nA){
			Ion_t   += delta_t;
			Ion_sum += delta_t*last_Ibeam;
		}else{
			Ioff_t  += delta_t;
		}
		last_t     = b.t;
		last_Ibeam = b.Ibeam;
	}
	double Ibeam_avg = Ion_sum/Ion_t;

	cout << endl;
	cout << "--------------------------------------------" << endl;
	cout << "Electron Beam Trip Report" << endl;
	cout << endl;
	cout << "              Run: " << RUN << endl;
	cout << "total time of run: " << t_total << " sec" << endl;
	cout << "     beam ON time: " << Ion_t << " sec (" << 100.0*Ion_t/t_total << "%)" << endl;
	cout << "    fiducial time: " << t_fiducial << " sec (" << 100.0*t_fiducial/t_total << "%)" << endl;
	cout << "      trip buffer: " << BEAM_TRIP_MIN_T << " sec" <<endl;
	cout << "   beam ON thresh: " << BEAM_ON_MIN_nA << " nA" << endl;
	cout << "          t range: " << T_MIN << " - " << T_MAX << " sec (" << 100.0*(T_MAX-T_MIN)/t_total << "%)" << endl;
	cout << "      Nboundaries: " << Nboundaries << " (current changed by >3nA)" << endl;
	cout << "           Ntrips: " << Ntrips << endl;
	cout << "         trips/hr: "  << (float)Ntrips/(float)(T_MAX - T_MIN)*3600.0 << endl;
	cout << "    avg. Ibeam ON: " << Ibeam_avg << " nA" << endl;
	cout << "--------------------------------------------" << endl;
	cout << endl;

	return 0;
}

//-----------------------
// Usage
//-----------------------
void Usage(void)
{
	cout<<endl;
	cout<<"Usage:"<<endl;
	cout<<"   hdbeam_current [options] run"<<endl;
	cout<<endl;
	cout<<" options:"<<endl;
	cout<<"    -h, --help    Show this Usage statement"<<endl;
// 	cout<<"    -tfiducial    Report the fiducial time (sec)"<<endl;
// 	cout<<"    -ttotal       Report the total time (sec)"<<endl;
// 	cout<<"    -ttrips       Report the total tripped time (sec)"<<endl;
// 	cout<<"    -Ntrips       Report the total number of beam trips" << endl;
	cout<<"    -tbuffer   #  Set the buffer time to exclude near trips (sec)"<<endl;
	cout<<"    -tmin      #  Set the minimum of the time range (sec)"<<endl;
	cout<<"    -tmax      #  Set the maximum of the time range (sec)"<<endl;
	cout<<"    -Ion       #  Set the beam current threshold for \"ON\" (nA)"<<endl;
	cout<<endl;
	cout<<" "
			"This will query the CCDB for electron the beam trip map for the \n"
			"given run number. The CCDB is filled from the EPICS archive using \n"
			"the mySampler program to interpolate values from the IBCAD00CRCUR6 \n"
			"current monitor. Changes in the beam current greater than 3nA are \n"
			"written to the archive as \"boundaries\". The above options can be\n"
			"used to interpret these boundaries in terms of beam trips. A \n"
			"fiducial time range is made from the regions between the trips with \n"
			"excluding a buffer just before trip and just after a recovery. This \n"
			"is intended to cut out regions where the beam may be slightly astray. \n"
			"\n"
// 			"By default, all numbers are reported. If any of -tfiducial, -ttotal,\n"
// 			"or -ttrips or -Ntrips is given, then only that number is printed and\n"
// 			"nothing else. This is to make it easier to use this in scripts.\n"  
// 			"\n"
			"The -tmin and -tmax options allow one to specify a time range in which\n"
			"to integrate fiducial time. If no minimum or maximum time are specified\n"
			"then the full time range of the run is used.\n"
			"\n"
			"Note that this functionality also exists within sim-recon via the\n"
			"DBeamCurrent objects and their \"is_fudicial\" flag. The \n"
			"DBeamCurrent_factory class can also be used to integrate the fiducial\n"
			"time in a run or given time range via the IntegratedFiducialTime()\n"
			"method.\n"
			"\n"
			"This utility is meant to give quick, easy access to the same information\n"
			"without having to write a special sim-recon program.\n" << endl;

}

//-----------------------
// ParseCommandLineArgs
//-----------------------
void ParseCommandLineArgs(int narg, char* argv[])
{
	for(int i=1; i<narg; i++){
		string arg(argv[i]);
		string next(i<narg-1 ? argv[i+1]:"");
		float argf = atof(next.c_str());
		bool used_next = false; // keep track if "next" is used so we can have a single error check below

		if(arg=="-tfiducial"){PRINT_FIDUCIAL=true;  QUIET=true;}
		if(arg=="-ttotal"   ){PRINT_TOTAL=true;     QUIET=true;}
		if(arg=="-ttrips"   ){PRINT_TRIP_TIME=true; QUIET=true;;}
		if(arg=="-Ntrips"   ){PRINT_NTRIPS=true;    QUIET=true;}
		if(arg=="-tbuffer"  ){used_next=true; BEAM_TRIP_MIN_T = argf;}
		if(arg=="-tmin"     ){used_next=true; T_MIN = argf;}
		if(arg=="-tmax"     ){used_next=true; T_MAX = argf;}
		if(arg=="-Ion"      ){used_next=true; BEAM_ON_MIN_nA = argf;}
		
		if(arg=="-h" || arg=="--help"){Usage(); exit(0);}
		
		if(used_next){
			i++;
			continue;			
		}
		
		if(arg[0] != '-') RUN =atoi(argv[i]);
	}
	
	if(RUN == 0){
		Usage();
		exit(0);
	}
}

