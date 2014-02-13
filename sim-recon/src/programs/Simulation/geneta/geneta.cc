

#include "bg_hddm.h"
#include "c_cern.h"

#include <iostream>
using namespace std;

unsigned int Nevents_to_gen = 10368; // this must be changed in eta_p_gen.dat as well
bool use_interactive = false;
string output_fname="eta_primakoff.hddm";

void Usage(void);
void ParseCommandLineArgs(int narg, char *argv[]);


//------------------------
// main
//------------------------
int main(int narg, char *argv[])
{

	// Parse the command line arguments
	ParseCommandLineArgs(narg, argv);

	// Initialize cernlib stuff
   int nwgean=NWGEAN;
   int nwpaw=-NWPAW;
	quest_[9] = 65000;
	if(use_interactive){
		gpaw_(&nwgean, &nwpaw);
	}else{
		cout<<"Initializing. Please wait a moment ...."<<endl;
		gzebra_(&nwgean);
		hlimit_(&nwpaw);
		uginit_();
	}
	
	// Open output file
	open_hddm_output(output_fname);
	
	// Loop over events
	cout<<endl;
	cout<<endl;
	cout<<"Starting event generation ........."<<endl;
	cout<<"______________________________________________________"<<endl;
	int Nevents = 0;
	for(unsigned int i=0; i<Nevents_to_gen; i++){
		// Generate event
		gukine_();
		
		// Copy generated event parameters to Event object
		Event event;
		event.runNo = 1;
		event.eventNo = Nevents;
		event.reaction_type = 0;
		event.vertex.SetXYZ(0.0, 0.0, 0.0); // make hdgeant distribute vertex
		event.beam.type = Gamma;
		event.beam.p.SetXYZT(0.0, 0.0, kinem1_.EINI, kinem1_.EINI);
		event.target.type = Proton;
		event.target.p.SetXYZT(0.0, 0.0, 0.0, ParticleMass(Proton));
		
		Particle eta;
		float *p = kinem3_.PPI0LF;
		eta.type = Eta;
		eta.p.SetXYZT(p[0], p[1], p[2], kinem3_.EPI0LF);
		event.intermediate.push_back(eta);
		
		Particle gamma1;
		p = kinem3_.PG1LF;
		gamma1.type = Gamma;
		gamma1.parentid = 1; // parent is eta which is always particle 1
		gamma1.p.SetXYZT(p[0], p[1], p[2], kinem3_.EG1LF);
		event.final.push_back(gamma1);

		Particle gamma2;
		p = kinem3_.PG2LF;
		gamma2.type = Gamma;
		gamma2.parentid = 1; // parent is eta which is always particle 1
		gamma2.p.SetXYZT(p[0], p[1], p[2], kinem3_.EG2LF);
		event.final.push_back(gamma2);
		
		Particle proton;
		proton.type = Proton;
		proton.parentid = 0;
		proton.p = event.beam.p + event.target.p - eta.p; // elastic production is great!
		event.final.push_back(proton);
		
		// Write event to file
		write_hddm_event(event);

		Nevents++;
		if(Nevents%100==0){
			cout<<"  "<<Nevents<<" generated       \r";
			cout.flush();
		}
	}

	cout<<"______________________________________________________"<<endl;
	cout<<"Generated "<<Nevents<<" events total"<<endl;

	close_hddm_output();

	return 0;
}

//------------------------
// Usage
//------------------------
void Usage(void)
{
	cout<<endl;
	cout<<"Usage:"<<endl;
	cout<<"     geneta [options]"<<endl;
	cout<<endl;
	cout<<" Real Primakoff eta generator using a proton target."<<endl;
	cout<<"This is adapted from code originally authored by"<<endl;
	cout<<"Ashot Gasparian (gasparan@jlab.org)."<<endl;
	cout<<endl;
	cout<<" Output of this generator is HDDM format suitable"<<endl;
	cout<<"for use with hdgeant."<<endl;
	cout<<endl;
	cout<<" options:"<<endl;
	cout<<endl;
	cout<<" -h           print this usage statement"<<endl;
	cout<<" -o fname     set output filename to fname"<<endl;
	cout<<" -N Nevents   generate Nevents events"<<endl;
	cout<<" -i           run in interactive mode (default is batch)"<<endl;
	cout<<endl;
}

//------------------------
// ParseCommandLineArgs
//------------------------
void ParseCommandLineArgs(int narg, char *argv[])
{
	for(int i=1; i<narg; i++){
		string arg = argv[i];
		string next = (i+1)<narg ? argv[i+1]:"";
		if(arg=="-h"){
			Usage();
			exit(0);
		}
		
		if(arg=="-o"){
			if(next==""){
				Usage();
				cerr<<"ERROR: You must supply a filename when using \"-f\"!"<<endl;
				exit(-1);
			}
			output_fname = next;
			i++;
		}
		
		if(arg=="-N"){
			if(next==""){
				Usage();
				cerr<<"ERROR: You must a number of events to generate when using \"-N\"!"<<endl;
				exit(-2);
			}
			Nevents_to_gen = atoi(next.c_str());
			i++;
		}
		
		if(arg=="-i")use_interactive=true;
	}
}

