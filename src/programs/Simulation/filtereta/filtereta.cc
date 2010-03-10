// $Id: mcsmear.cc 2388 2007-01-10 16:46:03Z davidl $
//
// Created August 24, 2007  David Lawrence

#include <iostream>
#include <iomanip>
using namespace std;

#include <signal.h>
#include <time.h>

#include "HDDM/hddm_s.h"

bool Filter(s_HDDM_t *hddm_s);
bool Filter_eta_p(s_HDDM_t *hddm_s);
bool Filter_eta_p_pi0(s_HDDM_t *hddm_s);
bool Filter_eta_n_pip(s_HDDM_t *hddm_s);
bool Filter_eta_p_gamma(s_HDDM_t *hddm_s);
void ParseCommandLineArguments(int narg, char* argv[]);
void Usage(void);
void ctrlCHandle(int x);
s_iostream_t* NewHDDMFile(char *fname);
s_HDDM_t* CloneEvent(s_HDDM_t *hddm_s);


char *INFILENAME = NULL;
char *OUTFILENAME = NULL;
int QUIT = 0;
bool CREATE_SINGLE_CHANNEL_FILES = false;


#ifndef _DBG_
#define _DBG_ cerr<<__FILE__<<":"<<__LINE__<<" "
#define _DBG__ _DBG_<<endl
#endif

//-----------
// main
//-----------
int main(int narg,char* argv[])
{
	// Set up to catch SIGINTs for graceful exits
	signal(SIGINT,ctrlCHandle);

	ParseCommandLineArguments(narg, argv);
		
	// Open Input file
	cout<<" input file: "<<INFILENAME<<endl;
	s_iostream_t *fin = open_s_HDDM(INFILENAME);
	if(!fin){
		cout<<" Error opening input file \""<<INFILENAME<<"\"!"<<endl;
		exit(-1);
	}
	
	// Output file
	s_iostream_t *fout = NewHDDMFile(OUTFILENAME);
	
	// Optionally create outputfiles for single channels
	s_iostream_t *fout_eta_p = NULL;
	s_iostream_t *fout_eta_p_pi0 = NULL;
	s_iostream_t *fout_eta_n_pip = NULL;
	s_iostream_t *fout_eta_p_gamma = NULL;
	int NEvents_written_eta_p = 0;
	int NEvents_written_eta_p_pi0 = 0;
	int NEvents_written_eta_n_pip = 0;
	int NEvents_written_eta_p_gamma = 0;
	if(CREATE_SINGLE_CHANNEL_FILES){
		fout_eta_p = NewHDDMFile("eta_p.hddm");
		fout_eta_p_pi0 = NewHDDMFile("eta_p_pi0.hddm");
		fout_eta_n_pip = NewHDDMFile("eta_n_pip.hddm");
		fout_eta_p_gamma = NewHDDMFile("eta_p_gamma.hddm");
	}

	// Loop over events in input file
	s_HDDM_t *hddm_s;
	int NEvents_read = 0;
	int NEvents_written = 0;
	time_t last_time = time(NULL);
	while((hddm_s = read_s_HDDM(fin))){
		NEvents_read++;
		time_t now = time(NULL);
		if(now != last_time){
			cout<<" "<<NEvents_read<<" events read -- "<<NEvents_written<<" events written      \r";cout.flush();
			last_time = now;
		}
		
		// Write or don't depending on return value of Filter()
		if(Filter(hddm_s)){
			if(fout_eta_p && Filter_eta_p(hddm_s)){
				s_HDDM_t *hddm_s_copy = CloneEvent(hddm_s);
				flush_s_HDDM(hddm_s_copy, fout_eta_p);
				NEvents_written_eta_p++;
			}
			if(fout_eta_p_pi0 && Filter_eta_p_pi0(hddm_s)){
				s_HDDM_t *hddm_s_copy = CloneEvent(hddm_s);
				flush_s_HDDM(hddm_s_copy, fout_eta_p_pi0);
				NEvents_written_eta_p_pi0++;
			}
			if(fout_eta_n_pip && Filter_eta_n_pip(hddm_s)){
				s_HDDM_t *hddm_s_copy = CloneEvent(hddm_s);
				flush_s_HDDM(hddm_s_copy, fout_eta_n_pip);
				NEvents_written_eta_n_pip++;
			}
			if(fout_eta_p_gamma && Filter_eta_p_gamma(hddm_s)){
				s_HDDM_t *hddm_s_copy = CloneEvent(hddm_s);
				flush_s_HDDM(hddm_s_copy, fout_eta_p_gamma);
				NEvents_written_eta_p_gamma++;
			}

			flush_s_HDDM(hddm_s, fout);
			NEvents_written++;
			
		}
		
		if(QUIT)break;
	}
	
	// close input and output files
	close_s_HDDM(fin);
	close_s_HDDM(fout);

	cout<<endl<<"FINAL:"<<endl;
	cout<<" "<<NEvents_read<<" events read -- "<<NEvents_written<<" events written"<<endl;
	double frac = (double)NEvents_written/(double)NEvents_read;
	cout<<"Output eta file has "<<100.0*frac<<"%.  (about "<<frac*122.0<<" microbarns)"<<endl;

	if(fout_eta_p){
		close_s_HDDM(fout_eta_p);
		double frac = (double)NEvents_written_eta_p/(double)NEvents_read;
		cout<<"Output eta p file has       "<<100.0*frac<<"%.  (about "<<frac*122000.0<<" nb)"<<endl;
	}
	if(fout_eta_p_pi0){
		close_s_HDDM(fout_eta_p_pi0);
		double frac = (double)NEvents_written_eta_p_pi0/(double)NEvents_read;
		cout<<"Output eta p pi0 file has   "<<100.0*frac<<"%.  (about "<<frac*122000.0<<" nb)"<<endl;
	}
	if(fout_eta_n_pip){
		close_s_HDDM(fout_eta_n_pip);
		double frac = (double)NEvents_written_eta_n_pip/(double)NEvents_read;
		cout<<"Output eta n pi+ file has   "<<100.0*frac<<"%.  (about "<<frac*122000.0<<" nb)"<<endl;
	}
	if(fout_eta_p_gamma){
		close_s_HDDM(fout_eta_p_gamma);
		double frac = (double)NEvents_written_eta_p_gamma/(double)NEvents_read;
		cout<<"Output eta p gamma file has "<<100.0*frac<<"%.  (about "<<frac*122000.0<<" nb)"<<endl;
	}
	

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
				case 'c': CREATE_SINGLE_CHANNEL_FILES=true;					break;
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
	sprintf(str, "%s_filtered.hddm", path_stripped);
	OUTFILENAME = strdup(str);
}


//-----------
// Usage
//-----------
void Usage(void)
{
	cout<<endl<<"Usage:"<<endl;
	cout<<"     filtergen [options] file.hddm"<<endl;
	cout<<endl;
	cout<<" Read filter events from the given HDDM file while copying it"<<endl;
	cout<<"to another HDDM file. This was written to provide an easy way"<<endl;
	cout<<"to filter events from a file produced by pythia that one does"<<endl;
	cout<<"not want to waste time tracking through the whole detector."<<endl;
	cout<<endl;
	cout<<"At this point, the filtering conditions are hardwired and"<<endl;
	cout<<"can only be changed by editing the source. In the future, we"<<endl;
	cout<<"will want command line options so the filter can be tuned"<<endl;
	cout<<"without recompiling."<<endl;
	cout<<endl;
	cout<<"  options:"<<endl;
	cout<<"    -h       Print this usage statement."<<endl;
	cout<<"    -c       Create single channel HDDM files as well."<<endl;
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

//-----------------------------------------------------------------
// NewHDDMFile
//-----------------------------------------------------------------
s_iostream_t* NewHDDMFile(char *fname)
{
	// Open a new HDDM output file with the specified filename.
	// Print an error message and exit immediately if there is
	// a problem.

	s_iostream_t *fout = init_s_HDDM(fname);
	if(!fout){
		cout<<" Error opening output file \""<<fname<<"\"!"<<endl;
		exit(-1);
	}
	cout<<" output file: "<<fname<<endl;
	
	return fout;
}

//-----------------------------------------------------------------
// CloneEvent
//-----------------------------------------------------------------
s_HDDM_t* CloneEvent(s_HDDM_t *hddm_s)
{
	// Clone the given event. This will only clone the thrown
	// particles in the event, not any of hits tree. It is the
	// caller's responsibility to call flush_s_HDDM(...) with the
	// returned pointer in order to free the memory.
	//
	// This routine is needed in order to write a single event to
	// multiple outputs because HDDM automatically frees the memory
	// of an hddm_s pointer whenever it writes it to a file. Also,
	// the HDDM library provides not routines to copy events so we
	// must do it "by hand" here.
	// Loop over Physics Events
	
	s_HDDM_t *hddm_s_copy = make_s_HDDM();

	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE){
		hddm_s_copy->physicsEvents = NULL;
		return hddm_s_copy;
	}

	s_PhysicsEvents_t* PE_copy = hddm_s_copy->physicsEvents = make_s_PhysicsEvents(PE->mult);
	PE_copy->mult = PE->mult;
	
	for(unsigned int i=0; i<PE->mult; i++){
		PE_copy->in[i].eventNo = PE->in[i].eventNo;
		PE_copy->in[i].runNo = PE->in[i].runNo;
		PE_copy->in[i].hitView = (s_HitView_t*)HDDM_NULL; // never copy hits
	
		// ------------ Reactions --------------
		s_Reactions_t *reactions=PE->in[i].reactions;
		if(!reactions){
			PE_copy->in[i].reactions = NULL;
			continue;
		}

		s_Reactions_t *reactions_copy = PE_copy->in[i].reactions = make_s_Reactions(reactions->mult);
		reactions_copy->mult = reactions->mult;

		for(unsigned int j=0; j<reactions->mult; j++){
			reactions_copy->in[j].type = reactions->in[j].type;
			reactions_copy->in[j].weight = reactions->in[j].weight;

			// ------------ Beam --------------
			s_Beam_t *beam = reactions->in[j].beam;
			s_Beam_t *beam_copy = reactions_copy->in[j].beam =  beam ?  make_s_Beam():NULL;
			if(beam_copy){
				beam_copy->type = beam->type;
				s_Momentum_t *mom = beam_copy->momentum = make_s_Momentum();
				s_Properties_t *prop = beam_copy->properties = make_s_Properties();
				*mom = *beam->momentum;
				*prop = *beam->properties;
			}

			// ------------ Target --------------
			s_Target_t *target = reactions->in[j].target;
			s_Target_t *target_copy = reactions_copy->in[j].target =  target ?  make_s_Target():NULL;
			if(target_copy){
				target_copy->type = target->type;
				s_Momentum_t *mom = target_copy->momentum = make_s_Momentum();
				s_Properties_t *prop = target_copy->properties = make_s_Properties();
				*mom = *target->momentum;
				*prop = *target->properties;
			}
			
			// ------------ Vertices --------------
			s_Vertices_t *vertices = reactions->in[j].vertices;
			if(vertices){
				s_Vertices_t* vertices_copy = reactions_copy->in[j].vertices = make_s_Vertices(vertices->mult);
				vertices_copy->mult = vertices->mult;
				for(unsigned int k=0; k<vertices->mult; k++){
					s_Origin_t *origin = vertices->in[k].origin;
					s_Products_t *products = vertices->in[k].products;
					
					s_Origin_t *origin_copy = vertices_copy->in[k].origin =  origin ?  make_s_Origin():NULL;
					s_Products_t *products_copy = vertices_copy->in[k].products =  products ?  make_s_Products(products->mult):NULL;
					
					if(origin && origin_copy){
						origin_copy->t = origin->t;
						origin_copy->vx = origin->vx;
						origin_copy->vy = origin->vy;
						origin_copy->vz = origin->vz;
					}
					
					if(products && products_copy){
						products_copy->mult = products->mult;
						for(unsigned int m=0;m<products->mult;m++){
							s_Product_t *product = &products->in[m];
							s_Product_t *product_copy = &products_copy->in[m];

							product_copy->decayVertex = product->decayVertex;
							product_copy->id = product->id;
							product_copy->mech = product->mech;
							product_copy->parentid = product->parentid;
							product_copy->pdgtype = product->pdgtype;
							product_copy->type = product->type;
							product_copy->momentum = (s_Momentum_t*)HDDM_NULL;
							product_copy->properties = (s_Properties_t*)HDDM_NULL;
							
							if((product->momentum!=NULL) && (product->momentum!=HDDM_NULL)){
								s_Momentum_t *mom = product_copy->momentum = make_s_Momentum();
								mom->E = product->momentum->E;
								mom->px = product->momentum->px;
								mom->py = product->momentum->py;
								mom->pz = product->momentum->pz;
							}

							if((product->properties!=NULL) && (product->properties!=HDDM_NULL)){
								s_Properties_t *prop = product_copy->properties = make_s_Properties();
								prop->charge = product->properties->charge;
								prop->mass = product->properties->mass;
							}
						}
					}
				}
			}else{
				reactions_copy->in[j].vertices = NULL;
			}
		}
	}

	return hddm_s_copy;
}

