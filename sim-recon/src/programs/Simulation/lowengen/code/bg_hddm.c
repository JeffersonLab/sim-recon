
#include <stdio.h>

#include "HDDM/hddm_s.h"

s_iostream_t* hddmOutputStream=NULL;

typedef struct{
	int geantid;
	int decayproducts;
}keve_t;

typedef struct{
	float px;
	float py;
	float pz;
	float en;
}peve_t;

/*-----------------
// open_hddm_output_
//-----------------*/
void open_hddm_output_(const char *outputfile, int len)
{
	/* Copy FORTRAN string into a C-style string */
	char outfile[256];
	strncpy(outfile, outputfile, len);
	outfile[len]=0;

	/* Open output file */
	hddmOutputStream = init_s_HDDM(outfile);
	if(!hddmOutputStream){
		fprintf(stderr, "Unable to open output file \"%s\" for writing.\n", outfile);
		exit(-3);
	}
	
	printf("Opened HDDM file \"%s\" for writing ...\n", outfile);
}

/*-----------------
// close_hddm_output_
//-----------------*/
void close_hddm_output_(void)
{
	/* Close output file */
	close_s_HDDM(hddmOutputStream);
	
	printf("Closed HDDM output file\n");
}

/*-----------------
// write_hddm_event_
//-----------------*/
void write_hddm_event_(int *iev, int *iproc, int *ntra, keve_t *keve, peve_t *peve)
{
	/* Loop over events */
	int i;
	static int Nevents = 0;
	static int Nevents_written = 0;
	int runNumber=2;
	float vertex[3]={0.0, 0.0, 65.0}; /* hardwired for center of target for now */

	/* z vertex uniform in target cell 15cm long*/
	double z=(drand48()-0.5)*30.;
	vertex[2] += (float) z;

	Nevents++;

	/* Start a new event */
	s_PhysicsEvents_t* pes;
	s_Reactions_t* rs;
	s_Beam_t* bs;
	s_Target_t* ts;
	s_Vertices_t* vs;
	s_Origin_t* origin;
	s_Products_t* ps;

	s_HDDM_t *thisOutputEvent = make_s_HDDM();
	thisOutputEvent->physicsEvents = pes = make_s_PhysicsEvents(1);
	pes->mult = 1;
	pes->in[0].runNo = runNumber;
	pes->in[0].eventNo = Nevents;
	pes->in[0].reactions = rs = make_s_Reactions(1);
	rs->mult = 1;
	rs->in[0].type = *iproc;

	rs->in[0].beam = bs = make_s_Beam();
        bs->type = keve[0].geantid;
        bs->momentum->px = peve[0].px;
        bs->momentum->py = peve[0].py;
        bs->momentum->pz = peve[0].pz;
        bs->momentum->E  = peve[0].en;
        
	rs->in[0].target = ts = make_s_Target();
        ts->type = keve[1].geantid;
        ts->momentum->px = peve[1].px;
        ts->momentum->py = peve[1].py;
        ts->momentum->pz = peve[1].pz;
        ts->momentum->E  = peve[1].en;
        
	rs->in[0].vertices = vs = make_s_Vertices(1);
	vs->mult = 1;
	vs->in[0].origin = origin = make_s_Origin();
	vs->in[0].products = ps = make_s_Products(*ntra);
	ps->mult = 0;
	
	origin->t = 0.0;
	origin->vx = vertex[0];
	origin->vy = vertex[1];
	origin->vz = vertex[2];
	
	for(i=2;i<*ntra; i++){
		double E2;
		//if(keve[i].geantid==0)continue;
		
		ps->in[ps->mult].type = keve[i].geantid;
		ps->in[ps->mult].pdgtype = 0;
		ps->in[ps->mult].id = i+1;
		ps->in[ps->mult].parentid = 0;
		ps->in[ps->mult].mech = 1;
		
		ps->in[ps->mult].momentum = make_s_Momentum();
		ps->in[ps->mult].momentum->px = peve[i].px;
		ps->in[ps->mult].momentum->py = peve[i].py;
		ps->in[ps->mult].momentum->pz = peve[i].pz;
		ps->in[ps->mult].momentum->E  = peve[i].en;
		ps->mult++;
	}
	
	if(*ntra>0){
		Nevents_written++;
		flush_s_HDDM(thisOutputEvent, hddmOutputStream);
		if(Nevents_written%1000 == 0)printf("Wrote event %d events (%d generated)\n", Nevents_written, Nevents);
	}	
}
