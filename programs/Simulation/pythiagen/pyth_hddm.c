
#include <stdio.h>

#include "HDDM/hddm_s.h"

s_iostream_t* hddmOutputStream=NULL;

typedef struct{
	int mech; /* what do the values of this correspond to */
	int kfid;
	int parent;
	int firstdaughter;
	int lastdaughter;
	int geantid;
}klund_t;

typedef struct{
	float px;
	float py;
	float pz;
	float mass;
}plund_t;

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
void write_hddm_event_(int *iev, float *beammom, int *nlnd, klund_t *klund, plund_t *plund)
{
	/* Loop over events */
	int i;
	static int Nevents = 0;
	static int Nevents_written = 0;
	int runNumber=2;
	float vertex[3]={0.0, 0.0, 65.0}; /* hardwired for center of target for now */

	/* hdgeant apparently has code to sample from an extended target if */
	/* the input event has a vertex position of 0,0,0. Here we set the */
	/* target position to 0,0,0 to trigger that. The original code for */
	/* setting the z-position is commented out below. */
	vertex[0] = vertex[1] = vertex[2] = 0.0;
	/* z vertex uniform in target cell 15cm long*/
	/*double z=(drand48()-0.5)*30.;*/
	/*vertex[2] += (float) z;*/

	Nevents++;

	/* Start a new event */
	s_PhysicsEvents_t* pes;
	s_Reactions_t* rs;
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
	rs->in[0].vertices = vs = make_s_Vertices(1);
	vs->mult = 1;
	vs->in[0].origin = origin = make_s_Origin();
	vs->in[0].products = ps = make_s_Products(*nlnd);
	ps->mult = 0;
	
	origin->t = 0.0;
	origin->vx = vertex[0];
	origin->vy = vertex[1];
	origin->vz = vertex[2];
	
	for(i=0;i<*nlnd; i++){
		double E2;
		//if(klund[i].geantid==0)continue;
		
		E2 = pow(plund[i].px,2.0) + pow(plund[i].py,2.0) + pow(plund[i].pz,2.0) + pow(plund[i].mass,2.0);
		ps->in[ps->mult].type = klund[i].geantid;
		ps->in[ps->mult].pdgtype = klund[i].kfid;
		ps->in[ps->mult].id = i+1;
		ps->in[ps->mult].parentid = klund[i].parent;
		ps->in[ps->mult].mech = klund[i].mech;
		
		ps->in[ps->mult].momentum = make_s_Momentum();
		ps->in[ps->mult].momentum->px = plund[i].px;
		ps->in[ps->mult].momentum->py = plund[i].py;
		ps->in[ps->mult].momentum->pz = plund[i].pz;
		ps->in[ps->mult].momentum->E  = sqrt(E2);
		ps->mult++;
	}
	
	if(*nlnd>0){
		Nevents_written++;
		flush_s_HDDM(thisOutputEvent, hddmOutputStream);
		if(Nevents_written%1000 == 0)printf("Wrote event %d events (%d generated)\n", Nevents_written, Nevents);
	}	
}
