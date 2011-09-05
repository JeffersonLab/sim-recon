// $Id$
//
//    File: JEventProcessor_extract_ptype_hddm.cc
// Created: Mon Sep  5 12:29:45 EDT 2011
// Creator: davidl (on Linux ifarm1101 2.6.18-128.7.1.el5 x86_64)
//


#include <iostream>
using namespace std;

#include "JEventProcessor_extract_ptype_hddm.h"
using namespace jana;

#include <particleType.h>
#include <TRACKING/DMCThrown.h>

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_extract_ptype_hddm());
}
} // "C"


float vertex[4]={0.0, 0.0, 65.0, 65.0};


//------------------
// JEventProcessor_extract_ptype_hddm (Constructor)
//------------------
JEventProcessor_extract_ptype_hddm::JEventProcessor_extract_ptype_hddm()
{
	pthread_mutex_init(&mutex, NULL);
	
	Nevents =0;
	
}

//------------------
// ~JEventProcessor_extract_ptype_hddm (Destructor)
//------------------
JEventProcessor_extract_ptype_hddm::~JEventProcessor_extract_ptype_hddm()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_extract_ptype_hddm::init(void)
{
	// Lock mutex
	pthread_mutex_lock(&mutex);

	// Get type of particle to extract
	PTYPE = Neutron;
	gPARMS->SetDefaultParameter("PTYPE", PTYPE, "GEANT particle type to extract to separate HDDM file.");

	// Get output filename
	OUTFILENAME = string(ParticleType((Particle_t)PTYPE))+".hddm";
	gPARMS->SetDefaultParameter("OUTPUT_FILENAME", OUTFILENAME, "Filename of HDDM file to write particles to.");

	// Open output file
	hddmout = init_s_HDDM((char*)OUTFILENAME.c_str());
	if(!hddmout){
		cout<<" Error opening output file \""<<OUTFILENAME<<"\"!"<<endl;
		exit(-1);
	}
	cout<<" output file: "<<OUTFILENAME<<endl;

	// Get vertex info
	string vertex_str = "0 0 65 65";
	gPARMS->SetDefaultParameter("VERTEX", vertex_str, "Vertex to throw particles from (should be string of 4 numbers x y zmin zmax)");
	sscanf(vertex_str.c_str(), "%f %f %f %f", &vertex[0], &vertex[1], &vertex[2], &vertex[3]);
	if(vertex[2] > vertex[3]){
	 cerr<<"Invalid parameter: z_min > z_max"<< endl;
	 exit(-1);
	}
	
	// Print message to user
	jout<<endl;
	jout<<"----------------------------------------"<<endl;
	jout<<"extract_ptype_hddm  plugin:"<<endl;
	jout<<endl;
	jout<<"   particle type: "<<ParticleType((Particle_t)PTYPE)<<" (set via -PPTYPE=geantid)"<<endl;
	jout<<" output filename: "<<OUTFILENAME<<" (set via -POUTFILENAME=fname.hddm)"<<endl;
	jout<<"          vertex: x="<<vertex[0]<<" y="<<vertex[1]<<" zmin="<<vertex[2]<<" zmax="<<vertex[3]<<endl;
	jout<<"                  (set via -PVERTEX=\"X Y Zmin Zmax\")"<<endl;
	jout<<"----------------------------------------"<<endl;
	jout<<endl;

	// Unlock mutex
	pthread_mutex_unlock(&mutex);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_extract_ptype_hddm::brun(JEventLoop *loop, int runnumber)
{
	// This is called whenever the run number changes
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_extract_ptype_hddm::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DMCThrown*> mcthrowns;
	loop->Get(mcthrowns);

	pthread_mutex_lock(&mutex);

	for(unsigned int i=0; i<mcthrowns.size(); i++){
		const DMCThrown *thrown = mcthrowns[i];
		
		if(thrown->type != (int)PTYPE)continue;
		
		// Start a new event
		s_PhysicsEvents_t* pes;
		s_Reactions_t* rs;
		s_Vertices_t* vs;
		s_Origin_t* origin;
		s_Products_t* ps;
		s_HDDM_t *thisOutputEvent = make_s_HDDM();
      thisOutputEvent->physicsEvents = pes = make_s_PhysicsEvents(1);
		pes->mult = 1;
		pes->in[0].runNo = 1;
		pes->in[0].eventNo = ++Nevents;
		pes->in[0].reactions = rs = make_s_Reactions(1);
		rs->mult = 1;
		rs->in[0].vertices = vs = make_s_Vertices(1);
		vs->mult = 1;
		vs->in[0].origin = origin = make_s_Origin();
		vs->in[0].products = ps = make_s_Products(1);
		ps->mult = 0;
		
		origin->t = 0.0;
		origin->vx = vertex[0];
		origin->vy = vertex[1];

		if(vertex[2]<vertex[3]){
		  origin->vz = randm(vertex[2],vertex[3]);
		}  else {
		  origin->vz = vertex[2];
		}

		DVector3 mom = thrown->momentum();

		ps->in[ps->mult].type = (Particle_t)thrown->type;
		ps->in[ps->mult].pdgtype = thrown->pdgtype;
		ps->in[ps->mult].id = thrown->myid;
		ps->in[ps->mult].parentid = thrown->parentid;
		ps->in[ps->mult].mech = thrown->mech;
		ps->in[ps->mult].momentum = make_s_Momentum();
		ps->in[ps->mult].momentum->px = mom.X();
		ps->in[ps->mult].momentum->py = mom.Y();
		ps->in[ps->mult].momentum->pz = mom.Z();
		ps->in[ps->mult].momentum->E  = thrown->energy();
		ps->mult++;
		
		flush_s_HDDM(thisOutputEvent, hddmout);
	}


	pthread_mutex_unlock(&mutex);

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_extract_ptype_hddm::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_extract_ptype_hddm::fini(void)
{
	pthread_mutex_lock(&mutex);
	if(hddmout)close_s_HDDM(hddmout);
	hddmout = NULL;
	pthread_mutex_unlock(&mutex);

	return NOERROR;
}

