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
void InitPlugin(JApplication *app) {
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
   gPARMS->SetDefaultParameter("PTYPE", PTYPE, 
           "GEANT particle type to extract to separate HDDM file.");

   // Get output filename
   OUTFILENAME = string(ParticleType((Particle_t)PTYPE)) + ".hddm";
   gPARMS->SetDefaultParameter("OUTPUT_FILENAME", OUTFILENAME,
           "Filename of HDDM file to write particles to.");

   // Open output file
   if (hddmout == 0) {
      std::cout << " Error opening output file \"" << OUTFILENAME << "\"!"
                << std::endl;
      exit(-1);
   }
   std::cout << " output file: " << OUTFILENAME << std::endl;
   ofsout = new ofstream(OUTFILENAME.c_str());
   hddmout = new hddm_s::ostream(*ofsout);

   // Get vertex info
   string vertex_str = "0 0 65 65";
   gPARMS->SetDefaultParameter("VERTEX", vertex_str, 
           "Vertex to throw particles from (should be"
           " string of 4 numbers x y zmin zmax)");
   sscanf(vertex_str.c_str(), "%f %f %f %f", &vertex[0], &vertex[1],
                                             &vertex[2], &vertex[3]);
   if (vertex[2] > vertex[3]) {
      std::cerr << "Invalid parameter: z_min > z_max" << std::endl;
      exit(-1);
   }
   
   // Print message to user
   jout << std::endl;
   jout << "----------------------------------------" << std::endl;
   jout << "extract_ptype_hddm  plugin:" << std::endl;
   jout << std::endl;
   jout << "   particle type: " << ParticleType((Particle_t)PTYPE)
        << " (set via -PPTYPE=geantid)" << std::endl;
   jout << " output filename: " << OUTFILENAME 
        << " (set via -POUTFILENAME=fname.hddm)" << std::endl;
   jout << "          vertex: x=" << vertex[0] << " y=" << vertex[1]
        << " zmin=" << vertex[2] << " zmax=" << vertex[3] << std::endl;
   jout << "                  (set via -PVERTEX=\"X Y Zmin Zmax\")"
        << std::endl;
   jout << "----------------------------------------" << std::endl;
   jout << std::endl;

   // Unlock mutex
   pthread_mutex_unlock(&mutex);

   return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_extract_ptype_hddm::brun(JEventLoop *loop, 
                                                  int runnumber)
{
   // This is called whenever the run number changes
   return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_extract_ptype_hddm::evnt(JEventLoop *loop,
                                                  int eventnumber)
{
   vector<const DMCThrown*> mcthrowns;
   loop->Get(mcthrowns);

   pthread_mutex_lock(&mutex);

   for (unsigned int i=0; i < mcthrowns.size(); i++) {
      const DMCThrown *thrown = mcthrowns[i];
      
      if (thrown->type != (int)PTYPE)
         continue;
      
      // Start a new event
      hddm_s::HDDM record;
      hddm_s::PhysicsEventList pes = record.addPhysicsEvents();
      pes().setRunNo(1);
      pes().setEventNo(++Nevents);
      hddm_s::ReactionList rs = pes().addReactions();
      hddm_s::VertexList vs = rs().addVertices();
      hddm_s::OriginList os = vs().addOrigins();
      hddm_s::ProductList ps = vs().addProducts();
      os().setT(0.0);
      os().setVx(vertex[0]);
      os().setVy(vertex[1]);
      if (vertex[2] < vertex[3]) {
        os().setVz(randm(vertex[2],vertex[3]));
      }
      else {
        os().setVz(vertex[2]);
      }

      DVector3 mom = thrown->momentum();
      ps().setType((Particle_t)thrown->type);
      ps().setPdgtype(thrown->pdgtype);
      ps().setId(thrown->myid);
      ps().setParentid(thrown->parentid);
      ps().setMech(thrown->mech);
      hddm_s::MomentumList pmoms = ps().addMomenta();
      pmoms().setPx(mom.X());
      pmoms().setPy(mom.Y());
      pmoms().setPz(mom.Z());
      pmoms().setE(thrown->energy());
      
      *hddmout << record;
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
   if (hddmout) {
      delete hddmout;
      hddmout = NULL;
   }
   if (ofsout) {
      delete ofsout;
      ofsout = NULL;
   }
   pthread_mutex_unlock(&mutex);

   return NOERROR;
}

