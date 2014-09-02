#include <stdlib.h>
#include <stdio.h>
#include <fstream>

#include "HDDM/hddm_s.hpp"

std::ofstream *fileOutputStream = NULL;
hddm_s::ostream *hddmOutputStream = NULL;

typedef struct {
        int geantid;
        int mech; /* what do the values of this correspond to */
        int kfid;
        int parent;
        int firstdaughter;
        int lastdaughter;
} keve_t;

typedef struct {
   float px;
   float py;
   float pz;
   float en;
} peve_t;

/*-----------------
// open_hddm_output_
//-----------------*/
void open_hddm_output(std::string outputfile)
{
   /* Open output file */
   fileOutputStream = new std::ofstream(outputfile.c_str());
   if (! fileOutputStream->is_open()) {
      fprintf(stderr, "Unable to open output file \"%s\" for writing.\n", 
              outputfile.c_str());
      exit(-3);
   }
   hddmOutputStream = new hddm_s::ostream(*fileOutputStream);
   printf("Opened HDDM file \"%s\" for writing ...\n", outputfile.c_str());
}

/*-----------------
// close_hddm_output_
//-----------------*/
void close_hddm_output(void)
{
   /* Close output file */
   delete hddmOutputStream;
   delete fileOutputStream;
   
   printf("Closed HDDM output file\n");
}

/*-----------------
// write_hddm_event_
//-----------------*/
void write_hddm_event(int *iev, int *iproc,
                       keve_t *kin,  peve_t *pin,   
            int *ntra, keve_t *keve, peve_t *peve)
{
   /* Loop over events */
   int i;
   static int Nevents = 0;
   static int Nevents_written = 0;
   int runNumber=2;
   float vertex[3]={0.0, 0.0, 65.0}; 

   Nevents++;

   /* Start a new event */
   hddm_s::HDDM record;
   hddm_s::PhysicsEventList pes = record.addPhysicsEvents();
   pes().setRunNo(runNumber);
   pes().setEventNo(Nevents);
   hddm_s::ReactionList rs = pes().addReactions();
   rs().setType(*iproc);
   hddm_s::BeamList bs = rs().addBeams();
   bs().setType((Particle_t)kin[0].geantid);
   hddm_s::MomentumList bmoms = bs().addMomenta();
   bmoms().setPx(pin[0].px);
   bmoms().setPy(pin[0].py);
   bmoms().setPz(pin[0].pz);
   bmoms().setE(pin[0].en);
   hddm_s::PropertiesList bpros = bs().addPropertiesList();
   bpros().setCharge(0.0);
   bpros().setMass(0.0);
        
   hddm_s::TargetList ts = rs().addTargets();
   ts().setType((Particle_t)kin[1].geantid);
   hddm_s::MomentumList tmoms = ts().addMomenta();
   tmoms().setPx(pin[1].px);
   tmoms().setPy(pin[1].py);
   tmoms().setPz(pin[1].pz);
   tmoms().setE(pin[1].en);
   hddm_s::PropertiesList tpros = ts().addPropertiesList();
   tpros().setCharge(+1);
   tpros().setMass(0.938272); /* this should be derived from type ... */
        
   hddm_s::VertexList vs = rs().addVertices();
   hddm_s::OriginList os = vs().addOrigins();
   hddm_s::ProductList ps = vs().addProducts(*ntra);
   
   os().setT(0.0);
   os().setVx(vertex[0]);
   os().setVy(vertex[1]);
   os().setVz(vertex[2]);
   
   for (i=0; i < *ntra; i++) {
      //double E2;
      //if (keve[i].geantid == 0)
      //   continue;
      
      ps(i).setType((Particle_t)keve[i].geantid);
      ps(i).setMech(keve[i].mech);
      ps(i).setPdgtype(keve[i].kfid);
      ps(i).setId(i+1);
      ps(i).setParentid(keve[i].parent);
      hddm_s::MomentumList pmoms = ps(i).addMomenta();
      pmoms().setPx(peve[i].px);
      pmoms().setPy(peve[i].py);
      pmoms().setPz(peve[i].pz);
      pmoms().setE(peve[i].en);
   }
   
   if (*ntra > 0) {
      Nevents_written++;
      *hddmOutputStream << record;
      if (Nevents_written%10000 == 0)
         printf("Wrote event %d events (%d generated)\n", 
                Nevents_written, Nevents);
   }   
}
