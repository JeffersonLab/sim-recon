#include <fstream>
#include <iostream>
#include <map>
using namespace std;

#include <stdio.h>
#include <stdlib.h>

#include <particleType.h>
#include "bg_hddm.hpp"

std::ofstream *hddmOutputFile = NULL;
hddm_s::ostream *hddmOutputStream = NULL;

map<int, Particle_t> PDG_to_GEANT_map;
map<Particle_t, int> GEANT_to_PDG_map;
bool PDG_GEANT_maps_initialized = false;

void InitializePDGGEANTmaps(void);
void AlignParticleTypes(Particle &part);

//-----------------
// open_hddm_output
//-----------------
void open_hddm_output(string fname)
{
   // Open output file
   hddmOutputFile = new std::ofstream(fname.c_str());
   if (! hddmOutputFile->is_open()) {
      std::cerr << "Unable to open output file \"" << fname.c_str()
                << "\" for writing." << std::endl;
      exit(-3);
   }
   hddmOutputStream = new hddm_s::ostream(*hddmOutputFile);
   
   std::cerr << "Opened output file \"" << fname.c_str() 
             << "\" for writing." << std::endl;
}

//-----------------
// close_hddm_output
//-----------------
void close_hddm_output(void)
{
   // Close output file
   delete hddmOutputStream;
   delete hddmOutputFile;
   
   std::cout << "Closed HDDM output file." << std::endl;
}

//-----------------
// write_hddm_event
//-----------------
void write_hddm_event(Event &event)
{
   // User may specify either GEANT type or PDG type for each particle
   // When one of the types is zero and the other isn't, the non-zero
   // value is used to overwrite the zero value (if possible) os that
   // both are set properly.
   AlignParticleTypes(event.beam);
   AlignParticleTypes(event.target);

   // physicsEvent
   hddm_s::HDDM record;
   hddm_s::PhysicsEventList pes = record.addPhysicsEvents();
   pes().setRunNo(event.runNo);
   pes().setEventNo(event.eventNo);
   
   // reaction
   hddm_s::ReactionList rs = pes().addReactions();
   rs().setType(event.reaction_type);
   rs().setWeight(event.reaction_weight);

   // beam
   hddm_s::BeamList bs = rs().addBeams();
   bs().setType(event.beam.type);
   hddm_s::MomentumList bmoms = bs().addMomenta();
   bmoms().setPx(event.beam.p.Px());
   bmoms().setPy(event.beam.p.Py());
   bmoms().setPz(event.beam.p.Pz());
   bmoms().setE(event.beam.p.E());
   hddm_s::PropertiesList bpros = bs().addPropertiesList();
   bpros().setCharge(ParticleCharge(bs->type));
   bpros().setMass(ParticleMass(bs->type));

   // target
   hddm_s::TargetList ts = rs().addTargets();
   ts().setType(event.target.type);
   hddm_s::MomentumList tmoms = ts().addMomenta();
   tmoms().setPx(event.target.p.Px());
   tmoms().setPy(event.target.p.Py());
   tmoms().setPz(event.target.p.Pz());
   tmoms().setE(event.target.p.E());
   hddm_s::PropertiesList tpros = ts().addPropertiesList();
   tpros().setCharge(ParticleCharge(ts->type));
   tpros().setMass(ParticleMass(ts->type));

   // vertex
   hddm_s::VertexList vs = rs().addVertices();
   hddm_s::OriginList os = vs().addOrigins();
   os().setT(0.0);
   os().setVx(event.vertex.X());
   os().setVy(event.vertex.Y());
   os().setVz(event.vertex.Z());

   // product
   unsigned int Npart = event.intermediate.size() + event.final.size();
   hddm_s::ProductList ps = vs().addProducts(Npart);

   // Add intermediate particles (i.e. ones that GEANT does not track)
   // These will have the "type" explicitly set to 0 to tell hdgeant
   // to ignore them. They will, however, be kept in the list of thrown
   // particles passed to the output data stream of hdgeant.
   int nprod = 0;
   for (unsigned int i=0; i < event.intermediate.size(); i++) {
      AlignParticleTypes(event.intermediate[i]);
      CopyParticleToProduct(nprod+1, event.intermediate[i], ps.end());
      ps(nprod).setType((Particle_t)0);
      ++nprod;
   }

   // Add final state particles (i.e. ones that GEANT does track)
   for (unsigned int i=0; i < event.final.size(); i++) {
      AlignParticleTypes(event.final[i]);
      CopyParticleToProduct(nprod+1, event.final[i], ps.end());
      ++nprod;
   }

   // Write event to output stream
   *hddmOutputStream << record;
}

//-----------------
// CopyParticleToProduct
//-----------------
void CopyParticleToProduct(int id, const Particle &part, PlistIter &prod)
{
   prod.decayVertex = part.decayVertex;
   prod.id = id;
   prod.mech = part.mech;
   prod.parentid = part.parentid;
   prod.pdgtype = part.pdgtype;
   prod.type = part.type;

   // momentum
   hddm_s::MomentumList pmoms = prod->addMomenta();
   pmoms().setPx(part.p.Px());
   pmoms().setPy(part.p.Py());
   pmoms().setPz(part.p.Pz());
   pmoms().setE(part.p.E());

   // properties
   hddm_s::PropertiesList ppros = prod->addPropertiesList();
   ppros().setCharge(ParticleCharge(prod.type));
   ppros().setMass(ParticleMass(prod.type));
}

//-----------------
// PDG_to_GEANT
//-----------------
Particle_t PDG_to_GEANT(int pdgtype)
{
   if (! PDG_GEANT_maps_initialized)
      InitializePDGGEANTmaps();

   map<int, Particle_t>::iterator iter = PDG_to_GEANT_map.find(pdgtype);
   if (iter == PDG_to_GEANT_map.end())
      pdgtype = 0;

   return PDG_to_GEANT_map[pdgtype];
}

//-----------------
// GEANT_to_PDG
//-----------------
int GEANT_to_PDG(Particle_t type)
{
   if (! PDG_GEANT_maps_initialized)
      InitializePDGGEANTmaps();

   map<Particle_t, int>::iterator iter = GEANT_to_PDG_map.find(type);
   if (iter == GEANT_to_PDG_map.end())
      type = Unknown;

   return GEANT_to_PDG_map[type];
}

//-----------------
// InitializePDGGEANTmap
//-----------------
void InitializePDGGEANTmaps(void)
{
   // Set values in GEANT_to_PDG_map first, then copy them into PDG_to_GEANT_map
   GEANT_to_PDG_map[Unknown]   = 0;
   GEANT_to_PDG_map[Gamma]     = 22;
   GEANT_to_PDG_map[Positron]  = -11;
   GEANT_to_PDG_map[Electron]  = 11;
   GEANT_to_PDG_map[Neutrino]  = 12;
   GEANT_to_PDG_map[MuonPlus]  = -13;
   GEANT_to_PDG_map[MuonMinus] = 13;
   GEANT_to_PDG_map[Pi0]       = 111;
   GEANT_to_PDG_map[PiPlus]    = 211;
   GEANT_to_PDG_map[PiMinus]   = -211;
   GEANT_to_PDG_map[KLong]     = 130;
   GEANT_to_PDG_map[KPlus]     = 321;
   GEANT_to_PDG_map[KMinus]    = -321;
   GEANT_to_PDG_map[Neutron]   = 2112;
   GEANT_to_PDG_map[Proton]    = 2212;
   GEANT_to_PDG_map[AntiProton]= -2212;
   GEANT_to_PDG_map[KShort]    = 310;
   GEANT_to_PDG_map[Eta]       = 221;
   GEANT_to_PDG_map[Lambda]    = 3122;

   // Copy values into PDG_to_GEANT_map
   map<Particle_t, int>::iterator iter = GEANT_to_PDG_map.begin();
   for (; iter != GEANT_to_PDG_map.end(); iter++)
      PDG_to_GEANT_map[iter->second] = iter->first;

   PDG_GEANT_maps_initialized = true;
}

//-----------------
// InitializePDGGEANTmap
//-----------------
void AlignParticleTypes(Particle &part)
{
   // If one of type or pdgtype is zero and the other isn't
   // then use the non-zero value to set the zero one (if
   // possible).
   if (part.type == 0 && part.pdgtype != 0) {
      part.type = PDG_to_GEANT(part.pdgtype);
   }
   else if (part.type != 0 && part.pdgtype == 0) {
      part.pdgtype = GEANT_to_PDG(part.type);
   }
}
