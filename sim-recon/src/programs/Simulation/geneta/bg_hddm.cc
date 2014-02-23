
#include <iostream>
#include <map>
using namespace std;

#include <stdio.h>

#include <particleType.h>
#include "bg_hddm.h"

s_iostream_t* hddmOutputStream=NULL;

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
	hddmOutputStream = init_s_HDDM((char*)fname.c_str());
	if(!hddmOutputStream){
		cerr<<"Unable to open output file \""<<fname.c_str()<<"\" for writing."<<endl;
		exit(-3);
	}
	
	cerr<<"Opened output file \""<<fname.c_str()<<"\" for writing."<<endl;
}

//-----------------
// close_hddm_output
//-----------------
void close_hddm_output(void)
{
	// Close output file
	close_s_HDDM(hddmOutputStream);
	
	cout<<"Closed HDDM output file."<<endl;
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
	s_HDDM_t *thisOutputEvent = make_s_HDDM();
	s_PhysicsEvents_t* pes = thisOutputEvent->physicsEvents = make_s_PhysicsEvents(1);
	pes->mult = 1;
	pes->in[0].runNo   = event.runNo;
	pes->in[0].eventNo = event.eventNo;
	
		// reaction
		s_Reactions_t* rs = pes->in[0].reactions = make_s_Reactions(1);
		rs->mult = 1;
		rs->in[0].type = event.reaction_type;
		rs->in[0].weight = event.reaction_weight;

			// beam
			s_Beam_t* bs = rs->in[0].beam = make_s_Beam();
			bs->type = event.beam.type;
				s_Momentum_t *mom = bs->momentum = make_s_Momentum();
				mom->px = event.beam.p.Px();
				mom->py = event.beam.p.Py();
				mom->pz = event.beam.p.Pz();
				mom->E  = event.beam.p.E();
					s_Properties_t *prop = bs->properties = make_s_Properties();
					prop->charge = ParticleCharge(bs->type);
					prop->mass = ParticleMass(bs->type);
        
			// target
			s_Target_t* ts = rs->in[0].target = make_s_Target();
			ts->type = event.target.type;
				mom = ts->momentum = make_s_Momentum();
				mom->px = event.target.p.Px();
				mom->py = event.target.p.Py();
				mom->pz = event.target.p.Pz();
				mom->E  = event.target.p.E();
					prop = ts->properties = make_s_Properties();
					prop->charge = ParticleCharge(ts->type);
					prop->mass = ParticleMass(ts->type);

			// vertex
			s_Vertices_t* vs = rs->in[0].vertices = make_s_Vertices(1);
			vs->mult = 1;
				
				// origin
				s_Origin_t* origin = vs->in[0].origin = make_s_Origin();
				origin->t = 0.0;
				origin->vx = event.vertex.X();
				origin->vy = event.vertex.Y();
				origin->vz = event.vertex.Z();

				// product
				unsigned int Npart = event.intermediate.size() + event.final.size();
				s_Products_t* ps = vs->in[0].products = make_s_Products(Npart);
				ps->mult = 0; // set to 0 so we can increment particles are added
	
				// Add intermediate particles (i.e. ones that GEANT does not track)
				// These will have the "type" explicitly set to 0 to tell hdgeant
				// to ignore them. They will, however, be kept in the list of thrown
				// particles passed to the output data stream of hdgeant.
				for(unsigned int i=0; i<event.intermediate.size(); i++){
					AlignParticleTypes(event.intermediate[i]);
					CopyParticleToProduct(ps->mult+1, event.intermediate[i], ps->in[ps->mult]);
					ps->in[ps->mult].type = (Particle_t)0;
					ps->mult++;
				}

				// Add final state particles (i.e. ones that GEANT does track)
				for(unsigned int i=0; i<event.final.size(); i++){
					AlignParticleTypes(event.final[i]);
					CopyParticleToProduct(ps->mult+1, event.final[i], ps->in[ps->mult]);
					ps->mult++;
				}
	
	// Write event to output stream
	flush_s_HDDM(thisOutputEvent, hddmOutputStream);
}

//-----------------
// CopyParticleToProduct
//-----------------
void CopyParticleToProduct(int id, const Particle &part, s_Product_t &prod)
{
	prod.decayVertex = part.decayVertex;
	prod.id = id;
	prod.mech = part.mech;
	prod.parentid = part.parentid;
	prod.pdgtype = part.pdgtype;
	prod.type = part.type;

		// momentum
		prod.momentum = make_s_Momentum();
		prod.momentum->px = part.p.Px();
		prod.momentum->py = part.p.Py();
		prod.momentum->pz = part.p.Pz();
		prod.momentum->E  = part.p.E();

		// properties
		prod.properties = make_s_Properties();
		prod.properties->charge = ParticleCharge(prod.type);
		prod.properties->mass = ParticleMass(prod.type);
}

//-----------------
// PDG_to_GEANT
//-----------------
Particle_t PDG_to_GEANT(int pdgtype)
{
	if(!PDG_GEANT_maps_initialized)InitializePDGGEANTmaps();

	map<int, Particle_t>::iterator iter = PDG_to_GEANT_map.find(pdgtype);
	if(iter==PDG_to_GEANT_map.end())pdgtype = 0;

	return PDG_to_GEANT_map[pdgtype];
}

//-----------------
// GEANT_to_PDG
//-----------------
int GEANT_to_PDG(Particle_t type)
{
	if(!PDG_GEANT_maps_initialized)InitializePDGGEANTmaps();

	map<Particle_t, int>::iterator iter = GEANT_to_PDG_map.find(type);
	if(iter==GEANT_to_PDG_map.end())type = Unknown;

	return GEANT_to_PDG_map[type];
}

//-----------------
// InitializePDGGEANTmap
//-----------------
void InitializePDGGEANTmaps(void)
{
	// Set values in GEANT_to_PDG_map first, then copy them into PDG_to_GEANT_map
	GEANT_to_PDG_map[Unknown]	= 0;
	GEANT_to_PDG_map[Gamma]		= 22;
	GEANT_to_PDG_map[Positron] = -11;
	GEANT_to_PDG_map[Electron] = 11;
	GEANT_to_PDG_map[Neutrino] = 12;
	GEANT_to_PDG_map[MuonPlus] = -13;
	GEANT_to_PDG_map[MuonMinus] = 13;
	GEANT_to_PDG_map[Pi0]		= 111;
	GEANT_to_PDG_map[PiPlus]	= 211;
	GEANT_to_PDG_map[PiMinus]	= -211;
	GEANT_to_PDG_map[KLong]		= 130;
	GEANT_to_PDG_map[KPlus]		= 321;
	GEANT_to_PDG_map[KMinus]	= -321;
	GEANT_to_PDG_map[Neutron]	= 2112;
	GEANT_to_PDG_map[Proton]	= 2212;
	GEANT_to_PDG_map[AntiProton] = -2212;
	GEANT_to_PDG_map[KShort]	= 310;
	GEANT_to_PDG_map[Eta]		= 221;
	GEANT_to_PDG_map[Lambda]	= 3122;

	// Copy values into PDG_to_GEANT_map
	map<Particle_t, int>::iterator iter = GEANT_to_PDG_map.begin();
	for(; iter!=GEANT_to_PDG_map.end(); iter++)PDG_to_GEANT_map[iter->second] = iter->first;

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
	if(part.type==0 && part.pdgtype!=0){
		part.type = PDG_to_GEANT(part.pdgtype);
	}else if(part.type!=0 && part.pdgtype==0){
		part.pdgtype = GEANT_to_PDG(part.type);
	}
}
