// $Id$
//
//    File: DParticleSet.h
// Created: Tue Mar 15 11:17:35 EDT 2011
// Creator: staylor (on Linux ifarml1 2.6.18-128.el5 x86_64)
//

#ifndef _DParticleSet_
#define _DParticleSet_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include <PID/DVertex.h>
#include <PID/DChargedTrack.h>
#include <PID/DNeutralParticle.h>

class DParticleSet : public jana::JObject{
	public:
		JOBJECT_PUBLIC(DParticleSet);
  
		const DVertex *vertex;

		vector<const DChargedTrack*> pip; // list of pi pluses
		vector<const DChargedTrack*> pim; // list of pi minuses
		vector<const DChargedTrack*> Kp; // list of K pluses
		vector<const DChargedTrack*> Km; // list of K minuses
		vector<const DChargedTrack*> proton; // list of protons
		vector<const DChargedTrack*> otherp; // unidentified positively charged particles
		vector<const DChargedTrack*> othern; // unidentified negatively charged particles

		vector<const DNeutralParticle*> photon; //list of photons
		vector<const DNeutralParticle*> neutron; //list of photons
		vector<const DNeutralParticle*> otherz; //unidentified neutrally charged particles

		vector<const DChargedTrack*> dChargedTracks; // list of all charged tracks associated with this DVertex
		vector<const DNeutralParticle*> dNeutralParticles; // list of all neutral tracks associated with this DVertex

  // Print out some summary information about the contents of this class
  void toStrings(vector<pair<string,string> > &items)const{
    AddString(items, "x", "%3.2f", vertex->dSpacetimeVertex.X());
    AddString(items, "y", "%3.2f", vertex->dSpacetimeVertex.Y());
    AddString(items, "z", "%3.2f", vertex->dSpacetimeVertex.Z());
    AddString(items, "t", "%3.2f", vertex->dSpacetimeVertex.T());
    AddString(items, "Nphoton",      "%d", photon.size());
    AddString(items, "Nneutron",      "%d", neutron.size());
    AddString(items, "Npi_plus",     "%d", pip.size());
    AddString(items, "Npi_minus",    "%d", pim.size());
    AddString(items, "Nproton",      "%d", proton.size());
    AddString(items, "NK_plus",      "%d", Kp.size());
    AddString(items, "NK_minus",     "%d", Km.size());    
    AddString(items, "Notherp",      "%d",  otherp.size());
    AddString(items, "Nothern",      "%d",  othern.size());
    AddString(items, "Notherz",      "%d",  otherz.size());
	
  }
};

#endif // _DParticleSet_

