// $Id$
//
//    File: DParticleSet_factory.cc
// Created: Tue Mar 15 11:17:35 EDT 2011
// Creator: staylor (on Linux ifarml1 2.6.18-128.el5 x86_64)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DParticleSet_factory.h"
#include <TRACKING/DTrackTimeBased.h>
using namespace jana;

//------------------
// init
//------------------
jerror_t DParticleSet_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DParticleSet_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DParticleSet_factory::evnt(JEventLoop *loop, int eventnumber)
{
  // Get vertices
  vector<const DVertex *>vertices;
  loop->Get(vertices);

  // Sort according to particle type
  for (unsigned int i=0;i<vertices.size();i++){
    DParticleSet *particle_set = new DParticleSet;
    particle_set->vertex=vertices[i];

    // Charged particles:  take the first entry in each hypothesis list; this entry has the 
    // highest figure of merit.
    for (unsigned int k=0;k<vertices[i]->hypotheses.size();k++){
      const DVertex::track_info_t *track_info=&vertices[i]->hypotheses[k][0];
      const DTrackTimeBased *track=track_info->track;
      double mass=track->mass();
      
      // Create a list of pointers to the results for the various hypotheses for this track
      vector<const DVertex::track_info_t*>hypothesis_list;
      for (unsigned int m=0;m<vertices[i]->hypotheses[k].size();m++){
	hypothesis_list.push_back(&vertices[i]->hypotheses[k][m]);
      }

      // Deal with positive particles
      if (track->charge()>0){
	if (mass<0.3) particle_set->pip.push_back(hypothesis_list);
	else if (mass>=0.3 && mass<0.7) particle_set->Kp.push_back(hypothesis_list);
	else if (mass>=0.7 && mass<1.0) particle_set->proton.push_back(hypothesis_list);
	
      }
      else{ // Deal with negative particles
	if (mass<0.3) particle_set->pim.push_back(hypothesis_list);
	else if (mass>=0.3 && mass<0.7) particle_set->Km.push_back(hypothesis_list);
      }
    }
    // Now deal with photons
    for (unsigned int k=0;k<vertices[i]->showers.size();k++){
      if (vertices[i]->showers[k].matched_track==NULL){
	particle_set->photon.push_back(&vertices[i]->showers[k]);
      }
    }


    _data.push_back(particle_set);
  }

  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DParticleSet_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DParticleSet_factory::fini(void)
{
  return NOERROR;
}

