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

class DParticleSet:public jana::JObject{
 public:
  JOBJECT_PUBLIC(DParticleSet);
  
  const DVertex *vertex;
  vector<const DVertex::shower_info_t *>photon;
  vector<const DVertex::track_info_t *>pip;  // list of pi pluses
  vector<const DVertex::track_info_t *>pim;  // list of pi minuses
  vector<const DVertex::track_info_t *>Kp; // list of K pluses
  vector<const DVertex::track_info_t *>Km;  // list of K minuses
  vector<const DVertex::track_info_t *>proton; // list of protons

  // Print out some summary information about the contents of this class
  void toStrings(vector<pair<string,string> > &items)const{
    AddString(items, "x", "%3.2f", vertex->x.X());
    AddString(items, "y", "%3.2f", vertex->x.Y());
    AddString(items, "z", "%3.2f", vertex->x.Z());
    AddString(items, "t", "%3.2f", vertex->x.T());
    AddString(items, "Nphoton",      "%d", photon.size());
    AddString(items, "Npi_plus",     "%d", pip.size());
    AddString(items, "Npi_minus",    "%d", pim.size());
    AddString(items, "Nproton",      "%d", proton.size());
    AddString(items, "NK_plus",      "%d", Kp.size());
    AddString(items, "NK_minus",     "%d", Km.size());    
	
  }
};

#endif // _DParticleSet_

