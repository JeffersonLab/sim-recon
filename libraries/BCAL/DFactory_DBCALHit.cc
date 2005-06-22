// $Id$
//
//    File: DFactory_DBCALHit.cc
// Created: Thu Jun  9 10:14:35 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include "DFactory_DBCALHit.h"

//------------------
// evnt
//------------------
derror_t DFactory_DBCALHit::evnt(DEventLoop *eventLoop, int eventnumber)
{
	/// Place holder for now. 

	return NOERROR;
}

//------------------
// Extract_HDDM
//------------------
derror_t DFactory_DBCALHit::Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v)
{
	/// Copies the data from the given hddm_s structure. This is called
	/// from DEventSourceHDDM::GetObjects.
	
	v.clear();

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	for(unsigned int i=0; i<PE->mult; i++){
		s_Mods_t *mods = NULL;
		if(PE->in[i].hitView)
			if(PE->in[i].hitView->barrelEMcal)
				mods = PE->in[i].hitView->barrelEMcal->mods;
		if(!mods)continue;

		for(unsigned int j=0;j<mods->mult;j++){
			int module = mods->in[j].module;
			s_Shells_t *shells = mods->in[j].shells;
			for(unsigned int k=0;k<shells->mult;k++){
				int layer = shells->in[k].layer;
				s_Cones_t *cones = shells->in[k].cones;
				for(unsigned int m=0;j<cones->mult;m++){
					int sector = cones->in[m].sector;

					s_Upstream_t *upstream = cones->in[m].upstream;
					if(upstream){
						s_Showers_t *showers = upstream->showers;
						if(!showers)continue;
				
						for(unsigned int m=0;m<showers->mult;m++){
							float E = showers->in[m].E;
							float t = showers->in[m].t;
					
							DBCALHit *bcalhit = new DBCALHit();
							bcalhit->module = module;
							bcalhit->layer = layer;
							bcalhit->sector = sector;
							bcalhit->end = DBCALHit::UPSTREAM;
							bcalhit->E = E;
							bcalhit->t = t;
							v.push_back(bcalhit);
						}
					}
					s_Downstream_t *downstream = cones->in[m].downstream;
					if(downstream){
						s_Showers_t *showers = downstream->showers;
						if(!showers)continue;
				
						for(unsigned int m=0;m<showers->mult;m++){
						float E = showers->in[m].E;
							float t = showers->in[m].t;
					
							DBCALHit *bcalhit = new DBCALHit();
							bcalhit->module = module;
							bcalhit->layer = layer;
							bcalhit->sector = sector;
							bcalhit->end = DBCALHit::DOWNSTREAM;
							bcalhit->E = E;
							bcalhit->t = t;
							v.push_back(bcalhit);
						}
					}
				}
			}
		}
	}

	return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DBCALHit::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	GetNrows();
	if(_data.size()==0)return string(); // don't print anything if we have no data!

	printheader("row:   module:  layer:  sector:   end:     E(GeV):   t(ns):");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DBCALHit *bcalhit = _data[i];
		
		printnewrow();
		printcol("%d",	i);
		printcol("%d",	bcalhit->module);
		printcol("%d",	bcalhit->layer);
		printcol("%d",	bcalhit->sector);
		printcol(bcalhit->end==DBCALHit::UPSTREAM ? "upstream":"downstream");
		printcol("%2.3f",	bcalhit->E);
		printcol("%4.0f",	bcalhit->t);
		printrow();
	}

	return _table;
}
