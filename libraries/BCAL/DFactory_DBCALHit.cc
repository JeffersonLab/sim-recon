// $Id$
//
//    File: DFactory_DBCALHit.cc
// Created: Sun Apr  3 10:49:22 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include "DFactory_DBCALHit.h"

//------------------
// evnt
//------------------
derror_t DFactory_DBCALHit::evnt(int enventnumber)
{
	/// This just copies the data from the hddm_s structure for now

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	for(int i=0; i<PE->mult; i++){
		s_Modules_t *modules = NULL;
		s_HitView_t *HV = PE->in[i].hitView;
		if(PE->in[i].hitView)
			if(PE->in[i].hitView->barrelEMcal)
				modules = PE->in[i].hitView->barrelEMcal->modules;
		if(!modules)continue;

		for(int j=0;j<modules->mult;j++){
			float phim = modules->in[j].phim;
			s_Upstream_t *upstream = modules->in[j].upstream;
			if(upstream){
				s_Showers_t *showers = upstream->showers;
				if(!showers)continue;
				
				for(int m=0;m<showers->mult;m++){
					float E = showers->in[m].E;
					float t = showers->in[m].t;
					
					DBCALHit *bcalhit = new DBCALHit();
					bcalhit->phim = phim;
					bcalhit->end = DBCALHit::UPSTREAM;
					bcalhit->E = E;
					bcalhit->t = t;
					_data.push_back(bcalhit);
				}
			}
			s_Downstream_t *downstream = modules->in[j].downstream;
			if(downstream){
				s_Showers_t *showers = downstream->showers;
				if(!showers)continue;
				
				for(int m=0;m<showers->mult;m++){
					float E = showers->in[m].E;
					float t = showers->in[m].t;
					
					DBCALHit *bcalhit = new DBCALHit();
					bcalhit->phim = phim;
					bcalhit->end = DBCALHit::DOWNSTREAM;
					bcalhit->E = E;
					bcalhit->t = t;
					_data.push_back(bcalhit);
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
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("row:   phim(rad):      end:     E(GeV):   t(ns):");
	
	for(int i=0; i<_data.size(); i++){
		DBCALHit *bcalhit = _data[i];
		
		printnewrow();
		printcol("%d",	i);
		printcol("%1.3f",	bcalhit->phim);
		printcol(bcalhit->end==DBCALHit::UPSTREAM ? "upstream":"downstream");
		printcol("%2.3f",	bcalhit->E);
		printcol("%4.0f",	bcalhit->t);
		printrow();
	}

	return _table;
}
