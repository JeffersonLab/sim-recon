// $Id$
//
//    File: DFactory_DCDCHit.cc
// Created: Sun Apr  3 10:46:28 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include "DFactory_DCDCHit.h"

//------------------
// evnt
//------------------
derror_t DFactory_DCDCHit::evnt(int enventnumber)
{
	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	for(int i=0; i<PE->mult; i++){
		// ------------ CdcPoints, Hits --------------
		s_Rings_t *rings=NULL;
		s_HitView_t *HV = PE->in[i].hitView;
		if(PE->in[i].hitView)
			if(PE->in[i].hitView->centralDC)
				rings = PE->in[i].hitView->centralDC->rings;
		if(rings){
			for(int j=0;j<rings->mult;j++){
				float radius = rings->in[j].radius;
				s_Straws_t *straws = rings->in[j].straws;
				if(straws){
					for(int k=0;k<straws->mult;k++){
						float phim = straws->in[k].phim;
						s_Hits_t *hits = straws->in[k].hits;
						if(hits){
							for(int m=0;m<hits->mult;m++){
								float dE = hits->in[m].dE;
								float t = hits->in[m].t;
								
								// Add a row to the factory data
								DCDCHit *cdchit = new DCDCHit;
								cdchit->radius = radius;
								cdchit->phim = phim;
								cdchit->dE = dE;
								cdchit->t = t;
								_data.push_back(cdchit);
							}
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
const string DFactory_DCDCHit::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("row: radius(cm):  phim(rad):   dE(MeV):   t(ns):");
	
	for(int i=0; i<_data.size(); i++){
		DCDCHit *cdchit = _data[i];

		printnewrow();
		printcol("%d",	i);
		printcol("%3.1f", cdchit->radius);
		printcol("%1.3f", cdchit->phim);
		printcol("%2.3f", cdchit->dE);
		printcol("%4.0f", cdchit->t);
		printrow();
	}

	return _table;
}
