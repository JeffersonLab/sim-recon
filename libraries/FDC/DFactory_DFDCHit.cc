// $Id$
//
//    File: DFactory_DFDCHit.cc
// Created: Thu Jun  9 10:25:22 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include "DFactory_DFDCHit.h"

// z-positions of modules (3 chambers each)
float fdc_zmodule[8]={231.0,237.0,286.0,292.0,341.0,347.0,396.0,402.0};


//------------------
// evnt
//------------------
derror_t DFactory_DFDCHit::evnt(DEventLoop *eventLoop, int eventnumber)
{
	/// Place holder for now. 

	return NOERROR;
}

//------------------
// Extract_HDDM
//------------------
derror_t DFactory_DFDCHit::Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v)
{
	/// Copies the data from the given hddm_s structure. This is called
	/// from DEventSourceHDDM::GetObjects.
	
	v.clear();

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;

	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL ||
			hits->forwardDC == HDDM_NULL ||
			hits->forwardDC->fdcChambers == HDDM_NULL)continue;

		s_FdcChambers_t* fdcChambers = hits->forwardDC->fdcChambers;
		for(unsigned int j=0; j<fdcChambers->mult; j++){
			s_FdcChamber_t *chamber = &fdcChambers->in[j];
			int  layer = chamber->layer;
			int  module = chamber->module;
			
			if(chamber->fdcAnodeWires != HDDM_NULL){
				for(unsigned int k=0; k<chamber->fdcAnodeWires->mult; k++){ 
				  s_FdcAnodeHits_t *hits=chamber->fdcAnodeWires->in[k].fdcAnodeHits;
				  for (unsigned int n=0;n<hits->mult;n++){
				    DFDCHit *fdchit = new DFDCHit;
				    
				    fdchit->layer=layer;
				    fdchit->module=module;
				    fdchit->plane=2;
				    fdchit->u=chamber->fdcAnodeWires->in[k].wire;
				    fdchit->dE=hits->in[n].dE;
				    fdchit->t=hits->in[n].t;
				    fdchit->z=fdc_zmodule[module-1]+2.0*float(layer-2);
				    fdchit->tau=60.0*float(layer-2);
				    fdchit->type = 0;

				    v.push_back(fdchit);
				  }//  n (anode hits) 
				} //  k  (fdcAnodeWires)
			} // fdcAnodeWires!=NULL

			if(chamber->fdcCathodeStrips != HDDM_NULL){
			  for(unsigned int k=0; k<chamber->fdcCathodeStrips->mult; k++){ 
			    s_FdcCathodeHits_t *hits=chamber->fdcCathodeStrips->in[k].fdcCathodeHits;
			    for (unsigned int n=0;n<hits->mult;n++){
			      DFDCHit *fdchit = new DFDCHit;
			      int plane=chamber->fdcCathodeStrips->in[k].plane;
			      
			      fdchit->layer=layer;
			      fdchit->module=module;
			      fdchit->plane=plane;
			      fdchit->u=chamber->fdcCathodeStrips->in[k].strip;
			      fdchit->dE=hits->in[n].dE;
			      fdchit->t=hits->in[n].t;
			      fdchit->z=fdc_zmodule[module-1]+2.0*float(layer-2);
			      fdchit->tau=60.0*float(layer-2)+45.0*float(plane-2);
			      fdchit->type = 1;

			      v.push_back(fdchit);
			    }//  n (cathode hits) 
			  } //  k  (fdcCathodeStrips)
			} // fdcAnodeWires!=NULL

		} // j (fdcChambers)
	} // i  (physicsEvents)

	return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DFDCHit::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("row: layer: module: tau(rad):    z(cm):  local index:  dE(MeV):   t(ns):   type:");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DFDCHit *fdchit = _data[i];
		
		printnewrow();
		printcol("%d",	i);
		printcol("%d", fdchit->layer);
		printcol("%d", fdchit->module);
		printcol("%3.1f", fdchit->tau);
		printcol("%3.1f", fdchit->z);
		printcol("%d", fdchit->u);
		if(!fdchit->type){
			printcol("%1.3f", fdchit->dE*1000.0);
			printcol("%4.0f", fdchit->t);
		}else{
			printcol("");
			printcol("");
		}
		printcol("%s", fdchit->type ? "cathode":"anode");
		printrow();
	}

	return _table;
}
