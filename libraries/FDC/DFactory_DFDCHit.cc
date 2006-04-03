///
///	DFactory_DFDCHit.cc - Implementation for the basic FDC hit factory.
/// Author:		Craig Bookwalter (craigb at jlab.org)
/// Date:		March 2006
///

#include "DFactory_DFDCHit.h"
#include "DStreamLog.h"
#include <fstream>

DFactory_DFDCHit::DFactory_DFDCHit() {}

DFactory_DFDCHit::~DFactory_DFDCHit() {}

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
	s_PhysicsEvents_t* allEvents = hddm_s->physicsEvents;
	if(!allEvents) {
		throw DException("Attempt to get physics events from HDDM source failed.");
		return NOERROR;
	}
	
	for (unsigned int i=0; i < allEvents->mult; i++) {
		s_HitView_t *hits = allEvents->in[i].hitView;
		
		if (hits == HDDM_NULL) {
			throw DException("HDDM source has no hits.");
			return NOERROR;
		}
		
		if (hits->forwardDC == HDDM_NULL) {
			throw DException("HDDM source has no forwardDC information.");
			return NOERROR;
		}
		
		if (hits->forwardDC->fdcChambers == HDDM_NULL) {
			throw DException("HDDM source has no hits in the FDC.");		
			return NOERROR;
		}
		
		s_FdcChambers_t* fdcChamberSet = hits->forwardDC->fdcChambers;
		
		for (unsigned int i=0; i < fdcChamberSet->mult; i++) {
			s_FdcChamber_t fdcChamber 		= fdcChamberSet->in[i];		
			s_FdcAnodeWires_t* wireSet 		= fdcChamber.fdcAnodeWires;
			s_FdcCathodeStrips_t* stripSet 	= fdcChamber.fdcCathodeStrips;
		
			// Pull out anode data.
			for (unsigned int j=0; j < wireSet->mult; j++) {
				s_FdcAnodeWire_t anodeWire		= wireSet->in[j];
				s_FdcAnodeHits_t* wireHitSet	= anodeWire.fdcAnodeHits;
				for (unsigned int k=0; k < wireHitSet->mult; k++) {
					s_FdcAnodeHit_t wireHit		= wireHitSet->in[k];
					DFDCHit* newHit				= new DFDCHit();
					newHit->layer		 		= fdcChamber.layer;
					newHit->module		 		= fdcChamber.module;
					newHit->element				= anodeWire.wire;
					newHit->dE					= wireHit.dE;
					newHit->t					= wireHit.t;
					newHit->plane				= 2;
					newHit->type				= 0;
					newHit->globalPlane			= _geo.globalPlane(newHit);
					newHit->globalLayer			= _geo.globalLayer(newHit);
					newHit->r					= _geo.getWireR(newHit); 
					v.push_back(newHit);
				}
			}
			
			// Pull out cathode data.
			for (unsigned int j=0; j < stripSet->mult; j++) {
				s_FdcCathodeStrip_t cathodeStrip = stripSet->in[j];
				s_FdcCathodeHits_t* stripHitSet = cathodeStrip.fdcCathodeHits;
				for (unsigned int k=0; k < stripHitSet->mult; k++) {
					s_FdcCathodeHit_t stripHit	= stripHitSet->in[k];
					DFDCHit* newHit				= new DFDCHit();
					newHit->layer				= fdcChamber.layer;
					newHit->module				= fdcChamber.module;
					newHit->element				= cathodeStrip.strip;
					newHit->plane				= cathodeStrip.plane;
					newHit->dE					= stripHit.dE;
					newHit->t					= stripHit.t;
					newHit->type				= 1;
					newHit->globalPlane			= _geo.globalPlane(newHit);	 
					newHit->globalLayer			= _geo.globalLayer(newHit);
					newHit->r					= _geo.getStripR(newHit);
					v.push_back(newHit);
				}
			}	
		}
	}


	return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DFDCHit::toString(void)
{
/*	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("row: layer: module: tau(rad):    z(cm):  u(cm):  dE(MeV):   t(ns):   type:");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DFDCHit *fdchit = _data[i];
		
		printnewrow();
		printcol("%d",	i);
		printcol("%d", fdchit->layer);
		printcol("%d", fdchit->module);
		printcol("%3.1f", fdchit->tau);
		printcol("%3.1f", fdchit->z);
		printcol("%2.3f", fdchit->u);
//		if(!fdchit->type){
			printcol("%1.3f", fdchit->dE*1000.0);
			printcol("%4.0f", fdchit->t);
		}else{
			printcol("");
			printcol("");
		}
//		printcol("%s", fdchit->type ? "cathode":"anode");
		printrow();
	}

	return _table;
*/
	return "";
}
