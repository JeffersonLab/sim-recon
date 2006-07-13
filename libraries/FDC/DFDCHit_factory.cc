//***********************************************************************
// DFDCHit_factory.cc - Implementation for the basic FDC hit factory.
// Author:		Craig Bookwalter (craigb at jlab.org)
// Date:		March 2006
//***********************************************************************

#include "DFDCHit_factory.h"

///
/// DFDCHit_factory::DFDCHit_factory():
/// default constructor -- empty, for now.
///
DFDCHit_factory::DFDCHit_factory() {}

///
/// DFDCHit_factory::~DFDCHit_factory()
/// default destructor -- also empty for now.
///
DFDCHit_factory::~DFDCHit_factory() {}

///
/// DFDCHit_factory::evnt():
/// This would be used if this factory was going to process data currently exisiting in
/// memory; since this factory's role is to read events from the data file and put them
/// into memory, this is not used. See DFDCHit_factory::Extract_HDDM().
///
jerror_t DFDCHit_factory::evnt(JEventLoop *eventLoop, int eventnumber)
{
	return NOERROR;
}

///
/// DFDCHit_factory::Extract_HDDM():
/// Reads an event from the data file and distills out the FDC information, creating a new
/// DFDCHit object for each FDC hit encountered in the data. If you wish to understand the 
/// s_Blah_t structures, see the documentation for hddm_s.h.
///
jerror_t DFDCHit_factory::Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v)
{		
	v.clear();
	
	// Acquire the pointer to the physics events
	s_PhysicsEvents_t* allEvents = hddm_s->physicsEvents;
	if(!allEvents) {
		throw JException("Attempt to get physics events from HDDM source failed.");
		return NOERROR;
	}
	
	for (unsigned int i=0; i < allEvents->mult; i++) {
		// Acquire the pointer to the overall hits section of the data
		s_HitView_t *hits = allEvents->in[i].hitView;
		
		if (hits == HDDM_NULL) {
			throw JException("HDDM source has no hits.");
			return NOERROR;
		}
		
		if (hits->forwardDC == HDDM_NULL) {
			throw JException("HDDM source has no forwardDC information.");
			return NOERROR;
		}
		
		if (hits->forwardDC->fdcChambers == HDDM_NULL) {
			throw JException("HDDM source has no hits in the FDC.");		
			return NOERROR;
		}
		
		// Acquire the pointer to the beginning of the FDC hit tree
		s_FdcChambers_t* fdcChamberSet = hits->forwardDC->fdcChambers;
		
		for (unsigned int i=0; i < fdcChamberSet->mult; i++) {
			// Each chamber in the ChamberSet has a wire set and a strip set
			s_FdcChamber_t fdcChamber 		= fdcChamberSet->in[i];		
			s_FdcAnodeWires_t* wireSet 		= fdcChamber.fdcAnodeWires;
			s_FdcCathodeStrips_t* stripSet 	= fdcChamber.fdcCathodeStrips;
		
			// Each set of wires has (obviously) wires inside of it, and each wire
			// may have one or more hits on it. Make a DFDCHit object for each one
			// of these hits.
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
					newHit->gPlane				= _geo.gPlane(newHit);
					newHit->gLayer				= _geo.gLayer(newHit);
					newHit->r					= _geo.getWireR(newHit); 
					v.push_back(newHit);
				}
			}
		
			// Ditto for the cathodes.
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
					newHit->gPlane				= _geo.gPlane(newHit);	 
					newHit->gLayer				= _geo.gLayer(newHit);
					newHit->r					= _geo.getStripR(newHit);
					v.push_back(newHit);
				}
			}	
		}
	}

	return NOERROR;
}

///
/// DFDCHit_factory::toString(): 
/// Provides a sensible std::string representation of all of the data in the factory.
///
const string DFDCHit_factory::toString(void)
{
	Get();
	if( _data.size() <= 0) 
		return ""; 
		
	stringstream s;
	
	s << (*_data.begin())->header() << endl;
	
	// Simply call the toString() method of each DFDCHit object and stream it into s.
	for (unsigned int i=0; i < _data.size(); ++i)
		s << _data[i]->toString() << endl;
	
	return s.str();
}
	
