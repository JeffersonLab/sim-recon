//*****************************************************************************
// DFDCTruth_factory.cc - implementation for pulling truth data out of HDDM files
// Author:		Craig Bookwalter (craigb at jlab.org)
// Date:		March 2006
//*****************************************************************************

#include "DFDCTruth_factory.h"

///
/// DFDCTruth_factory::DFDCTruth_factory():
/// default constructor -- empty for now
///
DFDCTruth_factory::DFDCTruth_factory() {}

///
/// DFDCTruth_factory::~DFDCTruth_factory():
/// default destructor -- empty for now as well
///
DFDCTruth_factory::~DFDCTruth_factory() {}

///
/// DFDCTruth_factory::evnt():
/// this function would be defined if this factory was processing data
/// already in memory; since this factory reads an HDDM file and puts
/// truth data into memory, the DFDCTruth_factory::Extract_HDDM() method
/// is used.
///
jerror_t DFDCTruth_factory::evnt(JEventLoop *eventLoop, int eventnumber)
{
	return NOERROR;
}

///
/// DFDCTruth_factory::Extract_HDDM():
/// reads an event from a data file and converts FDC truth information into 
/// DFDCTruth objects. For more detail on the s_Blah_t structures, see
/// the documentation for hddm_s.h.
///
jerror_t DFDCTruth_factory::Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v)
{	
	v.clear();
	
	// Acquire the pointer to the entire event
	s_PhysicsEvents_t* allEvents = hddm_s->physicsEvents;
	if(!allEvents) {
		throw JException("Attempt to get physics events from HDDM source failed.");
		return NOERROR;
	}
	
	for (unsigned int i=0; i < allEvents->mult; i++) {
		// Acquire the pointer to all of the hits for this event
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
		
		// Acquire the pointer to the top of the FDC hit tree
		s_FdcChambers_t* fdcChamberSet = hits->forwardDC->fdcChambers;
		
		for (unsigned int i=0; i < fdcChamberSet->mult; i++) {
			// For each chamber, extract its truth points and turn them into 
			// DFDCTruth objects.
			s_FdcChamber_t fdcChamber = fdcChamberSet->in[i];
			s_FdcTruthPoints_t* fdcTruthPoints = fdcChamber.fdcTruthPoints;
			if (fdcTruthPoints == HDDM_NULL)
				continue;
			for (unsigned int j=0; j < fdcTruthPoints->mult; j++) {
				s_FdcTruthPoint_t truthPoint = fdcTruthPoints->in[j];
				DFDCTruth* newPoint 	= new DFDCTruth();
				newPoint->x				= truthPoint.x;
				newPoint->y				= truthPoint.y;
				newPoint->z				= truthPoint.z;
				newPoint->track			= truthPoint.track;
				newPoint->primary		= truthPoint.primary;
				newPoint->dradius		= truthPoint.dradius;
				newPoint->dEdx			= truthPoint.dEdx;
				v.push_back(newPoint);
			}
		}
	}		
	return NOERROR;
}

///
/// DFDCTruth_factory::toString():
/// Print a sensible std::string representation of the contents of this factory.
///
const string DFDCTruth_factory::toString(void)
{
	if (_data.size() <= 0)
		return "";
	
	stringstream s;
	
	s << _data[0]->header() << endl;
	// Simply stream the string representation of each point into "s"	
	for (unsigned int i=0; i < _data.size(); i++)
		s << _data[i]->toString() << endl;
	
	return s.str();
}
