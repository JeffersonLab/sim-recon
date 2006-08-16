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
/// into memory, this is not used. See DEventSourceHDDM::Extract_DFDCHit().
///
jerror_t DFDCHit_factory::evnt(JEventLoop *eventLoop, int eventnumber)
{
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
	
