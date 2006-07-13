//**********************************************************************************
// DFDCHit_factory.h: Class definition for the factory responsible for extracting
// FDC hit data from the data file.
// Author: Craig Bookwalter
// Date: Mar 2006
//**********************************************************************************

#ifndef DFACTORY_DFDCHIT_H
#define DFACTORY_DFDCHIT_H

#include "JANA/JFactory.h"
#include "JANA/JEventLoop.h"
#include "HDDM/hddm_s.h"
#include "DFDCHit.h"
#include "DFDCGeometry.h"

#include <sstream>

///
/// class DFDCHit_factory: defines a JFactory that extracts FDC hit data out of the 
/// data file and fills a vector of DFDCHit objects.
///
class DFDCHit_factory:public JFactory<DFDCHit>{
	public:
		///
		/// DFDCHit_factory::DFDCHit_factory():
		/// default constructor -- empty, for now.
		///
		DFDCHit_factory();
		
		///
		/// DFDCHit_factory::~DFDCHit_factory()
		/// default destructor -- also empty for now.
		///
		~DFDCHit_factory();

		///
		/// DFDCHit_factory::Extract_HDDM():
		/// Reads an event from the data file and distills out the FDC information, creating a new
		/// DFDCHit object for each FDC hit encountered in the data. If you wish to understand the 
		/// s_Blah_t structures, see the documentation for hddm_s.h.
		///
		jerror_t Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v);

		///
		/// DFDCHit_factory::toString(): 
		/// Provides a sensible std::string representation of all of the data in the factory.
		///
		const string toString(void);
	
	protected:
		///
		/// DFDCHit_factory::evnt():
		/// This would be used if this factory was going to process data currently exisiting in
		/// memory; since this factory's role is to read events from the data file and put them
		/// into memory, this is not used. See DFDCHit_factory::Extract_HDDM().
		///
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);
		
	private:
		DFDCGeometry _geo;
};

#endif // DFACTORY_DFDCHIT_H

