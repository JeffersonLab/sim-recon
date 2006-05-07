//**********************************************************************************
// DFactory_DFDCHit.h: Class definition for the factory responsible for extracting
// FDC hit data from the data file.
// Author: Craig Bookwalter
// Date: Mar 2006
//**********************************************************************************

#ifndef DFACTORY_DFDCHIT_H
#define DFACTORY_DFDCHIT_H

#include "DFactory.h"
#include "DEventLoop.h"
#include "DFDCHit.h"
#include "DFDCGeometry.h"

#include <sstream>

///
/// class DFactory_DFDCHit: defines a DFactory that extracts FDC hit data out of the 
/// data file and fills a vector of DFDCHit objects.
///
class DFactory_DFDCHit:public DFactory<DFDCHit>{
	public:
		///
		/// DFactory_DFDCHit::DFactory_DFDCHit():
		/// default constructor -- empty, for now.
		///
		DFactory_DFDCHit();
		
		///
		/// DFactory_DFDCHit::~DFactory_DFDCHit()
		/// default destructor -- also empty for now.
		///
		~DFactory_DFDCHit();

		///
		/// DFactory_DFDCHit::Extract_HDDM():
		/// Reads an event from the data file and distills out the FDC information, creating a new
		/// DFDCHit object for each FDC hit encountered in the data. If you wish to understand the 
		/// s_Blah_t structures, see the documentation for hddm_s.h.
		///
		derror_t Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v);

		///
		/// DFactory_DFDCHit::toString(): 
		/// Provides a sensible std::string representation of all of the data in the factory.
		///
		const string toString(void);
	
	protected:
		///
		/// DFactory_DFDCHit::evnt():
		/// This would be used if this factory was going to process data currently exisiting in
		/// memory; since this factory's role is to read events from the data file and put them
		/// into memory, this is not used. See DFactory_DFDCHit::Extract_HDDM().
		///
		derror_t evnt(DEventLoop *eventLoop, int eventnumber);
		
	private:
		DFDCGeometry _geo;
};

#endif // DFACTORY_DFDCHIT_H

