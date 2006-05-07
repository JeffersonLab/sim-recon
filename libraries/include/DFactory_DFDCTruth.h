//***************************************************************************** 
// DFactory_DFDCTruth.h: DFactory for extracting truth points for the FDC from
// the HDDM file
// Author: Craig Bookwalter (craigb at jlab.org)
// Date: Mar 2006
//*****************************************************************************

#ifndef DFACTORY_DFDCTRUTH_H
#define DFACTORY_DFDCTRUTH_H

#include "DFactory.h"
#include "DEventLoop.h"
#include "DFDCTruth.h"
#include "DException.h"

#include <sstream>

///
/// class DFactory_DFDCTruth: definition for a DFactory that
/// extracts FDC truth points from an HDDM data file.
///
class DFactory_DFDCTruth : public DFactory<DFDCTruth>{
	public:
		///
		/// DFactory_DFDCTruth::DFactory_DFDCTruth():
		/// default constructor -- empty for now
		///
		DFactory_DFDCTruth();

		///
		/// DFactory_DFDCTruth::~DFactory_DFDCTruth():
		/// default destructor -- empty for now as well
		///
		~DFactory_DFDCTruth();

		///
		/// DFactory_DFDCTruth::Extract_HDDM():
		/// reads an event from a data file and converts FDC truth information into 
		/// DFDCTruth objects. For more detail on the s_Blah_t structures, see
		/// the documentation for hddm_s.h.
		///
		derror_t Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v);

		///
		/// DFactory_DFDCTruth::toString():
		/// Print a sensible std::string representation of the contents of this factory.
		///
		const string toString(void);
	
	protected:
		///
		/// DFactory_DFDCTruth::evnt():
		/// this function would be defined if this factory was processing data
		/// already in memory; since this factory reads an HDDM file and puts
		/// truth data into memory, the DFactory_DFDCTruth::Extract_HDDM() method
		/// is used.
		///
		derror_t evnt(DEventLoop *eventLoop, int eventnumber);
};

#endif // DFACTORY_DFDCTRUTH_H

