//***************************************************************************** 
// DFDCTruth_factory.h: JFactory for extracting truth points for the FDC from
// the HDDM file
// Author: Craig Bookwalter (craigb at jlab.org)
// Date: Mar 2006
//*****************************************************************************

#ifndef DFACTORY_DFDCTRUTH_H
#define DFACTORY_DFDCTRUTH_H

#include "JANA/JFactory.h"
#include "JANA/JEventLoop.h"
#include "JANA/JException.h"

#include "HDDM/hddm_s.h"
#include "DFDCTruth.h"

#include <sstream>

///
/// class DFDCTruth_factory: definition for a JFactory that
/// extracts FDC truth points from an HDDM data file.
///
class DFDCTruth_factory : public JFactory<DFDCTruth>{
	public:
		///
		/// DFDCTruth_factory::DFDCTruth_factory():
		/// default constructor -- empty for now
		///
		DFDCTruth_factory();

		///
		/// DFDCTruth_factory::~DFDCTruth_factory():
		/// default destructor -- empty for now as well
		///
		~DFDCTruth_factory();

		///
		/// DFDCTruth_factory::Extract_HDDM():
		/// reads an event from a data file and converts FDC truth information into 
		/// DFDCTruth objects. For more detail on the s_Blah_t structures, see
		/// the documentation for hddm_s.h.
		///
		jerror_t Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v);

		///
		/// DFDCTruth_factory::toString():
		/// Print a sensible std::string representation of the contents of this factory.
		///
		const string toString(void);
	
	protected:
		///
		/// DFDCTruth_factory::evnt():
		/// this function would be defined if this factory was processing data
		/// already in memory; since this factory reads an HDDM file and puts
		/// truth data into memory, the DFDCTruth_factory::Extract_HDDM() method
		/// is used.
		///
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);
};

#endif // DFACTORY_DFDCTRUTH_H

