//******************************************************************
// DFDCPseudo_factory.h: class definition for a factory creating
// pseudopoints from anode hits and cathode clusters.
// Author: Craig Bookwalter
// Date: Apr 2006
//******************************************************************
#ifndef DFACTORY_DFDCPSEUDO_H
#define DFACTORY_DFDCPSEUDO_H

#include "JANA/JFactory.h"
#include "JANA/JException.h"
#include "JANA/JStreamLog.h"

#include "DFDCPseudo.h"
#include "DFDCCathodeCluster.h"
#include "DFDCHit.h"
#include "DFDCGeometry.h"

#include <algorithm>
#include <map>
#include <cmath>

///
/// class DFDCPseudo_factory: definition for a JFactory that
/// produces pseudopoints from anode hits and DFDCCathodeClusters.
/// For now, it is purely geometry-based.
/// 
class DFDCPseudo_factory : public JFactory<DFDCPseudo> {
	public:
		
		///
		/// DFDCPseudo_factory::DFDCPseudo_factory():
		/// default constructor -- initializes log file
		///
		DFDCPseudo_factory();
		
		///
		/// DFDCPseudo_factory::~DFDCPseudo_factory():
		/// default destructor -- closes log file
		///
		~DFDCPseudo_factory();	
							
	protected:
		///
		/// DFDCPseudo_factory::evnt():
		/// this is the place that anode hits and DFDCCathodeClusters 
		/// are organized into pseudopoints.
		/// For now, this is done purely by geometry, with no drift
		/// information. See also
		/// DFDCPseudo_factory::makePseudo().
		///
		jerror_t evnt(JEventLoop *eventLoop, int eventNo);

		/// 
		/// DFDCPseudo_factory::makePseudo():
		/// performs UV+X matching to create pseudopoints
		///
		void makePseudo( vector<const DFDCHit*>& x,
				 vector<const DFDCCathodeCluster*>& u,
				 vector<const DFDCCathodeCluster*>& v,
				 float angle,
				 int layer);
		///
		/// DFDCPseudo_factory::FindCentroid()
		/// Calculates the centroids of groups of three adjacent strips
		/// containing a peak.
		///
		jerror_t FindCentroid(const vector<const DFDCHit*>& H, 
				   vector<const DFDCHit *>::const_iterator peak,
				      vector<centroid_t> &centroids);

	private:
				
		std::vector<centroid_t>upeaks;
		std::vector<centroid_t>vpeaks;
		DFDCGeometry _geo;
		JStreamLog* _log;
		ofstream* logFile;
};

#endif // DFACTORY_DFDCPSEUDO_H

