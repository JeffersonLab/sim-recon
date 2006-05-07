//********************************************************************
// DFactory_DFDCCathodeCluster.h: Class definition for a factory that
// associates cathode hits into clusters.
// Author: Craig Bookwalter
//********************************************************************

#ifndef DFACTORY_DFDCCATHODECLUSTER_H
#define DFACTORY_DFDCCATHODECLUSTER_H

#include "DFactory.h"
#include "DFDCCathodeCluster.h"
#include "DFDCHit.h"
#include "DException.h"
#include "DFDCGeometry.h"
#include "DStreamLog.h"

#include <algorithm>
#include <map>
#include <cmath>

///
/// class DFactory_DFDCCathodeCluster: 
/// defines a DFactory for producing groups of cathode strips that form a cluster
///  
class DFactory_DFDCCathodeCluster : public DFactory<DFDCCathodeCluster> {
	public:
		///
		/// DFactory_DFDCCathodeCluster::DFactory_DFDCCathodeCluster():
		///	default constructor--initializes log file
		///
		DFactory_DFDCCathodeCluster();
		
		///
		/// DFactory_DFDCCathodeCluster::~DFactory_DFDCCathodeCluster():
		/// default destructor--closes log file.
		///
		~DFactory_DFDCCathodeCluster();
		
		///
		/// DFactory_DFDCCathodeCluster::pique():
		/// takes a single layer's worth of cathode hits and attempts to 
		/// create DFDCCathodeClusters
		/// by grouping together hits with consecutive strip numbers.
		///
		void pique(vector<const DFDCHit*>& h);
		
		///
		/// DFactory_DFDCCathodeCluster::toString():
		/// returns a sensible std::string representation of the data contained in this 
		/// factory.
		///		
		const string toString();
	
	protected:
		///
		/// DFactory_DFDCCathodeCluster::evnt():
		/// This (along with DFactory_DFDCCathodeCluster::pique()) 
		/// is the place cathode hits are associated into cathode clusters. This function 
		/// should eventually be modified to do more sophisticated peak finding. 
		///
		derror_t evnt(DEventLoop *eventLoop, int eventNo);	
		
	private:
		DStreamLog* _log;
		ofstream* _logFile;
};

#endif // DFACTORY_DFDCCATHODECLUSTER_H

