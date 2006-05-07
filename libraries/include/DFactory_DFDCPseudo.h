//******************************************************************
// DFactory_DFDCPseudo.h: class definition for a factory creating
// pseudopoints from anode hits and cathode clusters.
// Author: Craig Bookwalter
// Date: Apr 2006
//******************************************************************
#ifndef DFACTORY_DFDCPSEUDO_H
#define DFACTORY_DFDCPSEUDO_H

#include "DFactory.h"
#include "DFDCPseudo.h"
#include "DFDCCathodeCluster.h"
#include "DFDCHit.h"
#include "DException.h"
#include "DFDCGeometry.h"
#include "DStreamLog.h"

#include <algorithm>
#include <map>
#include <cmath>

///
/// class DFactory_DFDCPseudo: definition for a DFactory that
/// produces pseudopoints from anode hits and DFDCCathodeClusters.
/// For now, it is purely geometry-based.
/// 
class DFactory_DFDCPseudo : public DFactory<DFDCPseudo> {
	public:
		
		///
		/// DFactory_DFDCPseudo::DFactory_DFDCPseudo():
		/// default constructor -- initializes log file
		///
		DFactory_DFDCPseudo();
		
		///
		/// DFactory_DFDCPseudo::~DFactory_DFDCPseudo():
		/// default destructor -- closes log file
		///
		~DFactory_DFDCPseudo();	
							
	protected:
		///
		/// DFactory_DFDCPseudo::evnt():
		/// this is the place that anode hits and DFDCCathodeClusters are organized into pseudopoints.
		/// For now, this is done purely by geometry, with no drift or peak-finding. See also
		/// DFactory_DFDCPseudo::makePseudo().
		///
		derror_t evnt(DEventLoop *eventLoop, int eventNo);

		/// 
		/// DFactory_DFDCPseudo::makePseudo():
		/// performs UV+X matching to create pseudopoints
		///
		void makePseudo(	map<const int, const DFDCHit*>& x,
							vector<const DFDCCathodeCluster*>& u,
							vector<const DFDCCathodeCluster*>& v,
							float angle,
							int layer);

		///
		/// DFactory_DFDCPseudo::intersectX():
		/// finds the X coordinate of a U-V intersection
		///
		float intersectX(int u, int v);
		
		///
		/// DFactory_DFDCPseudo::intersectY():
		/// finds the Y coordinate of a U-V intersection
		///
		float intersectY(int u, int v);

	private:
		DFDCGeometry _geo;
		DStreamLog* _log;
		ofstream* logFile;
};

#endif // DFACTORY_DFDCPSEUDO_H

