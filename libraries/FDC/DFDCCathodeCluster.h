//********************************************************************* 
// DFDCCathodeCluster.h: defintion for a cathode strip distribution
// Author: Craig Bookwalter (craigb at jlab.org)
// Date: Apr 2006
//*********************************************************************

#ifndef DFDCCATHODECLUSTER_H
#define DFDCCATHODECLUSTER_H

#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

#include "DFDCHit.h"
#include "JANA/JObject.h"

#define HIT_TIME_DIFF_MIN 10.0

class DFDCCathodeCluster : public JObject {
	public:
		vector<const DFDCHit*> members;		/// DFDCHits that make up this cluster
		HDCLASSDEF(DFDCCathodeCluster);		/// DANA identifier
		int gPlane;							/// #1-74, which plane out of all FDC modules
		int plane;							/// V=1, U=3
		int gLayer;							/// #1-24, which detection layer of all FDC
		int	beginStrip;						/// Strip # of first cluster member
		int endStrip;						/// Strip # of last cluster member
		int maxStrip;						/// Strip # number with highest dE
		int width;							/// beginStrip - endStrip
		float q_tot;						/// total energy/charge deposited in the cluster
		
		/// Return a sensible string representation of this object
		const string toString() const {
			stringstream s;
			s.precision(4);
			s.width(7);
			s << gPlane << " "; 
			s.width(7);
			s << gLayer << " ";
			s.width(7); 
			s << beginStrip << " ";
			s.width(7);
			s << endStrip << " ";
			s.width(7);
			s << maxStrip << " ";
			s.width(7);
			s << width << " ";
			s.width(7);
			s << q_tot;
			return s.str();
		}
		
		/// Print a sensible header for a list of stringified DFDCCathodeClusters.
		const string header() const {
			stringstream s;
			s.width(7);
			s << "gPla:" << " ";	// gLayer
			s.width(7);
			s << "gLay:" << " ";	// gPlane
			s.width(7);
			s << "bStrp:" << " ";	// beginStrip
			s.width(7);
			s << "eStrp:" << " ";	// endStrip
			s.width(7);
			s << "mxStrp:" << " "; 	// maxStrip
			s.width(7);
			s << "width:" << " ";	// duh
			s.width(7);
			s << "q_tot:" << " ";	// duh
			
			return s.str();
		}	
};

#endif //DFDCCATHODECLUSTER_H
