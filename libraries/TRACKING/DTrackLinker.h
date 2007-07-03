//***********************************************************************
// DTrackLinker.h : header file for forward linking factory  
//***********************************************************************

#ifndef DTRACKLINKER_H
#define DTRACKLINKER_H

#include "FDC/DFDCPseudo.h"
#include "JANA/JObject.h"
#include <DMatrix.h>
#include <sstream>


///
/// class DTrackLinker: definition for a linked set of segments
/// 
class DTrackLinker : public JObject {
	public :
		HDCLASSDEF(DTrackLinker);			/// DANA identifier
		
		/// 
		/// DTrackLinker::DTrackLinker():
		/// Default constructor
		///
		DTrackLinker(){}

	        DMatrix S;   // state vector
		DMatrix cov; // ... and its covariance matrix
		double chisq;

		// List of points belonging to this track
		vector<DFDCPseudo *>points;	

};

#endif //DTRACKLINKER_H
