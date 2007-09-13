//***********************************************************************
// DFDCSegment.h : definition for a track segment built from pseudopoints
//***********************************************************************

#ifndef DFDCSEGMENT_H
#define DFDCSEGMENT_H

#include "DFDCHit.h"
#include "DFDCWire.h"
#include "DFDCPseudo.h"
#include "JANA/JObject.h"
#include <DMatrix.h>
#include <sstream>

// Residuals and arc length
typedef struct fdc_track_t{
  unsigned int hit_id;
  double sign; // Sign of left-right ambiguity resolution for this hit
  double dx,dy; //residuals
  double chi2; // chi2 contribution of this point
  double s; // path length
};


///
/// class DFDCSegment: definition for a track segment in the FDC
/// 
class DFDCSegment : public JObject {
	public :
		HDCLASSDEF(DFDCSegment);			/// DANA identifier
		
		/// 
		/// DFDCSegment::DFDCSegment():
		/// Default constructor
		///
		DFDCSegment(){}

	        DMatrix S;   // state vector
		DMatrix cov; // ... and its covariance matrix
		double chisq;

		// circle parameters
		double xc,yc,rc;  
		// azimuthal angle of the intersection of the helical path to 
		// the most downstream plane in a given package containing a
		// hit
		double Phi1;		                

		// List of pseudopoints belonging to this track segment
		vector<DFDCPseudo *>hits;	
		// Supplementary track info
		vector<fdc_track_t>track;
};

#endif //DFDCSEGMENT_H
