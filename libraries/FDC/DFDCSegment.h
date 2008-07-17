//***********************************************************************
// DFDCSegment.h : definition for a track segment built from pseudopoints
//***********************************************************************

#ifndef DFDCSEGMENT_H
#define DFDCSEGMENT_H

#include <JANA/JObject.h>
using namespace jana;

#include "DFDCHit.h"
#include "DFDCWire.h"
#include "DFDCPseudo.h"

#include <DMatrix.h>
#include <sstream>

// Residuals and arc length
struct fdc_track_t{
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
		JOBJECT_PUBLIC(DFDCSegment);			/// DANA identifier
		
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

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "xc", "%3.2f", xc);
			AddString(items, "yc", "%3.2f", yc);
			AddString(items, "rc", "%3.2f", rc);
			AddString(items, "Phi1(rad)", "%3.2f", Phi1);
			AddString(items, "Nhits", "%d", hits.size());
		}
};

#endif //DFDCSEGMENT_H
