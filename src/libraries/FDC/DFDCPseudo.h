//***********************************************************************
// DFDCPseudo.h : definition for a set of FDCHits that have gone 
// through first-order reconstruction. 
// Author: Craig Bookwalter (craigb at jlab.org)
// Date:	March 2006
//***********************************************************************

#ifndef DFDCPSEUDO_H
#define DFDCPSEUDO_H

#include <JANA/JObject.h>
using namespace jana;

#include "DFDCHit.h"
#include "DFDCWire.h"
#include <DMatrix.h>

#include <sstream>

typedef struct {
  float pos;
  float q;
  int numstrips;
  float t; // mean time of strips in peak
  float t_rms; // rms of strips in peak
}centroid_t;

///
/// class DFDCPseudo: definition for a reconstructed point in the FDC
/// 
class DFDCPseudo : public JObject {
	public :
		JOBJECT_PUBLIC(DFDCPseudo);			/// DANA identifier
		
		/// 
		/// DFDCPseudo::DFDCPseudo():
		/// Default constructor-- provide the X, Y, global layer #, and resolution
		///
		DFDCPseudo(){}
	

		float w,dw; //local coordinate of pseudopoint in the direction 
		            //perpendicular to the wires and its uncertainty
		float s,ds; //local coordinate of pseudopoint in the direction 
		            // along the wire and its uncertainty
		const DFDCWire* wire; ///< DFDCWire for this wire 
		float time; // time corresponding to this pseudopoint.
		float dist;	// drift distance from time
		int status; // status word for pseudopoint
		float x,y;  // coordinates rotated into lab coordinate system
		float covxx,covxy,covyy;
		float dE;

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "w", "%3.2f", w);
			AddString(items, "s", "%3.2f", s);
			AddString(items, "layer", "%d", wire->layer);
			AddString(items, "wire", "%d", wire->wire);
			AddString(items, "time", "%3.1f", time);
			AddString(items, "dist", "%d", dist);
			AddString(items, "status", "%d", status);
			AddString(items, "x", "%3.1f", x);
			AddString(items, "y", "%3.1f", y);
			AddString(items, "dE", "%3.1f", dE);
		}

};

#endif //DFDCPSEUDO_H
