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

#include "DFDCWire.h"
#include <sstream>
#include <DVector2.h>


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
	

		double u,v; ///< centroid positions in the two cathode views
		double w,dw; ///< local coordinate of pseudopoint in the direction perpendicular to the wires and its uncertainty
		double s,ds; ///< local coordinate of pseudopoint in the direction along the wire and its uncertainty
		const DFDCWire* wire; ///< DFDCWire for this wire 
		double time; ///< time corresponding to this pseudopoint.
		int status; ///< status word for pseudopoint
		double covxx,covxy,covyy; ///< Covariance terms for (x,y) 
		double dE; ///< 
		int itrack;
		DVector2 xy; ///< rough x,y coordinates in lab coordinate system

		void toStrings(vector<pair<string,string> > &items)const{ 
		  AddString(items,"u","%3.2f",u);
		  AddString(items,"v","%3.2f",v);
		  AddString(items, "w", "%3.2f", w);
		  AddString(items, "s", "%3.2f", s);
		  AddString(items, "layer", "%d", wire->layer);
		  AddString(items, "wire", "%d", wire->wire);
		  AddString(items, "time", "%3.1f", time);
		  AddString(items, "status", "%d", status);
		  AddString(items, "x", "%3.1f", xy.X());
		  AddString(items, "y", "%3.1f", xy.Y());
		  AddString(items, "dE", "%3.1f", dE);
		}

};

#endif //DFDCPSEUDO_H
