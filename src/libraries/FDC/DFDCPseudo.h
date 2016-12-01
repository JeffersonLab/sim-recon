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
#include <vector>
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
		double t_u,t_v; ///< time of the two cathode clusters
		double phi_u,phi_v; ///< rotation angles for cathode planes
		double w,dw; ///< local coordinate of pseudopoint in the direction perpendicular to the wires and its uncertainty
		double w_c; /// < wire position computed from cathode data, assuming the avalanche occurs at the wire
		double s,ds; ///< local coordinate of pseudopoint in the direction along the wire and its uncertainty
		const DFDCWire* wire; ///< DFDCWire for this wire 
		double time; ///< time corresponding to this pseudopoint.
		int status; ///< status word for pseudopoint
		double covxx,covxy,covyy; ///< Covariance terms for (x,y) 
		double dE; ///< 
		double q; ///< anode charge deduced from cathode strips
		int itrack;
		DVector2 xy; ///< rough x,y coordinates in lab coordinate system

		void toStrings(vector<pair<string,string> > &items)const{ 
		  AddString(items,"u","%3.2f",u);
		  AddString(items,"v","%3.2f",v);
		  AddString(items,"t_u","%3.2f",t_u);
		  AddString(items,"t_v","%3.2f",t_v);
		  AddString(items,"phi_u","%3.2f",phi_u);
		  AddString(items,"phi_v","%3.2f",phi_v);
		  AddString(items, "w", "%3.4f", w);
        AddString(items, "w_c", "%3.4f", w_c);
		  AddString(items, "s", "%3.4f", s);
		  AddString(items, "layer", "%d", wire->layer);
		  AddString(items, "wire", "%d", wire->wire);
		  AddString(items, "time", "%3.1f", time);
		  AddString(items, "status", "%d", status);
		  AddString(items, "x", "%.4f", xy.X());
		  AddString(items, "y", "%.4f", xy.Y());
		  AddString(items, "dE", "%3.1f", dE);
		}

      // For alignment purposes the residuals wrt the alignment parameters are needed.
      // These routines calculate the derivatives of the "s" and "w" variables wrt the alignment parameters
      // Currently implemented alignment parameters are dU, dV, dPhiU, dPhiV, dX (Acts in direction w), dY (acts in direction s)

      enum FDCPseudoD {
         dWcddeltaU=0,
         dWcddeltaV,
         dWcddeltaPhiU,
         dWcddeltaPhiV,
         dSddeltaU,
         dSddeltaV,
         dSddeltaPhiU,
         dSddeltaPhiV,
         dWdX,
         dSdX,
      };

      vector<double> GetFDCPseudoAlignmentDerivatives(){
         // Create the storage vector for each of the derivatives
         size_t nDerivatives = 10;
         vector<double> derivatives(nDerivatives);

         // Useful numbers...
         double sinPhiU = sin(phi_u);
         double cosPhiU = cos(phi_u);
         double sinPhiV = sin(phi_v);
         double cosPhiV = cos(phi_v);
         double sinPhiUmPhiV = sin(phi_u-phi_v);
         double sinPhiUmPhiV2 = sinPhiUmPhiV*sinPhiUmPhiV;
         double cosPhiUmPhiV = cos(phi_u-phi_v);

         // Calculate the derivatives
         derivatives[dWcddeltaU] = sinPhiV/sinPhiUmPhiV;
         derivatives[dWcddeltaV] = -sinPhiU/sinPhiUmPhiV;
         derivatives[dWcddeltaPhiU] = (v-u*cosPhiUmPhiV)*sinPhiV/sinPhiUmPhiV2;
         derivatives[dWcddeltaPhiV] = (u-v*cosPhiUmPhiV)*sinPhiU/sinPhiUmPhiV2;

         derivatives[dSddeltaU] = -cosPhiV/sinPhiUmPhiV;
         derivatives[dSddeltaV] = cosPhiU/sinPhiUmPhiV;
         derivatives[dSddeltaPhiU] = -(v-u*cosPhiUmPhiV)*cosPhiV/sinPhiUmPhiV2;
         derivatives[dSddeltaPhiU] = -(u-v*cosPhiUmPhiV)*cosPhiU/sinPhiUmPhiV2;

         derivatives[dWdX]=1.0;
         derivatives[dSdX]=1.0;

         return derivatives;

      }

};

#endif //DFDCPSEUDO_H
