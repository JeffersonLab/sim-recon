//***********************************************************************************
// DFDCTruth.h - data type containing truth values from HDDM
// Author: Craig Bookwalter (craigb at jlab.org)
// Date:	March 2006
//
// TODO: Add provision for wire and layer number associated with truth point, so that
// comparisons between truth and reconstructed points can be easily made.
//***********************************************************************************

#ifndef DFDCTRUTH_H
#define DFDCTRUTH_H

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

#include <string>
#include <sstream>
#include <iostream>

///
/// class DFDCTruth: definition for FDC truth points from HDDM
/// 
class DFDCTruth : public JObject {
	public:
		HDCLASSDEF(DFDCTruth);		/// DANA identifier
		
		float dEdx;
		float dradius;
		bool primary;
		int track;
		
		float x;
		float y;
		float z;
		
		///
		/// DFDCTruth::header():
		/// Return a string that contains a sensible header for DFDCTruth data
		///
		const string header() const {
			stringstream s;
			s.width(7);
			s << "dEdx:" << " ";
			s.width(7);
			s << "dRadius:" << " ";
			s.width(7);
			s << "primary:" << " ";
			s.width(7);
			s << "track #:" << " ";
			s.width(12);
			s << "point:" << " ";
			return s.str();
		}
		
		///
		/// DFDCTruth::toString():
		/// Return a string representation of a single DFDCTruth object.
		/// 
		const string toString() const {
			stringstream s;
			s.precision(3);
			s.width(7);
			s << dEdx << " ";
			s.width(7);	
			s << dradius << " ";
			s.width(7);
			s << primary << " ";
			s.width(7);
			s << track << " ";
			s.width(7);
			s << "(" << x << "," << y << "," << z << ")" << " ";
			return s.str();
		}
};

#endif // DFDCTRUTH_H

