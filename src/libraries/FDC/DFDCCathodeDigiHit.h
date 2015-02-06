// $Id$
//
//    File: DFDCCathodeDigiHit.h
// Created: Wed Aug  7 11:53:57 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _DFDCCathodeDigiHit_
#define _DFDCCathodeDigiHit_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DFDCCathodeDigiHit:public jana::JObject{
	public:
		JOBJECT_PUBLIC(DFDCCathodeDigiHit);

		uint32_t package;
		uint32_t chamber;
		uint32_t view;
		uint32_t strip;
		uint32_t strip_type;		
		uint32_t pulse_integral;       ///< identified pulse integral as returned by FPGA algorithm
		uint32_t pulse_time;           ///< identified pulse time as returned by FPGA algorithm
		uint32_t pedestal;             ///< pedestal info used by FPGA (if any)
		uint32_t QF;                   ///< Quality Factor from FPGA algorithms
		uint32_t nsamples_integral;    ///< number of samples used in integral 
		uint32_t nsamples_pedestal;    ///< number of samples used in pedestal
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "package", "%d", package);
			AddString(items, "chamber", "%d", chamber);
			AddString(items, "view", "%d", view);
			AddString(items, "strip", "%d", strip);
			AddString(items, "strip_type", "%d", strip_type);
			AddString(items, "pulse_integral", "%d", pulse_integral);
			AddString(items, "pulse_time", "%d", pulse_time);
			AddString(items, "pedestal", "%d", pedestal);
			AddString(items, "QF", "%d", QF);
			AddString(items, "nsamples_integral", "%d", nsamples_integral);
			AddString(items, "nsamples_pedestal", "%d", nsamples_pedestal);
		}
		
};

#endif // _DFDCCathodeDigiHit_

