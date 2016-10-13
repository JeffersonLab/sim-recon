// $Id$
//
//    File: DBCALDigiHit.h
// Created: Tue Aug  6 09:14:41 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _DBCALDigiHit_
#define _DBCALDigiHit_

#include <BCAL/DBCALGeometry.h>

#include <JANA/JObject.h>
using namespace jana;

class DBCALDigiHit:public JObject{

	/// This class holds a single hit from a BCAL fADC250 module.
	/// The values are in the digitized form coming from the module
	/// and are therefore uncalibrated.

	public:
		JOBJECT_PUBLIC(DBCALDigiHit);

		int module;
		int layer;
		int sector;
		DBCALGeometry::End end;
		uint32_t pulse_integral; ///< identified pulse integral as returned by FPGA algorithm
		uint32_t pulse_peak;     ///< identified pulse height as returned by FPGA algorithm
		uint32_t pulse_time;     ///< identified pulse time as returned by FPGA algorithm
		uint32_t pedestal;       ///< pedestal info used by FPGA (if any)
		uint32_t QF;             ///< Quality Factor from FPGA algorithms
		uint32_t nsamples_integral;    ///< number of samples used in integral 
		uint32_t nsamples_pedestal;    ///< number of samples used in pedestal
		
		uint32_t datasource;           ///<  0=window raw data, 1=old(pre-Fall16) firmware, 2=Df250PulseData,  3=MC
		
		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "module", "%d", module);
			AddString(items, "layer", "%d", layer);
			AddString(items, "sector", "%d", sector);
			AddString(items, "end", "%s", end==0 ? "upstream":"downstream" );
			AddString(items, "pulse_integral", "%d", pulse_integral);
			AddString(items, "pulse_peak", "%d", pulse_peak);
			AddString(items, "pulse_time", "%d", pulse_time);
			AddString(items, "pedestal", "%d", pedestal);
			AddString(items, "QF", "%d", QF);
			AddString(items, "nsamples_integral", "%d", nsamples_integral);
			AddString(items, "nsamples_pedestal", "%d", nsamples_pedestal);
		}

};

#endif // _DBCALDigiHit_

