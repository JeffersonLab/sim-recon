#ifndef _DTPOLRingDigiHit_
#define _DTPOLRingDigiHit_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DTPOLRingDigiHit:public jana::JObject{
	public:
		JOBJECT_PUBLIC(DTPOLRingDigiHit);
		

		int      ring;                 // ring number 1-24
		uint32_t pulse_integral;       // identified pulse integral as returned by FPGA algorithm
		uint32_t pulse_time;           // identified pulse time as returned by FPGA algorithm
		uint32_t pedestal;             // pedestal info used by FPGA (if any)
		uint32_t QF;                   // Quality Factor from FPGA algorithms
		uint32_t nsamples_integral;    // number of samples used in integral 
		uint32_t nsamples_pedestal;    // number of samples used in pedestal
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "ring", "%d", ring);
			AddString(items, "pulse_integral", "%d", pulse_integral);
			AddString(items, "pulse_time", "%d", pulse_time);
			AddString(items, "pedestal", "%d", pedestal);
			AddString(items, "QF", "%d", QF);
			AddString(items, "nsamples_integral", "%d", nsamples_integral);
			AddString(items, "nsamples_pedestal", "%d", nsamples_pedestal);
		}
		
};

#endif // _DTPOLRingDigiHit_
