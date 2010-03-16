// $Id$
//
//    File: DTOFMCResponse.h
// Created: Mon Aug 15 11:33:45 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8 i686)
//

#ifndef _DTOFMCResponse_
#define _DTOFMCResponse_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DTOFMCResponse:public JObject{
    public:
        JOBJECT_PUBLIC(DTOFMCResponse);
	
        int orientation;  // 0: vertical,  1: horizontal
			int ptype;        // particle type
			int bar;          // bar number
        float y;          // x/y position of bar center
        float t_north;          // time of light at end of bar  (smeared) 
        float E_north;          // attenuated energy deposition  (smeared)
        float t_south;          // time of light at end of bar  (smeared) 
        float E_south;          // attenuated energy deposition  (smeared)
	int ADC_north;
	int ADC_south;
	int TDC_north;
	int TDC_south;
	
		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "orientation", "%d", orientation);
 			AddString(items, "ptype", "%d", ptype);
			AddString(items, "bar", "%d", bar);
			AddString(items, "y", "%2.3f", y);
			AddString(items, "t_north", "%1.3f", t_north);
			AddString(items, "E_north", "%1.3f", E_north);
			AddString(items, "t_south", "%1.3f", t_south);
			AddString(items, "E_south", "%1.3f", E_south);
			AddString(items, "ADC_north", "%d", ADC_north);
			AddString(items, "ADC_south", "%d", ADC_south);
			AddString(items, "TDC_north", "%d", TDC_north);
			AddString(items, "TDC_south", "%d", TDC_south);
		}
};

#endif // _DTOFMCResponse_

