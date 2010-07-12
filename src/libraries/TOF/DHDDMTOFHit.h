// $Id$
//
//    File: DHDDMTOFHit.h
// Created: Mon Oct 17 15:01:51 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//
// changes: Tue Jun 19 17:29:07 EDT 2007 B.Zihlmann
//          put north and south information to the same structure 

#ifndef _DHDDMTOFHit_
#define _DHDDMTOFHit_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DHDDMTOFHit:public JObject{

    public:
		JOBJECT_PUBLIC(DHDDMTOFHit);
	
		int plane;			// plane (0: vertical, 1: horizontal)
		int bar;				// bar number
		int ptype;			// GEANT particle type
		float t_north;		// time of light at end of bar
		float dE_north;	// attenuated energy deposition
		float t_south;		// time of light at end of bar
		float dE_south;	// attenuated energy deposition
		float x;          // hit location in global coordiantes
		float y;
		float z;
		float px;			// particle momentum
		float py;
		float pz;
		float E;				// particle Energy

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "bar", "%d", bar);
			AddString(items, "plane", "%d", plane);
			AddString(items, "t_north", "%1.3f", t_north);
			AddString(items, "dE_north", "%1.3f", dE_north);
			AddString(items, "t_south", "%1.3f", t_south);
			AddString(items, "dE_south", "%1.3f", dE_south);
			AddString(items, "x", "%12.4e", x);
			AddString(items, "y", "%12.4e", y);
			AddString(items, "z", "%12.4e", z);
			AddString(items, "px", "%12.4e", px);
			AddString(items, "py", "%12.4e", py);
			AddString(items, "pz", "%12.4e", pz);
			AddString(items, "E", "%12.4e", E);
			AddString(items, "ptype", "%d", ptype);
		}
};

#endif // _DHDDMTOFHit_

