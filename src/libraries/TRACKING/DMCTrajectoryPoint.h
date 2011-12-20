// $Id$
//
//    File: DMCTrajectoryPoint.h
// Created: Mon Jun 12 09:29:58 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.6.0 powerpc)
//

#ifndef _DMCTrajectoryPoint_
#define _DMCTrajectoryPoint_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DMCTrajectoryPoint:public jana::JObject{
	public:
		JOBJECT_PUBLIC(DMCTrajectoryPoint);
		
		float x,y,z,t;
		float px,py,pz;
		float E, dE;
		int primary_track;
		int track;
		int part;
		float radlen;
		float step;
		int mech;

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "x", "%1.3f", x);
			AddString(items, "y", "%1.3f", y);
			AddString(items, "z", "%1.3f", z);
			AddString(items, "t", "%1.3f", t/1.0E-9);
			AddString(items, "px", "%1.3f", px);
			AddString(items, "py", "%1.3f", py);
			AddString(items, "pz", "%1.3f", pz);
			AddString(items, "E", "%1.3f", E);
			AddString(items, "dE(MeV)", "%1.3f", 1000.0*dE);
			AddString(items, "primary", "%d", primary_track);
			AddString(items, "track", "%d", track);
			AddString(items, "part", "%d", part);
			AddString(items, "radlen", "%1.3f", radlen);
			AddString(items, "step", "%1.3f", step);
			AddString(items, "mech", "%d", mech);
		}
};

#endif // _DMCTrajectoryPoint_

