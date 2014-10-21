// $Id$
//
//    File: DTOFTruth.h
// Created: Mon Oct 17 13:58:02 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#ifndef _DTOFTruth_
#define _DTOFTruth_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

class DTOFTruth:public JObject{

    public:
        JOBJECT_PUBLIC(DTOFTruth);

        int track;         //  track index
        int itrack;        // MCThrown track index
        int primary;       //  0: secondary, 1: primary
        float x, y, z;     //  true point of intersection
	float px,py,pz;    //  momentum of the particle
        float t;           //  true time
        float E;           //  energy of the particle
	int ptype;         //  GEANT particle type

	void toStrings(vector<pair<string,string> > &items)const{
	   AddString(items, "track", "%d", track);
	   AddString(items, "itrack", "%d", itrack);
	   AddString(items, "primary", "%d", primary);
	   AddString(items, "ptype", "%d", ptype);
	   AddString(items, "x", "%1.3f", x);
	   AddString(items, "y", "%1.3f", y);
	   AddString(items, "z", "%1.3f", z);
	   AddString(items, "t", "%1.3f", t);
	   AddString(items, "px", "%1.3f", px);
	   AddString(items, "py", "%1.3f", py);
	   AddString(items, "pz", "%1.3f", pz);
	   AddString(items, "E", "%1.3f", E);
	}
};

#endif // _DTOFTruth_

