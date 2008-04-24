// $Id$
//
//    File: DTOFPoint.h
// Created: Tue Oct 18 09:50:52 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#ifndef _DTOFPoint_
#define _DTOFPoint_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

class DTOFPoint:public JObject{
    public:
        JOBJECT_PUBLIC(DTOFPoint);

        unsigned int trackid;  //index to DHDDMTOFTruth (temporary)
		  int ptype;             // GEANT particle type
        float x, y, z;         //reconstructed position
        float t;               //reconstructed time
        float dedx;            //reconstructed dedx
        unsigned int nhits;    //number of hits in this point
        unsigned int hits[16]; //indices to DTOFHit (temporary form)
        float chisq;           //chisquare of this point

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "trackid", "%d", trackid);
			AddString(items, "ptype", "%d", ptype);
			AddString(items, "x", "%1.3f", x);
			AddString(items, "y", "%1.3f", y);
			AddString(items, "z", "%1.3f", z);
			AddString(items, "t", "%1.3f", t);
			AddString(items, "dedx", "%1.3f", dedx);
			AddString(items, "nhits", "%d", nhits);
			AddString(items, "chisq", "%1.3f", chisq);
		}
};

#endif // _DTOFPoint_

