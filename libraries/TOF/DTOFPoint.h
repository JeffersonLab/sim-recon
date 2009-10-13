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
#include <DVector3.h>
using namespace jana;

class DTOFPoint:public JObject{
    public:
        JOBJECT_PUBLIC(DTOFPoint);

	DVector3 pos;   	//reconstructed position
        float t;               //reconstructed time
        float dedx;            //reconstructed dedx
        float chisq;           //chisquare of this point

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "x", "%1.3f", pos.x());
			AddString(items, "y", "%1.3f", pos.y());
			AddString(items, "z", "%1.3f", pos.z());
			AddString(items, "t", "%1.3f", t);
			AddString(items, "dedx", "%1.3f", dedx);
			AddString(items, "chisq", "%1.3f", chisq);
		}
};

#endif // _DTOFPoint_

