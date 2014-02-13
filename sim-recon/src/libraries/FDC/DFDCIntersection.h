// $Id$
//
//    File: DFDCIntersection.h
// Created: Tue Oct 30 11:24:53 EDT 2007
// Creator: davidl (on Darwin fwing-dhcp95.jlab.org 8.10.1 i386)
//

#ifndef _DFDCIntersection_
#define _DFDCIntersection_

#include <cmath>

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include <DVector3.h>
#include "FDC/DFDCHit.h"
#include "FDC/DFDCWire.h"

class DFDCIntersection:public JObject{
	public:
		JOBJECT_PUBLIC(DFDCIntersection);
		
		const DFDCHit *hit1;
		const DFDCHit *hit2;
		const DFDCWire *wire1;
		const DFDCWire *wire2;
		DVector3 pos;

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "layer1", "%d", wire1->layer);
			AddString(items, "wire1", "%d", wire1->wire);
			AddString(items, "angle1(deg)", "%3.1f", wire1->angle*180.0/M_PI);
			AddString(items, "layer2", "%d", wire2->layer);
			AddString(items, "wire2", "%d", wire2->wire);
			AddString(items, "angle2(deg)", "%3.1f", wire2->angle*180.0/M_PI);
		}
};

#endif // _DFDCIntersection_

