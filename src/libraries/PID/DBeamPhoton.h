// $Id$
//
//    File: DBeamPhoton.h
// Created: Thu Feb 14 10:11:52 EST 2008
// Creator: davidl (on Darwin fwing-dhcp13.jlab.org 8.10.1 i386)
//

#ifndef _DBeamPhoton_
#define _DBeamPhoton_

#include <PID/DKinematicData.h>

class DBeamPhoton: public DKinematicData
{
	public:
		JOBJECT_PUBLIC(DBeamPhoton);
		
		unsigned int dCounter = 0;
		DetectorSystem_t dSystem = SYS_NULL; //SYS_TAGM or SYS_TAGH (SYS_NULL if tag MCGEN and misses the tagger!)

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "E(GeV)", "%f", momentum().Mag());
			AddString(items, "System", "%s", SystemName(dSystem));
			AddString(items, "Counter", "%d", dCounter);
			AddString(items, "t(ns)", "%3.1f", time());
		}
};



#endif // _DBeamPhoton_

