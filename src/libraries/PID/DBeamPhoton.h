// $Id$
//
//    File: DBeamPhoton.h
// Created: Thu Feb 14 10:11:52 EST 2008
// Creator: davidl (on Darwin fwing-dhcp13.jlab.org 8.10.1 i386)
//

#ifndef _DBeamPhoton_
#define _DBeamPhoton_

#include <PID/DKinematicData.h>

class DBeamPhoton: public DKinematicData{
	public:
		JOBJECT_PUBLIC(DBeamPhoton);
		
		double t; ///< Time at which photon arrives at interaction vertex location
		
		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "E(GeV)", "%3.3f", momentum().Mag());
			AddString(items, "theta(deg)", "%3.3f", momentum().Theta()*180.0/M_PI);
			AddString(items, "phi(deg)", "%3.1f", momentum().Phi()*180.0/M_PI);
			AddString(items, "t(ns)", "%3.1f", t);
		}
};



#endif // _DBeamPhoton_

