// $Id$
//
//    File: Df125BORConfig.h
// Created: Tue Jan 26 13:04:46 EST 2016
// Creator: davidl (on Darwin harriet.jlab.org 13.4.0 i386)
//

#ifndef _Df125BORConfig_
#define _Df125BORConfig_

#include <JANA/JObject.h>

#include <DAQ/bor_roc.h>

// This class inherits both from JObject and f125config. The former
// so that it can be incorporated easily into the JANA framework.
// The latter so we can use the data struct defined in bor_roc.h.
// The file bor_roc.h exists in 2 places:
//
//  1. in the DAQ library of sim-recon
//  2. in the vme/src/rcm/monitor directory in the online
//


class Df125BORConfig:public jana::JObject, public f125config{
	public:
		JOBJECT_PUBLIC(Df125BORConfig);

		Df125BORConfig(){}
		virtual ~Df125BORConfig(){}
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "rocid"   , "%d", rocid);
			AddString(items, "slot"    , "%d", slot);
			AddString(items, "id"      , "0x%x", board_id);
			AddString(items, "version" , "0x%x", version);
			AddString(items, "proc_version" , "0x%x", proc_version);
			AddString(items, "ctrl1" , "0x%x", ctrl1);
			AddString(items, "proc_blocklevel" , "%d", proc_blocklevel);
		}

};

#endif // _Df125BORConfig_

