// $Id$
// $HeadURL$
//
//    File: DDAQConfig.h
// Created: Mon Sep  8 03:26:42 EDT 2014
// Creator: davidl (on Darwin harriet.local 13.3.0 i386)
//

#ifndef _DDAQConfig_
#define _DDAQConfig_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

/// This class is a base class used for classes that hold
/// DAQ module configuration parameters. A subclass for each
/// type of digitization module exists that has the attributes
/// appropriate for that type of module. (See Df250Config,
/// DF1TDCConfig, ...) This class only holds the rocid and
/// slot_mask fields which are common to all configurations.
/// One of the main purposes of this base class is to allow
/// configuration objects for all module types to be stored
/// in a single container used internally by the DAQ plugin.

class DDAQConfig:public jana::JObject{
	public:
		JOBJECT_PUBLIC(DDAQConfig);
		
		DDAQConfig(uint32_t rocid, uint32_t slot_mask):rocid(rocid),slot_mask(slot_mask){}

		uint32_t rocid;      // crate
		uint32_t slot_mask;  // slots
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "rocid"     , "%d", rocid);
			AddString(items, "slot_mask" , "0x%06x", slot_mask);
		}
		
};

#endif // _DDAQConfig_

