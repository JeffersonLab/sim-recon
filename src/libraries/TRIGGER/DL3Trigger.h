// $Id$
//
//    File: DL3Trigger.h
// Created: Wed Jul 31 14:34:24 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _DL3Trigger_
#define _DL3Trigger_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DL3Trigger:public jana::JObject{
	public:
		JOBJECT_PUBLIC(DL3Trigger);
		
		/// The DL3Trigger object is used to tell the level-3 trigger process whether
		/// or not to discard the event. The intent is for the algorithm to be contained
		/// in the DL3Trigger_factory class which creates a single DL3Trigger object
		/// for each event. 
		///
		/// The factory should set the 3 members of the DL3Trigger object based on its
		/// results. The L3_decision value must be set to one of the values in the
		/// L3_decision_t enum such as DL3Trigger::kDiscardEvent. An enum is used in
		/// case additional states are needed later. The "status" member is a 64bit
		/// word that will be written out with those events that are kept. The intent
		/// is for this to hold a set of flags indicating the cause for the decision.
		/// The meaning of the bits in therefore algorithm specific. Because of this,
		/// a 3rd member "algorithm" is used to record which L3 algorithm (i.e. which
		/// version) was used for this event. Recording this for each event is important
		/// in case a set of filtered events from multiple runs are ever concatentated
		/// into a single file.
		///
		/// Different L3 algorithms can be developed, each in its own plugin. (The
		/// DL3Trigger factory in the plugin takes precedence over one statically linked
		/// into the executable.) To implement a L3 trigger algorithm, one should create
		/// a plugin that contains a DL3Trigger factory object. You can do this using the
		/// mkplugin_factory script that comes with JANA e.g.:
		///
		/// mkfactory_plugin DL3Trigger
		///
		/// This will create directory named "DL3Trigger" with skeleton files (including
		/// a Makefile) that can be built into an appropriate plugin. Please note that
		/// the generated code will need to be modified to use the existing DL3Trigger
		/// class definitiong in DANA. Specifically:
		///
		///  1.) remove the skeleton DL3Trigger.h file
		///
		///  2.) Inside DL3Trigger_factory.h, replace
		///
		///           #include "DL3Trigger.h"
		///
		///               with
		///
		///           #include <TRIGGER/DL3Trigger.h>
		///
		/// The plugin name is taken as the name of the directory (in this case "DL3Trigger").
		/// To change the name of the plugin, simply rename the DL3Trigger directory.
		///
		///  Once the plugin is built use it with any DANA program in the standard way:
		///
		///   e.g.
		///
		///         hd_ana -PPLUGINS=DL3Trigger ....
		///
		
		enum L3_decision_t{
			kNO_DECISION,
			kKEEP_EVENT,
			kDISCARD_EVENT
		};
		
		DL3Trigger(L3_decision_t L3_decision=kNO_DECISION, uint64_t status=0L, uint32_t algorithm=0)
			:L3_decision(L3_decision),status(status),algorithm(algorithm),mva_response(-1E6){}
		
		L3_decision_t L3_decision;  // keep event or not?
		uint64_t status;            // algorithm specific status bits
		uint32_t algorithm;         // unique identifier for this algorithm
		double mva_response;        // response of MVA algorithm
		
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "L3_decision", "%d", L3_decision);
			AddString(items, "status", "0x%16x", status);
			AddString(items, "algorithm", "0x%8x", algorithm);
		}
		
};

#endif // _DL3Trigger_

