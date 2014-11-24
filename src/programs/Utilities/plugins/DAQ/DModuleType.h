// $Id$
// $HeadURL$
//
//    File: DModuleType.h
// Created: Thu Aug  9 08:56:08 EDT 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//

#ifndef _DModuleType_
#define _DModuleType_

#include <string>

class DModuleType{
	public:

		// NOTE: The following enum MUST be kept in alignment with the DAQ
		// group's definitions! The real data will include this module
		// type in the Block Header for JLab modules (bits 18-21).
		//
		// These were taken from mc2coda.h v2
		// (see /site/coda/3.0/extensions/mc2coda/v2.0/mc2coda.h)
		//
		// Add new module types to end of list *before*  N_MODULE_TYPES
		// To add a new module, one also needs to add a new case to the
		// the switch in the GetModule(type_id_t) method below.
		enum type_id_t {
			TID,               // =0
			FADC250,           // =1
			FADC125,           // =2
			F1TDC32,           // =3
			F1TDC48,           // =4
			JLAB_TS,           // =5
			TD,                // =6
			SSP,               // =7
			JLAB_DISC,         // =8
			MODULE_TYPE_RES1,  // =9
			MODULE_TYPE_RES2,  // =10
			MODULE_TYPE_RES3,  // =11
			MODULE_TYPE_RES4,  // =12
			MODULE_TYPE_RES5,  // =13
			MODULE_TYPE_RES6,  // =14
			MODULE_TYPE_RES7,  // =15

			// The following are not defined by the DAQ group (i.e.
			// they don't control the data format so can't encode
			// the module type in it)
			UNKNOWN,           // =16
			VMECPU,            // =17
			CAEN1190,          // =18
			CAEN1290,          // =19

			N_MODULE_TYPES     // Make sure this is the last thing in the enum!
		};

		//------------------------------
		// Constructor
		//
		DModuleType(type_id_t type, std::string name, std::string description):type(type),name(name),description(description){}

	
		//------------------------------
		// GetModule(type_id_t)
		//
		/// Given the type_id of a module type, return a DModuleType object.
		/// This is much more efficient than GetModule(string), but is still
		/// not really intended to be called every event.
		static DModuleType GetModule(type_id_t id){
			switch(id){
				case TID      : return DModuleType(TID,       "TID"      , "JLab Trigger Interface");
				case FADC250  : return DModuleType(FADC250,   "FADC250"  , "JLab Flash 250 MHz ADC");
				case FADC125  : return DModuleType(FADC125,   "FADC125"  , "JLab Flash 125 MHz ADC");
				case F1TDC32  : return DModuleType(F1TDC32,   "F1TDC32"  , "JLab F1 TDC (32 chan)");
				case F1TDC48  : return DModuleType(F1TDC48,   "F1TDC48"  , "JLab F1 TDC (48 chan)");
				case JLAB_TS  : return DModuleType(JLAB_TS,   "JLAB_TS"  , "JLab Trigger Supervisor");
				case TD       : return DModuleType(TD,        "TD"       , "JLab Trigger Distribution");
				case SSP      : return DModuleType(SSP,       "SSP"      , "JLab Sub-System Processor");
				case JLAB_DISC: return DModuleType(JLAB_DISC, "JLAB_DISC", "JLab Discriminator");

				case VMECPU   : return DModuleType(VMECPU,    "VMECPU"   , "VME CPU");
				case CAEN1190 : return DModuleType(CAEN1190,  "CAEN1190" , "CAEN 1190 TDC");
				case CAEN1290 : return DModuleType(CAEN1290,  "CAEN1290" , "CAEN 1290 TDC");
				
				default: return DModuleType(UNKNOWN, "UNKNOWN", "Unknown module type");
			}
		}
	
		//------------------------------
		// GetModule(string)
		//
		/// Given the name of a module type, return a DModuleType object.
		/// This is not a terribly efficient mechanism and is inteded to
		/// be called only at the beginning of event processing rather than
		/// every event.
		static DModuleType GetModule(std::string name){
			std::vector<DModuleType> modules;
			GetModuleList(modules);
			for(unsigned int i=UNKNOWN; i<N_MODULE_TYPES; i++){
				if(modules[i].name == name)return modules[i];
			}
			
			return GetModule(UNKNOWN);
		}	
	
		//------------------------------
		// GetModuleList
		//
		/// Get a list of all module types currently defined. This
		/// will append the full list to the given "modules" container.
		/// (The vector is not cleared on input).
		static void GetModuleList(std::vector<DModuleType> &modules){
			
			for(int id=UNKNOWN; id<N_MODULE_TYPES; id++){
				modules.push_back(GetModule((type_id_t)id));
			}
		}
		
		//------------------------------
		// GetName
		//
		/// Get the name of a module type based on its id. This can
		/// be called without an instance of the class if one already
		/// has the type.
		static string GetName(type_id_t id){
			return GetModule(id).GetName();
		}
	
		//------------------------------
		// GetDescription
		//
		/// Get the name of a module type based on its id. This can
		/// be called without an instance of the class if one already
		/// has the type.
		static string GetDescription(type_id_t id){
			return GetModule(id).GetDescription();
		}
	
		//------------------------------
		// non-static methods
		type_id_t GetType(void) const {return type;}
		string GetName(void) const {return name;}
		string GetDescription(void) const {return description;}
		
	protected:
		type_id_t type;
		std::string name;
		std::string description;
	
};

typedef DModuleType::type_id_t MODULE_TYPE;

#endif // _DModuleType_

