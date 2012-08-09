// $Id$
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

		// Add new module types to end of list *before*  N_MODULE_TYPES
		// To add a new module, one also needs to add a new case to the
		// the switch in the GetModule(type_id_t) method below.
		enum type_id_t{
			
			UNKNOWN,
			
			F250ADC,
			F125ADC,
			F1TDC,
			JLAB_TS,
			JLAB_TI,
			
			N_MODULE_TYPES
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
		DModuleType GetModule(type_id_t id){
			switch(id){
				case F250ADC: return DModuleType(F250ADC, "F250ADC", "JLab Flash 250 MHz ADC");
				case F125ADC: return DModuleType(F125ADC, "F125ADC", "JLab Flash 125 MHz ADC");
				case F1TDC:   return DModuleType(F1TDC  , "F1TDC"  , "JLab F1 TDC");
				case JLAB_TS: return DModuleType(JLAB_TS, "JLAB_TS", "JLab Trigger Supervisor");
				case JLAB_TI: return DModuleType(JLAB_TI, "JLAB_TI", "JLab Trigger Interface");
					
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
		DModuleType GetModule(std::string name){
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
		void GetModuleList(std::vector<DModuleType> &modules){
			
			for(type_id_t id=UNKNOWN; id<N_MODULE_TYPES; id++){
				modules.push_back(GetModule(id));
			}
		}
		
		
	protected:
		type_id_t type;
		std::string name;
		std::string description;
	
};

typedef DModuleType::type_id_t MODULE_TYPE;

#endif // _DModuleType_

