// $Id$
//
//    File: DParsedEvent.h
// Created: Mon Mar 28 11:07:41 EDT 2016
// Creator: davidl (on Darwin harriet.jlab.org 13.4.0 i386)
//

#ifndef _DParsedEvent_
#define _DParsedEvent_

#include <string>
#include <map>
using std::string;
using std::map;

#include <JANA/jerror.h>
#include <JANA/JObject.h>
#include <JANA/JEventLoop.h>
using namespace jana;

#include <DANA/DStatusBits.h>
#include <DAQ/daq_param_type.h>
#include <DAQ/DModuleType.h>


#include <DAQ/Df250Config.h>
#include <DAQ/Df250PulseIntegral.h>
#include <DAQ/Df250StreamingRawData.h>
#include <DAQ/Df250WindowSum.h>
#include <DAQ/Df250PulseRawData.h>
#include <DAQ/Df250TriggerTime.h>
#include <DAQ/Df250PulseTime.h>
#include <DAQ/Df250PulsePedestal.h>
#include <DAQ/Df250WindowRawData.h>
#include <DAQ/Df125Config.h>
#include <DAQ/Df125TriggerTime.h>
#include <DAQ/Df125PulseIntegral.h>
#include <DAQ/Df125PulseTime.h>
#include <DAQ/Df125PulsePedestal.h>
#include <DAQ/Df125PulseRawData.h>
#include <DAQ/Df125WindowRawData.h>
#include <DAQ/Df125CDCPulse.h>
#include <DAQ/Df125FDCPulse.h>
#include <DAQ/DF1TDCConfig.h>
#include <DAQ/DF1TDCHit.h>
#include <DAQ/DF1TDCTriggerTime.h>
#include <DAQ/DCAEN1290TDCConfig.h>
#include <DAQ/DCAEN1290TDCHit.h>
#include <DAQ/DCODAEventInfo.h>
#include <DAQ/DCODAROCInfo.h>
#include <DAQ/DTSscalers.h>
#include <DAQ/DEPICSvalue.h>
#include <DAQ/DEventTag.h>
#include <DAQ/Df250BORConfig.h>
#include <DAQ/Df125BORConfig.h>
#include <DAQ/DF1TDCBORConfig.h>
#include <DAQ/DCAEN1290TDCBORConfig.h>
#include <DAQ/DBORptrs.h>

// Here is some C++ macro script-fu. For each type of class the DParsedEvent
// can hold, we want to have a vector of pointers to that type of object. 
// There's also a number of other things we need to do for each of these types
// but don't want to have to write the entire list out multiple times. The
// "MyTypes" macro below defines all of the types and then is used multiple
// times in the DParsedEvent class. It would be nice if we could also generate
// the above #includes using this trick but alas, the C++ language
// prohibits using #includes in macros so it's not possible.
#define MyTypes(X) \
		X(Df250Config) \
		X(Df250PulseIntegral) \
		X(Df250StreamingRawData) \
		X(Df250WindowSum) \
		X(Df250PulseRawData) \
		X(Df250TriggerTime) \
		X(Df250PulseTime) \
		X(Df250PulsePedestal) \
		X(Df250WindowRawData) \
		X(Df125Config) \
		X(Df125TriggerTime) \
		X(Df125PulseIntegral) \
		X(Df125PulseTime) \
		X(Df125PulsePedestal) \
		X(Df125PulseRawData) \
		X(Df125WindowRawData) \
		X(Df125CDCPulse) \
		X(Df125FDCPulse) \
		X(DF1TDCConfig) \
		X(DF1TDCHit) \
		X(DF1TDCTriggerTime) \
		X(DCAEN1290TDCConfig) \
		X(DCAEN1290TDCHit) \
		X(DCODAEventInfo) \
		X(DCODAROCInfo) \
		X(DTSscalers) \
		X(DEPICSvalue) \
		X(DEventTag)


class DParsedEvent{
	public:		
		
		atomic<bool> in_use;
		uint64_t Nrecycled;     // Incremented in DEVIOWorkerThread::MakeEvents()
		uint64_t MAX_RECYCLED;
		bool copied_to_factories;
		
		uint64_t istreamorder;
		uint64_t run_number;
		uint64_t event_number;
		uint64_t event_status_bits;
		bool     sync_flag;
		
		DBORptrs *borptrs;
		
		// For each type defined in "MyTypes" above, define a vector of
		// pointers to it with a name made by prepending a "v" to the classname
		// The following expands to things like e.g.
		//
		//       vector<Df250Config*> vDf250Config;
		//
		#define makevector(A) vector<A*>  v##A;
		MyTypes(makevector)
		MyBORTypes(makevector)
		
		// DParsedEvent objects are recycled to save malloc/delete cycles. Do the
		// same for the objects they provide by creating a pool vector for each
		// object type. No need for locks here since this will only ever be accessed
		// by the same worker thread.
		#define makepoolvector(A) vector<A*>  v##A##_pool;
		MyTypes(makepoolvector)

		// Method to return all objects in vectors to their respective pools and 
		// clear the vectors to set up for processing the next event. Vectors
		// with BOR types are just cleared.
		// This is called from DEVIOWorkerThread::MakeEvents
		#define returntopool(A) if(!v##A.empty()){ v##A##_pool.insert(v##A##_pool.end(), v##A.begin(), v##A.end()); }
		#define clearvectors(A)     v##A.clear();
		void Clear(void){ 
			MyTypes(returntopool)
			MyTypes(clearvectors)
			MyBORTypes(clearvectors)
		}

		// Method to delete all objects in all vectors and all pools. This should
		// usually only be called from the DParsedEvent destructor
		#define deletevector(A)     for(auto p : v##A       ) delete p;
		#define deletepool(A)       for(auto p : v##A##_pool) delete p;
		#define clearpoolvectors(A) v##A##_pool.clear();
		void Delete(void){
			MyTypes(deletevector)
			MyTypes(deletepool)
			MyTypes(clearvectors)
			MyTypes(clearpoolvectors)
			MyBORTypes(clearvectors)
		}
		
		// This is used to occasionally delete extra pool objects to reduce the
		// average memory use. It is called from DEVIOWorkerThread::MakeEvents
		// every MAX_RECYCLED events processed by this DParsedEvent object.
		void Prune(void){
			MyTypes(deletepool)
			MyTypes(clearpoolvectors)
		}
		
		// Define a class that has pointers to factories for each data type.
		// One of these is instantiated for each JEventLoop encountered.
		// See comments below for CopyToFactories for details.
		#define makefactoryptr(A) JFactory<A> *fac_##A;
		#define copyfactoryptr(A) fac_##A = (JFactory<A>*)loop->GetFactory(#A);
		class DFactoryPointers{
			public:
				JEventLoop *loop;
				MyTypes(makefactoryptr)
				MyBORTypes(makefactoryptr)

				DFactoryPointers():loop(NULL){}
				~DFactoryPointers(){}

				void Init(JEventLoop *loop){
					this->loop = loop;
					MyTypes(copyfactoryptr)
					MyBORTypes(copyfactoryptr)
				}
		};
		
		// Copy objects to factories. For efficiency, we keep an object that
		// holds the relevant factory pointers for each JEventLoop we encounter.
		// This avoids having to look up the factory pointer for each data type
		// for every event. Note that only one processing thread at a time will
		// ever call this method for this DParsedEvent object so we don't need
		// to lock access to the factory_pointers map.
		#define copytofactory(A)    facptrs.fac_##A->CopyTo(v##A);
		#define copybortofactory(A) facptrs.fac_##A->CopyTo(borptrs->v##A);
		#define setevntcalled(A)    facptrs.fac_##A->Set_evnt_called();
		#define keepownership(A)    facptrs.fac_##A->SetFactoryFlag(JFactory_base::NOT_OBJECT_OWNER);
		void CopyToFactories(JEventLoop *loop){
			// Get DFactoryPointers for this JEventLoop, creating new one if necessary
			DFactoryPointers &facptrs = factory_pointers[loop];
			if(facptrs.loop == NULL) facptrs.Init(loop);

			// Copy all data vectors to appropriate factories
			MyTypes(copytofactory)
			MyTypes(setevntcalled)
			MyTypes(keepownership)
			if(borptrs){
				MyBORTypes(copybortofactory)
				MyBORTypes(setevntcalled)
				MyBORTypes(keepownership)
			}
			copied_to_factories=true;
		}
		
		// Method to check class name against each classname in MyTypes returning
		// true if found and false if not.
		#define checkclassname(A) if(classname==#A) return true;
		bool IsParsedDataType(string &classname){
			MyTypes(checkclassname)
			MyBORTypes(checkclassname)
			return false;
		}
		
		// The following is pretty complicated to understand. What it does is
		// create a templated method for each data type that is used to
		// replace the use of "new" to allocate an object of that type.
		// What makes it complicated is the use of C++11 variadic functions
		// to allow passing in (and through) variable length argument lists.
		// This is needed since each type's constructor takes a different
		// set of arguments and we don't want to have to encode all of that
		// here.
		//
		// For each data type, a method called NEW_XXX is defined that looks
		// to see if an object already exists in the corresponding pool vector.
		// If so, it returns a pointer to it after doing an in-place constructor
		// call with the given arguments. If no objects are available in the
		// pool, then a new one is created in the normal way with the given
		// arguments.
		//
		// This will also automatically add the created/recycled object to
		// the appropriate vXXX vector as part of the current event. It
		// returns of pointer to the object so that it can be accessed by the
		// caller if needed.
		//
		// This will provide a method that looks something like this:
		//
		//  Df250TriggerTime* NEW_Df250TriggerTime(...);
		//
		// where the "..." represents whatever arguments are passed into the
		// Df250TriggerTime constructor.
		//
		#define makeallocator(A) template<typename... Args> \
		A* NEW_##A(Args&&... args){ \
			A* t = NULL; \
			if(v##A##_pool.empty()){ \
				t = new A(std::forward<Args>(args)...); \
			}else{ \
				t = v##A##_pool.back(); \
				v##A##_pool.pop_back(); \
				t->ClearAssociatedObjects(); \
				new(t) A(std::forward<Args>(args)...); \
			} \
			v##A.push_back(t); \
			return t; \
		}
		MyTypes(makeallocator);

		// Constructor and destructor
		DParsedEvent(void):in_use(false),Nrecycled(0),MAX_RECYCLED(1000),borptrs(NULL){}
		#define printcounts(A) if(!v##A.empty()) cout << v##A.size() << " : " << #A << endl;
		#define printpoolcounts(A) if(!v##A##_pool.empty()) cout << v##A##_pool.size() << " : " << #A << "_pool" << endl;
		virtual ~DParsedEvent(){
//			cout << "----- DParsedEvent (" << this << ") -------" << endl;
//			MyTypes(printcounts);
//			MyBORTypes(printcounts);
//			MyTypes(printpoolcounts);
			Delete();
		}

	protected:
		map<JEventLoop*, DFactoryPointers> factory_pointers;


};

#undef MyTypes

#endif // _DParsedEvent_

