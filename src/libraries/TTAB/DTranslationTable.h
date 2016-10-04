// $Id$
//
//    File: DTranslationTable.h
// Created: Thu Jun 27 15:33:38 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _DTranslationTable_
#define _DTranslationTable_

#include <set>
#include <string>
using namespace std;


#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>
#include <JANA/JCalibration.h>
#include <JANA/JException.h>
#include <JANA/JEventSource.h>
using namespace jana;

#include <DAQ/DModuleType.h>
#include <DAQ/Df250PulseData.h>
#include <DAQ/Df250PulseIntegral.h>
#include <DAQ/Df250StreamingRawData.h>
#include <DAQ/Df250WindowSum.h>
#include <DAQ/Df250PulseRawData.h>
#include <DAQ/Df250TriggerTime.h>
#include <DAQ/Df250PulseTime.h>
#include <DAQ/Df250PulsePedestal.h>
#include <DAQ/Df250WindowRawData.h>
#include <DAQ/Df125PulseIntegral.h>
#include <DAQ/Df125PulseTime.h>
#include <DAQ/Df125PulsePedestal.h>
#include <DAQ/Df125TriggerTime.h>
#include <DAQ/Df125CDCPulse.h>
#include <DAQ/Df125FDCPulse.h>
#include <DAQ/DF1TDCHit.h>
#include <DAQ/DF1TDCTriggerTime.h>
#include <DAQ/DCAEN1290TDCHit.h>

#include <BCAL/DBCALDigiHit.h>
#include <BCAL/DBCALTDCDigiHit.h>
#include <CDC/DCDCDigiHit.h>
#include <FCAL/DFCALDigiHit.h>
#include <FDC/DFDCCathodeDigiHit.h>
#include <FDC/DFDCWireDigiHit.h>
#include <RF/DRFDigiTime.h>
#include <RF/DRFTDCDigiTime.h>
#include <START_COUNTER/DSCDigiHit.h>
#include <START_COUNTER/DSCTDCDigiHit.h>
#include <TOF/DTOFDigiHit.h>
#include <TOF/DTOFTDCDigiHit.h>
#include <TAGGER/DTAGMDigiHit.h>
#include <TAGGER/DTAGMTDCDigiHit.h>
#include <TAGGER/DTAGHDigiHit.h>
#include <TAGGER/DTAGHTDCDigiHit.h>
#include <PAIR_SPECTROMETER/DPSDigiHit.h>
#include <PAIR_SPECTROMETER/DPSCDigiHit.h>
#include <PAIR_SPECTROMETER/DPSCTDCDigiHit.h>
#include <TPOL/DTPOLSectorDigiHit.h>

// (See comments in DParsedEvent.h for enlightenment)
#define MyTypes(X) \
		X(DBCALDigiHit) \
		X(DBCALTDCDigiHit) \
		X(DCDCDigiHit) \
		X(DFCALDigiHit) \
		X(DFDCCathodeDigiHit) \
		X(DFDCWireDigiHit) \
		X(DRFDigiTime) \
		X(DRFTDCDigiTime) \
		X(DSCDigiHit) \
		X(DSCTDCDigiHit) \
		X(DTOFDigiHit) \
		X(DTOFTDCDigiHit) \
		X(DTAGMDigiHit) \
		X(DTAGMTDCDigiHit) \
		X(DTAGHDigiHit) \
		X(DTAGHTDCDigiHit) \
		X(DPSDigiHit) \
		X(DPSCDigiHit) \
		X(DPSCTDCDigiHit) \
		X(DTPOLSectorDigiHit)


#include "GlueX.h"

class DTranslationTable:public jana::JObject{
	public:
		JOBJECT_PUBLIC(DTranslationTable);
		
		DTranslationTable(JEventLoop *loop);
		~DTranslationTable();
		
		// Each detector system has its own native indexing scheme.
		// Here, we define a class for each of them that has those
		// indexes. These are then used below in the DChannelInfo
		// class to relate them to the DAQ indexing scheme of
		// crate, slot, channel.
		struct csc_t{
			uint32_t rocid;
			uint32_t slot;
			uint32_t channel;

			inline bool operator==(const struct csc_t &rhs) const {
			    return (rocid==rhs.rocid) && (slot==rhs.slot) && (channel==rhs.channel);
			}
		};
		
		enum Detector_t{
			UNKNOWN_DETECTOR,
			BCAL,
			CDC,
			FCAL,
			FDC_CATHODES,
			FDC_WIRES,
			PS,
			PSC,
			RF,
			SC,
			TAGH,
			TAGM,
			TOF,
			TPOLSECTOR,
			NUM_DETECTOR_TYPES
		};

		string DetectorName(Detector_t type) const {
			switch(type){
				case BCAL: return "BCAL";
				case CDC: return "CDC";
				case FCAL: return "FCAL";
				case FDC_CATHODES: return "FDC_CATHODES";
				case FDC_WIRES: return "FDC_WIRES";
				case PS: return "PS";
				case PSC: return "PSC";
				case RF: return "RF";
				case SC: return "SC";
				case TAGH: return "TAGH";
				case TAGM: return "TAGM";
				case TOF: return "TOF";
			        case TPOLSECTOR: return "TPOL"; // is set to TPOL to match what is in CCDB, fix later
				case UNKNOWN_DETECTOR:
				default:
					return "UNKNOWN";
			}
		}

		class BCALIndex_t{
			public:
			uint32_t module;
			uint32_t layer;
			uint32_t sector;
			uint32_t end;

			inline bool operator==(const BCALIndex_t &rhs) const {
			    return (module==rhs.module) && (layer==rhs.layer) 
			      && (sector==rhs.sector) && (end==rhs.end);
			}
		};
		
		class CDCIndex_t{
			public:
			uint32_t ring;
			uint32_t straw;

			inline bool operator==(const CDCIndex_t &rhs) const {
			    return (ring==rhs.ring) && (straw==rhs.straw);
			}
		};
		
		class FCALIndex_t{
			public:
			uint32_t row;
			uint32_t col;

			inline bool operator==(const FCALIndex_t &rhs) const {
			    return (row==rhs.row) && (col==rhs.col);
			}
		};

		class FDC_CathodesIndex_t{
			public:
			uint32_t package;
			uint32_t chamber;
			uint32_t view;
			uint32_t strip;
			uint32_t strip_type;

			inline bool operator==(const FDC_CathodesIndex_t &rhs) const {
			    return (package==rhs.package) && (chamber==rhs.chamber) && (view==rhs.view)
			      && (strip==rhs.strip) && (strip_type==rhs.strip_type);
			}
		};

		class FDC_WiresIndex_t{
			public:
			uint32_t package;
			uint32_t chamber;
			uint32_t wire;

			inline bool operator==(const FDC_WiresIndex_t &rhs) const {
			    return (package==rhs.package) && (chamber==rhs.chamber) && (wire==rhs.wire);
			}
		};

		class PSIndex_t{
			public:
			uint32_t side;
			uint32_t id;

			inline bool operator==(const PSIndex_t &rhs) const {
			    return (side==rhs.side) && (id==rhs.id);
			}
		};

		class PSCIndex_t{
			public:
			uint32_t id;

			inline bool operator==(const PSCIndex_t &rhs) const {
			    return (id==rhs.id);
			}
		};

		class RFIndex_t{
			public:
			DetectorSystem_t dSystem;

			inline bool operator==(const RFIndex_t &rhs) const {
				return (dSystem == rhs.dSystem);
			}
		};

		class SCIndex_t{
			public:
			uint32_t sector;

			inline bool operator==(const SCIndex_t &rhs) const {
			    return (sector==rhs.sector);
			}
		};

		class TAGHIndex_t{
			public:
			uint32_t id;

			inline bool operator==(const TAGHIndex_t &rhs) const {
			    return (id==rhs.id);
			}
		};

		class TAGMIndex_t{
			public:
			uint32_t col;
			uint32_t row;

			inline bool operator==(const TAGMIndex_t &rhs) const {
			    return (col==rhs.col) && (row==rhs.row);
			}
		};

		class TOFIndex_t{
			public:
			uint32_t plane;
			uint32_t bar;
			uint32_t end;

			inline bool operator==(const TOFIndex_t &rhs) const {
			    return (plane==rhs.plane) && (bar==rhs.bar) && (end==rhs.end);
			}
		};
		
		class TPOLSECTORIndex_t{
			public:
			uint32_t sector;

			inline bool operator==(const TPOLSECTORIndex_t &rhs) const {
			    return (sector==rhs.sector);
			}
		};

		// DChannelInfo holds translation between indexing schemes
		// for one channel.
		class DChannelInfo{
			public:
				csc_t CSC;
				DModuleType::type_id_t module_type;
				Detector_t det_sys;
				union{
					BCALIndex_t bcal;
					CDCIndex_t cdc;
					FCALIndex_t fcal;
					FDC_CathodesIndex_t fdc_cathodes;
					FDC_WiresIndex_t fdc_wires;
					PSIndex_t ps;
					PSCIndex_t psc;
					RFIndex_t rf;
					SCIndex_t sc;
					TAGHIndex_t tagh;
					TAGMIndex_t tagm;
					TOFIndex_t tof;
					TPOLSECTORIndex_t tpolsector;
				};
		};

		//-----------------------------------------------------------------------
		//
		// Pre-processor macro monkey shines using the MyTypes define above to repeat
		// lots of code for all types in a compact and robust way.
		//
		// For each type defined in "MyTypes" above, define a vector of
		// pointers to it with a name made by prepending a "v" to the classname
		// The following expands to things like e.g.
		//
		//       vector<DBCALDigiHit*> vDBCALDigiHit;
		//
		#define makevector(A) mutable vector<A*>  v##A;
		MyTypes(makevector)
		
		// Similarly, define a pointer to the factory for each type.
		#define makefactoryptr(A) JFactory<A> *fac_##A;
		MyTypes(makefactoryptr)
		
		// Method to initialize factory pointers
		#define copyfactoryptr(A) fac_##A = (JFactory<A>*)loop->GetFactory(#A);
		void InitFactoryPointers(JEventLoop *loop){ MyTypes(copyfactoryptr) }

		// Method to clear each of the vectors at beginning of event
		#define clearvector(A) v##A.clear();
		void ClearVectors(void) const { MyTypes(clearvector) }

		// Method to copy all produced objects to respective factories
		#define copytofactory(A) fac_##A->CopyTo(v##A);
		void CopyToFactories(void) const { MyTypes(copytofactory) }
		
		// Method to check class name against each classname in MyTypes returning
		// true if found and false if not.
		#define checkclassname(A) if(classname==#A) return true;
		bool IsSuppliedType(string &classname) const {
			MyTypes(checkclassname)
			return false;
		}
		
		// Method to print sizes of all vectors (for debugging)
		#define printvectorsize(A) ttout << "     v" #A ".size() = " << v##A.size() << std::endl;
		void PrintVectorSizes(void) const { MyTypes(printvectorsize) }

		//-----------------------------------------------------------------------

		// Methods
		void ApplyTranslationTable(jana::JEventLoop *loop) const;
		
		// fADC250 -- Fall 2016 -> ?
		DBCALDigiHit*       MakeBCALDigiHit(       const BCALIndex_t &idx,       const Df250PulseData *pd) const;
		DFCALDigiHit*       MakeFCALDigiHit(       const FCALIndex_t &idx,       const Df250PulseData *pd) const;
		DSCDigiHit*         MakeSCDigiHit(         const SCIndex_t &idx,         const Df250PulseData *pd) const;
		DTOFDigiHit*        MakeTOFDigiHit(        const TOFIndex_t &idx,        const Df250PulseData *pd) const;
		DTAGMDigiHit*       MakeTAGMDigiHit(       const TAGMIndex_t &idx,       const Df250PulseData *pd) const;
		DTAGHDigiHit*       MakeTAGHDigiHit(       const TAGHIndex_t &idx,       const Df250PulseData *pd) const;
		DPSDigiHit*         MakePSDigiHit(         const PSIndex_t &idx,         const Df250PulseData *pd) const;
		DPSCDigiHit*        MakePSCDigiHit(        const PSCIndex_t &idx,        const Df250PulseData *pd) const;
		DRFDigiTime*        MakeRFDigiTime(        const RFIndex_t &idx,         const Df250PulseData *pd) const;
		DTPOLSectorDigiHit* MakeTPOLSectorDigiHit( const TPOLSECTORIndex_t &idx, const Df250PulseData *pd) const;

		// fADC250 -- commissioning -> Fall 2016
		DBCALDigiHit*       MakeBCALDigiHit(const BCALIndex_t &idx, const Df250PulseIntegral *pi, const Df250PulseTime *pt, const Df250PulsePedestal *pp) const;
		DFCALDigiHit*       MakeFCALDigiHit(const FCALIndex_t &idx, const Df250PulseIntegral *pi, const Df250PulseTime *pt, const Df250PulsePedestal *pp) const;
		DSCDigiHit*         MakeSCDigiHit(  const SCIndex_t &idx,   const Df250PulseIntegral *pi, const Df250PulseTime *pt, const Df250PulsePedestal *pp) const;
		DTOFDigiHit*        MakeTOFDigiHit( const TOFIndex_t &idx,  const Df250PulseIntegral *pi, const Df250PulseTime *pt, const Df250PulsePedestal *pp) const;
		DTAGMDigiHit*       MakeTAGMDigiHit(const TAGMIndex_t &idx, const Df250PulseIntegral *pi, const Df250PulseTime *pt, const Df250PulsePedestal *pp) const;
		DTAGHDigiHit*       MakeTAGHDigiHit(const TAGHIndex_t &idx, const Df250PulseIntegral *pi, const Df250PulseTime *pt, const Df250PulsePedestal *pp) const;
		DPSDigiHit*         MakePSDigiHit(  const PSIndex_t &idx,   const Df250PulseIntegral *pi, const Df250PulseTime *pt, const Df250PulsePedestal *pp) const;
		DPSCDigiHit*        MakePSCDigiHit( const PSCIndex_t &idx,  const Df250PulseIntegral *pi, const Df250PulseTime *pt, const Df250PulsePedestal *pp) const;
		DRFDigiTime*        MakeRFDigiTime( const RFIndex_t &idx,   const Df250PulseTime *hit) const;
		DTPOLSectorDigiHit* MakeTPOLSectorDigiHit(const TPOLSECTORIndex_t &idx, const Df250PulseIntegral *pi, const Df250PulseTime *pt, const Df250PulsePedestal *pp) const;

		// fADC125
		DCDCDigiHit* MakeCDCDigiHit(const CDCIndex_t &idx, const Df125PulseIntegral *pi, const Df125PulseTime *pt, const Df125PulsePedestal *pp) const;
		DCDCDigiHit* MakeCDCDigiHit(const CDCIndex_t &idx, const Df125CDCPulse *p) const;
		DFDCCathodeDigiHit* MakeFDCCathodeDigiHit(const FDC_CathodesIndex_t &idx, const Df125PulseIntegral *pi, const Df125PulseTime *pt, const Df125PulsePedestal *pp) const;
		DFDCCathodeDigiHit* MakeFDCCathodeDigiHit(const FDC_CathodesIndex_t &idx, const Df125FDCPulse *p) const;

		// F1TDC
		DBCALTDCDigiHit* MakeBCALTDCDigiHit(const BCALIndex_t &idx,      const DF1TDCHit *hit) const;
		DFDCWireDigiHit* MakeFDCWireDigiHit(const FDC_WiresIndex_t &idx, const DF1TDCHit *hit) const;
		DRFTDCDigiTime*  MakeRFTDCDigiTime( const RFIndex_t &idx,        const DF1TDCHit *hit) const;
		DSCTDCDigiHit*   MakeSCTDCDigiHit(  const SCIndex_t &idx,        const DF1TDCHit *hit) const;
		DTAGMTDCDigiHit* MakeTAGMTDCDigiHit(const TAGMIndex_t &idx,      const DF1TDCHit *hit) const;
		DTAGHTDCDigiHit* MakeTAGHTDCDigiHit(const TAGHIndex_t &idx,      const DF1TDCHit *hit) const;
		DPSCTDCDigiHit*  MakePSCTDCDigiHit( const PSCIndex_t &idx,       const DF1TDCHit *hit) const;
		
		// CAEN1290TDC
		DTOFTDCDigiHit*  MakeTOFTDCDigiHit(const TOFIndex_t &idx,        const DCAEN1290TDCHit *hit) const;
		DRFTDCDigiTime*  MakeRFTDCDigiTime(const RFIndex_t &idx,         const DCAEN1290TDCHit *hit) const;

		void Addf250ObjectsToCallStack(JEventLoop *loop, string caller) const;
		void Addf125CDCObjectsToCallStack(JEventLoop *loop, string caller, bool addpulseobjs) const;
		void Addf125FDCObjectsToCallStack(JEventLoop *loop, string caller, bool addpulseobjs) const;
		void AddF1TDCObjectsToCallStack(JEventLoop *loop, string caller) const;
		void AddCAEN1290TDCObjectsToCallStack(JEventLoop *loop, string caller) const;
		void AddToCallStack(JEventLoop *loop, string caller, string callee) const;

		void ReadOptionalROCidTranslation(void);
		void SetSystemsToParse(string systems, JEventSource *eventsource);
		void SetSystemsToParse(JEventSource *eventsource){SetSystemsToParse(SYSTEMS_TO_PARSE, eventsource);}
		void ReadTranslationTable(JCalibration *jcalib=NULL);
		
		template<class T> void CopyDf250Info(T *h, const Df250PulseIntegral *pi, const Df250PulseTime *pt, const Df250PulsePedestal *pp) const;
		template<class T> void CopyDf250Info(T *h, const Df250PulseData *pd) const;
		template<class T> void CopyDf125Info(T *h, const Df125PulseIntegral *pi, const Df125PulseTime *pt, const Df125PulsePedestal *pp) const;
		template<class T> void CopyDF1TDCInfo(T *h, const DF1TDCHit *hit) const;
		template<class T> void CopyDCAEN1290TDCInfo(T *h, const DCAEN1290TDCHit *hit) const;

		
		// methods for others to search the Translation Table
		const DChannelInfo &GetDetectorIndex(const csc_t &in_daq_index) const;
		const csc_t &GetDAQIndex(const DChannelInfo &in_channel) const;

		//public so that StartElement can access it
		static map<DTranslationTable::Detector_t, set<uint32_t> >& Get_ROCID_By_System(void); //this is static so that StartElement can access it

		// This was left over from long ago and is not currently used. It seems
		// potentially useful though in the future so I don't want to get rid
		// of it. However, it was causing "defined but not used" warnings when
		// defined as a simple C function. Moved it here just to prevent those
		// warnings.
		static int ModuleStr2ModID(string &type)
		{
   		if (type == "vmecpu") {
      		return(DModuleType::VMECPU);
   		} else if (type == "tid") {
      		return(DModuleType::TID);
   		} else if (type == "fadc250") {
      		return(DModuleType::FADC250);
   		} else if (type == "fadc125") {
      		return(DModuleType::FADC125);
   		} else if (type == "f1tdcv2") {
      		return(DModuleType::F1TDC32);
   		} else if (type == "f1tdcv3") {
      		return(DModuleType::F1TDC48);
   		} else if (type == "jldisc") {
      		return(DModuleType::JLAB_DISC);
   		} else if (type == "vx1290a") {
      		return(DModuleType::CAEN1290);
   		} else {
      		return(DModuleType::UNKNOWN);
   		}
		}

	protected:
		string XML_FILENAME;
		bool NO_CCDB;
		set<string> supplied_data_types;
		int VERBOSE;
		string SYSTEMS_TO_PARSE;
		string ROCID_MAP_FILENAME;
		bool CALL_STACK;
		
		mutable JStreamLog ttout;

		string Channel2Str(const DChannelInfo &in_channel) const;

	private:

		/****************************************** STATIC-VARIABLE-ACCESSING PRIVATE MEMBER FUNCTIONS ******************************************/

		//Some variables needs to be shared amongst threads (e.g. the memory used for the branch variables)
		//However, you cannot make them global/extern/static/static-member variables in the header file:
			//The header file is included in the TTAB library AND in each plugin that uses it
				//When a header file is included in a src file, it's contents are essentially copied directly into it
			//Thus there are two instances of each static variable: one in each translation unit (library)
			//Supposedly(?) they are linked together during runtime when loading, so there is (supposedly) no undefined behavior.
			//However, this causes a double free (double-deletion) when these libraries are closed at the end of the program, crashing it.
		//Thus the variables must be in a single source file that is compiled into a single library
		//However, you (somehow?) cannot make them global/extern variables in the cpp function
			//This also (somehow?) causes the double-free problem above for (at least) stl containers
			//It works for pointers-to-stl-containers and fundamental types, but I dunno why.
			//It's not good encapsulation anyway though.
		//THE SOLUTION:
			//Define the variables as static, in the source file, WITHIN A PRIVATE MEMBER FUNCTION.
			//Thus the static variables themselves only have function scope.
			//Access is only available via the private member functions, thus access is fully controlled.
			//They are shared amongst threads, so locks are necessary, but since they are private this class can handle it internally

		pthread_mutex_t& Get_TT_Mutex(void) const;
		bool& Get_TT_Initialized(void) const;
		map<DTranslationTable::csc_t, DTranslationTable::DChannelInfo>& Get_TT(void) const;
		map<uint32_t, uint32_t>& Get_ROCID_Map(void) const;
		map<uint32_t, uint32_t>& Get_ROCID_Inv_Map(void) const;
};

//---------------------------------
// CopyDf250Info
//---------------------------------
template<class T>
void DTranslationTable::CopyDf250Info(T *h, const Df250PulseIntegral *pi, const Df250PulseTime *pt, const Df250PulsePedestal *pp) const
{
	/// Copy info from the fADC250 into a hit object.
	h->pulse_integral    = pi->integral;
	h->pulse_time        = pt==NULL ? 0:pt->time;
	h->pulse_peak        = pp==NULL ? 0:pp->pulse_peak;
	h->pedestal          = pi->pedestal;
	h->QF                = pi->quality_factor;
	h->nsamples_integral = pi->nsamples_integral;
	h->nsamples_pedestal = pi->nsamples_pedestal;
	h->datasource        = pi->emulated ? 0:1;
	
	h->AddAssociatedObject(pi);
	if(pt) h->AddAssociatedObject(pt);
	if(pp) h->AddAssociatedObject(pp);
}

//---------------------------------
// CopyDf250Info
//---------------------------------
template<class T>
void DTranslationTable::CopyDf250Info(T *h, const Df250PulseData *pd) const
{
	/// Copy info from the fADC250 into a hit object.
	h->pulse_integral          = pd->integral;
	h->pulse_time              = (pd->course_time<<6) + pd->fine_time;
	h->pulse_peak              = pd->pulse_peak;
	h->pedestal                = pd->pedestal;
	h->nsamples_integral       = pd->nsamples_integral;
	h->nsamples_pedestal       = pd->nsamples_pedestal;
	h->datasource              = pd->emulated ? 0:2;

	uint32_t QF = 0; // make single quality factor number for compactness
	if( pd->QF_pedestal         ) QF |= (1<<0);
	if( pd->QF_NSA_beyond_PTW   ) QF |= (1<<1);
	if( pd->QF_overflow         ) QF |= (1<<2);
	if( pd->QF_underflow        ) QF |= (1<<3);
	if( pd->QF_vpeak_beyond_NSA ) QF |= (1<<4);
	if( pd->QF_vpeak_not_found  ) QF |= (1<<5);
	if( pd->QF_bad_pedestal     ) QF |= (1<<6);
	h->QF = QF;
	
	h->AddAssociatedObject(pd);
}

//---------------------------------
// CopyDf125Info
//---------------------------------
template<class T>
void DTranslationTable::CopyDf125Info(T *h, const Df125PulseIntegral *pi, const Df125PulseTime *pt, const Df125PulsePedestal *pp) const
{
	/// Copy info from the fADC125 into a hit object.
	h->pulse_integral    = pi->integral;
	h->pulse_time        = pt==NULL ? 0:pt->time;
	h->pedestal          = pi->pedestal;
	h->QF                = pi->quality_factor;
	h->nsamples_integral = pi->nsamples_integral;
	h->nsamples_pedestal = pi->nsamples_pedestal;
	
	h->AddAssociatedObject(pi);
	if(pt) h->AddAssociatedObject(pt);
	if(pp) h->AddAssociatedObject(pp);
}

//---------------------------------
// CopyDF1TDCInfo
//---------------------------------
template<class T>
void DTranslationTable::CopyDF1TDCInfo(T *h, const DF1TDCHit *hit) const
{
	/// Copy info from the f1tdc into a hit object.

	h->time = hit->time;
	h->AddAssociatedObject(hit);
}

//---------------------------------
// CopyDCAEN1290TDCInfo
//---------------------------------
template<class T>
void DTranslationTable::CopyDCAEN1290TDCInfo(T *h, const DCAEN1290TDCHit *hit) const
{
	/// Copy info from the CAEN1290 into a hit object.
	h->time = hit->time;
	
	h->AddAssociatedObject(hit);
}

#undef MyTypes
#undef makevector
#undef makefactoryptr
#undef copyfactoryptr
#undef clearvector
#undef copytofactory
#undef checkclassname
#undef printvectorsize

#endif // _DTranslationTable_

