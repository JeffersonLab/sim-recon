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
using namespace jana;

#include <DAQ/DModuleType.h>
#include <DAQ/Df250PulseIntegral.h>
#include <DAQ/Df250StreamingRawData.h>
#include <DAQ/Df250WindowSum.h>
#include <DAQ/Df250PulseRawData.h>
#include <DAQ/Df250TriggerTime.h>
#include <DAQ/Df250PulseTime.h>
#include <DAQ/Df250WindowRawData.h>
#include <DAQ/Df125PulseIntegral.h>
#include <DAQ/Df125PulseTime.h>
#include <DAQ/Df125TriggerTime.h>
#include <DAQ/DF1TDCHit.h>
#include <DAQ/DF1TDCTriggerTime.h>
#include <DAQ/DCAEN1290TDCHit.h>

#include <BCAL/DBCALDigiHit.h>
#include <BCAL/DBCALTDCDigiHit.h>
#include <CDC/DCDCDigiHit.h>
#include <FCAL/DFCALDigiHit.h>
#include <FDC/DFDCCathodeDigiHit.h>
#include <FDC/DFDCWireDigiHit.h>
#include <START_COUNTER/DSCDigiHit.h>
#include <START_COUNTER/DSCTDCDigiHit.h>
#include <TOF/DTOFDigiHit.h>
#include <TOF/DTOFTDCDigiHit.h>
#include <TAGGER/DTAGMDigiHit.h>
#include <TAGGER/DTAGMTDCDigiHit.h>
#include <TAGGER/DTAGHDigiHit.h>
#include <TAGGER/DTAGHTDCDigiHit.h>


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
			SC,
			TAGH,
			TAGM,
			TOF
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
				case SC: return "SC";
				case TAGH: return "TAGH";
				case TAGM: return "TAGM";
				case TOF: return "TOF";
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
					SCIndex_t sc;
					TAGHIndex_t tagh;
					TAGMIndex_t tagm;
					TOFIndex_t tof;
				};
		};
		
		// Full translation table is collection of DChannelInfo objects
		//map<csc_t, DChannelInfo> TT;
		
		// Methods
		bool IsSuppliedType(string dataClassName) const;
		void ApplyTranslationTable(jana::JEventLoop *loop) const;
		
		// fADC250
		DBCALDigiHit* MakeBCALDigiHit(const BCALIndex_t &idx, const Df250PulseIntegral *pi, const Df250PulseTime *pt) const;
		DFCALDigiHit* MakeFCALDigiHit(const FCALIndex_t &idx, const Df250PulseIntegral *pi, const Df250PulseTime *pt) const;
		DSCDigiHit*   MakeSCDigiHit(  const SCIndex_t &idx,   const Df250PulseIntegral *pi, const Df250PulseTime *pt) const;
		DTOFDigiHit*  MakeTOFDigiHit( const TOFIndex_t &idx,  const Df250PulseIntegral *pi, const Df250PulseTime *pt) const;
		DTAGMDigiHit* MakeTAGMDigiHit(const TAGMIndex_t &idx, const Df250PulseIntegral *pi, const Df250PulseTime *pt) const;
		DTAGHDigiHit* MakeTAGHDigiHit(const TAGHIndex_t &idx, const Df250PulseIntegral *pi, const Df250PulseTime *pt) const;

		// fADC125
		DCDCDigiHit* MakeCDCDigiHit(const CDCIndex_t &idx, const Df125PulseIntegral *pi, const Df125PulseTime *pt) const;
		DFDCCathodeDigiHit* MakeFDCCathodeDigiHit(const FDC_CathodesIndex_t &idx, const Df125PulseIntegral *pi, const Df125PulseTime *pt) const;

		// F1TDC
		DBCALTDCDigiHit* MakeBCALTDCDigiHit(const BCALIndex_t &idx,      const DF1TDCHit *hit) const;
		DFDCWireDigiHit* MakeFDCWireDigiHit(const FDC_WiresIndex_t &idx, const DF1TDCHit *hit) const;
		DSCTDCDigiHit*   MakeSCTDCDigiHit(  const SCIndex_t &idx,        const DF1TDCHit *hit) const;
		DTAGMTDCDigiHit* MakeTAGMTDCDigiHit(const TAGMIndex_t &idx,      const DF1TDCHit *hit) const;
		DTAGHTDCDigiHit* MakeTAGHTDCDigiHit(const TAGHIndex_t &idx,      const DF1TDCHit *hit) const;
		
		// CAEN1290TDC
		DTOFTDCDigiHit*  MakeTOFTDCDigiHit(const TOFIndex_t &idx,        const DCAEN1290TDCHit *hit) const;

		void AddToCallStack(JEventLoop *loop, string caller, string callee) const;

		void ReadOptionalROCidTranslation(void);
		void ReadTranslationTable(JCalibration *jcalib=NULL);
		
		template<class T> void CopyToFactory(JEventLoop *loop, vector<T*> &v) const;
		template<class T> void CopyDf250Info(T *h, const Df250PulseIntegral *pi, const Df250PulseTime *pt) const;
		template<class T> void CopyDf125Info(T *h, const Df125PulseIntegral *pi, const Df125PulseTime *pt) const;
		template<class T> void CopyDF1TDCInfo(T *h, const DF1TDCHit *hit) const;
		template<class T> void CopyDCAEN1290TDCInfo(T *h, const DCAEN1290TDCHit *hit) const;

		
		// methods for others to search the Translation Table
		const DChannelInfo &GetDetectorIndex(const csc_t &in_daq_index) const;
		const csc_t &GetDAQIndex(const DChannelInfo &in_channel) const;


	protected:
		string XML_FILENAME;
		bool NO_CCDB;
		set<string> supplied_data_types;
		int VERBOSE;
		
		mutable JStreamLog ttout;

		string Channel2Str(const DChannelInfo &in_channel) const;
};

//---------------------------------
// CopyToFactory
//---------------------------------
template<class T>
void DTranslationTable::CopyToFactory(JEventLoop *loop, vector<T*> &v) const
{
	/// Template method for copying values from local containers into
	/// factories. This is done in a template since the type appears
	/// in at least 3 places below. It makes the code calling this
	/// more succinct and therefore easier to add new types.

	// It would be a little safer to use a dynamic_cast here, but
	// all documentation seems to discourage using that as it is
	// inefficient.
	JFactory<T> *fac = (JFactory<T> *)loop->GetFactory(T::static_className());
	if(VERBOSE>8) ttout << " Copying " << T::static_className() << " objects to factory: " << hex << fac << dec << endl;
	if(fac) fac->CopyTo(v);
}

//---------------------------------
// CopyDf250Info
//---------------------------------
template<class T>
void DTranslationTable::CopyDf250Info(T *h, const Df250PulseIntegral *pi, const Df250PulseTime *pt) const
{
	/// Copy info from the fADC250 into a hit object.
	h->pulse_integral    = pi->integral;
	h->pulse_time        = pt==NULL ? 0:pt->time;
	h->pedestal          = pi->pedestal;
	h->QF                = pi->quality_factor;
	h->nsamples_integral = pi->nsamples_integral;
	h->nsamples_pedestal = pi->nsamples_pedestal;
	
	h->AddAssociatedObject(pi);
	if(pt) h->AddAssociatedObject(pt);
}

//---------------------------------
// CopyDf125Info
//---------------------------------
template<class T>
void DTranslationTable::CopyDf125Info(T *h, const Df125PulseIntegral *pi, const Df125PulseTime *pt) const
{
	/// Copy info from the fADC125 into a hit object.
	h->pulse_integral = pi->integral;
	h->pulse_time     = pt==NULL ? 0:pt->time;
	h->pedestal       = 0; // for future development
	h->QF             = pi->quality_factor;
	
	h->AddAssociatedObject(pi);
	if(pt) h->AddAssociatedObject(pt);
}

//---------------------------------
// CopyDF1TDCInfo
//---------------------------------
template<class T>
void DTranslationTable::CopyDF1TDCInfo(T *h, const DF1TDCHit *hit) const
{
	/// Copy info from the fADC125 into a hit object.
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


#endif // _DTranslationTable_

