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
using namespace jana;

#include <DAQ/DModuleType.h>
#include <DAQ/Df250PulseIntegral.h>
#include <DAQ/Df250StreamingRawData.h>
#include <DAQ/Df250WindowSum.h>
#include <DAQ/Df250PulseRawData.h>
#include <DAQ/Df250TriggerTime.h>
#include <DAQ/Df250PulseTime.h>
#include <DAQ/Df250WindowRawData.h>
#include <DAQ/DF1TDCHit.h>
#include <DAQ/DF1TDCTriggerTime.h>

#include <BCAL/DBCALHit.h>
#include <CDC/DCDCHit.h>
#include <FCAL/DFCALHit.h>
#include <FDC/DFDCHit.h>
#include <START_COUNTER/DSCHit.h>
#include <TOF/DTOFRawHit.h>


class DTranslationTable:public jana::JObject{
	public:
		JOBJECT_PUBLIC(DTranslationTable);
		
		DTranslationTable(JEventLoop *loop);
		~DTranslationTable();
		
		// Each detector system has its own native indexing scheme.
		// Here, we define a class for each of them that has those
		// indexes. These are then used below in the DChannelInfo
		// class to related them to the DAQ indexing scheme of
		// crate, slot, channel.
		typedef struct{
			uint32_t rocid;
			uint32_t slot;
			uint32_t channel;
		}csc_t;
		
		enum Detector_t{
			BCAL,
			CDC,
			FCAL,
			FDC_CATHODES,
			FDC_WIRES,
			PS,
			PSC,
			ST,
			TAGH,
			TAGM,
			TOF
		};
		
		class BCALIndex_t{
			public:
			uint32_t module;
			uint32_t layer;
			uint32_t sector;
			uint32_t end;
		};
		
		class CDCIndex_t{
			public:
			uint32_t ring;
			uint32_t straw;
		};
		
		class FCALIndex_t{
			public:
			uint32_t row;
			uint32_t col;
		};

		class FDC_CathodesIndex_t{
			public:
			uint32_t package;
			uint32_t chamber;
			uint32_t view;
			uint32_t strip;
			uint32_t strip_type;
		};

		class FDC_WiresIndex_t{
			public:
			uint32_t package;
			uint32_t chamber;
			uint32_t wire;
		};

		class PSIndex_t{
			public:
			uint32_t side;
			uint32_t id;
		};

		class PSCIndex_t{
			public:
			uint32_t id;
		};

		class STIndex_t{
			public:
			uint32_t sector;
		};

		class TAGHIndex_t{
			public:
			uint32_t id;
		};

		class TAGMIndex_t{
			public:
			uint32_t col;
			uint32_t row;
		};

		class TOFIndex_t{
			public:
			uint32_t plane;
			uint32_t bar;
			uint32_t end;
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
					STIndex_t st;
					TAGHIndex_t tagh;
					TAGMIndex_t tagm;
					TOFIndex_t tof;
				};
		};
		
		// Full translation table is collection of DChannelInfo objects
		map<csc_t, DChannelInfo> TT;
		
		// Methods
		bool IsSuppliedType(string dataClassName) const;
		void ApplyTranslationTable(jana::JEventLoop *loop) const;
		
		DBCALHit* MakeBCALHit(const Df250PulseIntegral *hit, const BCALIndex_t &idx) const;
		DCDCHit* MakeCDCHit(const Df250PulseIntegral *hit, const CDCIndex_t &idx) const;
		DBCALHit* MakeBCALHit(const Df250PulseIntegral *hit, const BCALIndex_t &idx) const;
		DBCALHit* MakeBCALHit(const Df250PulseIntegral *hit, const BCALIndex_t &idx) const;
		DBCALHit* MakeBCALHit(const Df250PulseIntegral *hit, const BCALIndex_t &idx) const;
		DBCALHit* MakeBCALHit(const Df250PulseIntegral *hit, const BCALIndex_t &idx) const;

	
	protected:
		set<string> supplied_data_types;
};

#endif // _DTranslationTable_

