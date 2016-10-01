// $Id:$
//
//    File: DMapEVIOWords.h
// Created: Sat May 28 19:22:47 EDT 2016
// Creator: davidl (on Linux gluon104.jlab.org 2.6.32-358.23.2.el6.x86_64)
//

#ifndef _DMapEVIOWords_
#define _DMapEVIOWords_

#include <TDirectory.h>

#include <JANA/JEventProcessor.h>
#include <JANA/JEvent.h>

class DMapEVIOWords{
	public:

		enum EVIOWordType{
			kSpacerBefore,

			kUnknown,
			kEVIOHeader,
			kEVIOEventNumber,
			kEVIOTimestamp,
			
			kSpacer0,
			
			kBORData,
			
			kSpacer1,
			
			kf250BlockHeader,
			kf250BlockTrailer,
			kf250EventHeader,
			kf250TriggerTime,
			kf250WindowRawData,
			kf250WindowSum,
			kf250PulseRawData,
			kf250PulseData,
			kf250PulseIntegral,
			kf250PulseTime,
			kf250PulsePedestal,
			kf250EventTrailer,
			kf250DataNotValid,
			kf250Filler,
			kf250Unknown,
			
			kSpacer2,

			kf125BlockHeader,
			kf125BlockTrailer,
			kf125EventHeader,
			kf125TriggerTime,
			kf125WindowRawData,
			kf125CDCPulse,
			kf125FDCPulse6,
			kf125FDCPulse9,
			kf125PulseIntegral,
			kf125PulseTime,
			kf125PulsePedestal,
			kf125EventTrailer,
			kf125DataNotValid,
			kf125Filler,
			kf125Unknown,
			
			kSpacer3,
			
			kF1v2BlockHeader,
			kF1v2BLockTrailer,
			kF1v2EventHeader,
			kF1v2TriggerTime,
			kF1v2ChipHeader,
			kF1v2Data,
			kF1v2Filler,
			kF1v2BreakWord,
			kF1v2Unknown,
			
			kSpacer4,

			kF1v3BlockHeader,
			kF1v3BLockTrailer,
			kF1v3EventHeader,
			kF1v3TriggerTime,
			kF1v3ChipHeader,
			kF1v3Data,
			kF1v3Filler,
			kF1v3BreakWord,
			kF1v3Unknown,
			
			kSpacer5,

			kCAEN1190GlobalHeader,
			kCAEN1190GlobalTrailer,
			kCAEN1190GlobalTriggerTime,
			kCAEN1190TDCHeader,
			kCAEN1190TDCData,
			kCAEN1190TDCError,
			kCAEN1190TDCTrailer,
			kCAEN1190Filler,
			kCAEN1190Unknown,
			
			kSpacer6,
			
			kConfig,
			kConfigf250,
			kConfigf125,
			kConfigF1,
			kConfigCAEN1190,
			
			kSpacer7,

			kEPICSheader,
			kEPICSdata,
			
			kSpacer8,

			kF800FAFA,
			kD00DD00D,
			
			kTotWords,
			kNevents,
			
			kNEVIOWordTypes
		};

		DMapEVIOWords();
		~DMapEVIOWords();
		void AddROCIDLabels(void);
		void ParseEvent(uint32_t *buff);
		void DataWordStats(uint32_t *iptr, uint32_t *iend, uint32_t *word_stats);

		uint32_t TotWordCount(uint32_t *word_stats);

		void ParseJLabModuleData(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats);
		void Parsef250Bank(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats);
		void Parsef125Bank(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats);
		void ParseF1v2TDCBank(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats);
		void ParseF1v3TDCBank(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats);
		void ParseCAEN1190(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats);
		void ParseModuleConfiguration(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats);

		set<uint64_t> ts_history;

	private:
};

#endif // _DMapEVIOWords_

