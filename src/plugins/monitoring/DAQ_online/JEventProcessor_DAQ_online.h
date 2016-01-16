// $Id$
//
//    File: JEventProcessor_DAQ_online.h
// Created: Thu Aug  7 09:37:01 EDT 2014
// Creator: dalton (on Linux gluon05.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_DAQ_online_
#define _JEventProcessor_DAQ_online_

#include <TDirectory.h>

#include <JANA/JEventProcessor.h>
#include <JANA/JEvent.h>

class JEventProcessor_DAQ_online:public jana::JEventProcessor{
	public:

		enum EVIOWordType{
			kSpacerBefore,

			kUnknown,			
			kEVIOEventNumber,
			kEVIOTimestamp,
			
			kSpacer0,
			
			kf250BlockHeader,
			kf250BlockTrailer,
			kf250EventHeader,
			kf250TriggerTime,
			kf250WindowRawData,
			kf250WindowSum,
			kf250PulseRawData,
			kf250PulseIntegral,
			kf250PulseTime,
			kf250PulsePedestal,
			kf250EventTrailer,
			kf250DataNotValid,
			kf250Filler,
			
			kSpacer1,

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
			
			kSpacer2,
			
			kF1v2BlockHeader,
			kF1v2BLockTrailer,
			kF1v2EventHeader,
			kF1v2TriggerTime,
			kF1v2ChipHeader,
			kF1v2Data,
			kF1v2Filler,
			kF1v2BreakWord,
			
			kSpacer3,

			kF1v3BlockHeader,
			kF1v3BLockTrailer,
			kF1v3EventHeader,
			kF1v3TriggerTime,
			kF1v3ChipHeader,
			kF1v3Data,
			kF1v3Filler,
			kF1v3BreakWord,
			
			kSpacer4,

			kCAEN1190GlobalHeader,
			kCAEN1190GlobalTrailer,
			kCAEN1190GlobalTriggerTime,
			kCAEN1190TDCHeader,
			kCAEN1190TDCData,
			kCAEN1190TDCError,
			kCAEN1190TDCTrailer,
			kCAEN1190Filler,
			
			kSpacer5,
			
			kConfig,
			kConfigf250,
			kConfigf125,
			kConfigF1,
			kConfigCAEN1190,
			
			kSpacer6,

			kEPICSheader,
			kEPICSdata,
			
			kSpacer7,

			kF800FAFA,
			kD00DD00D,
			
			kTotWords,
			kNevents,
			
			kNEVIOWordTypes
		};

		JEventProcessor_DAQ_online();
		~JEventProcessor_DAQ_online();
		const char* className(void){return "JEventProcessor_DAQ_online";}

		TDirectory *maindir;
		TDirectory *daqdir;
		
		void AddROCIDLabels(jana::JEventLoop *loop);
		void ParseEventSize(jana::JEvent &event);
		void DataWordStats(uint32_t *iptr, uint32_t *iend, uint32_t *word_stats);

		void ParseJLabModuleData(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats);
		void Parsef250Bank(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats);
		void Parsef125Bank(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats);
		void ParseF1v2TDCBank(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats);
		void ParseF1v3TDCBank(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats);
		void ParseCAEN1190(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats);
		void ParseModuleConfiguration(uint32_t rocid, uint32_t *&iptr, uint32_t *iend, uint32_t *word_stats);

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_DAQ_online_

