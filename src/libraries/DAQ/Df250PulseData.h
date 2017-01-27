// $Id$
// $HeadURL$
//
//    File: Df250PulseData.h
// Created: Thu Sep  1 09:30:13 EDT 2016
// Creator: davidl (on Darwin harriet.jlab.org 15.6.0 x86_64)
//
// This data type corresponds to a new firware version deployed
// in Fall 2016. It replaces the f250PulseIntegral, Df250PulsePedestal,
// and Df250PulseTime objects. Those are left in place to support the
// data taken with the earlier firmware.
//
// This corresponds to data type 9 in the new firmware (previously,
// data type 9 was reserved for "streaming raw data", but that was
// never implemented in the firware).
//
// Note that the new data format ties the integral, pedestal, and
// time together (unlike the original). This allows us to put all
// of these in a single object without requiring additional objects
// and object associations.
//


#ifndef _Df250PulseData_
#define _Df250PulseData_

#include <DAQ/DDAQAddress.h>

class Df250PulseData:public DDAQAddress{
	
	/// Holds pulse integral,pedestal,time data for one identified
	/// pulse in one event in one channel of a single
	/// f250 Flash ADC module.
	
	public:
		JOBJECT_PUBLIC(Df250PulseData);

		Df250PulseData(uint32_t rocid=0, uint32_t slot=0, uint32_t channel=0, uint32_t itrigger=0 
			, uint32_t event_within_block=0
			, bool     QF_pedestal=false
			, uint32_t pedestal=0
			, uint32_t integral=0
			, bool     QF_NSA_beyond_PTW=false
			, bool     QF_overflow=false
			, bool     QF_underflow=false
			, uint32_t nsamples_over_threshold=0
			, uint32_t course_time=0
			, uint32_t fine_time=0
			, uint32_t pulse_peak=0
			, bool     QF_vpeak_beyond_NSA=false
			, bool     QF_vpeak_not_found=false
			, bool     QF_bad_pedestal=false
			, uint32_t pulse_number=0
			, uint32_t nsamples_integral=0
			, uint32_t nsamples_pedestal=0
			, bool     emulated=false
        	, uint32_t integral_emulated=0
        	, uint32_t pedestal_emulated=0
			, uint32_t time_emulated=0
			, uint32_t pulse_peak_emulated=0)
			:DDAQAddress(rocid, slot, channel, itrigger)
			, event_within_block(event_within_block)
			, QF_pedestal(QF_pedestal)
			, pedestal(pedestal)
			, integral(integral)
			, QF_NSA_beyond_PTW(QF_NSA_beyond_PTW)
			, QF_overflow(QF_overflow)
			, QF_underflow(QF_underflow)
			, nsamples_over_threshold(nsamples_over_threshold)
			, course_time(course_time)
			, fine_time(fine_time)
			, pulse_peak(pulse_peak)
			, QF_vpeak_beyond_NSA(QF_vpeak_beyond_NSA)
			, QF_vpeak_not_found(QF_vpeak_not_found)
			, QF_bad_pedestal(QF_bad_pedestal)
			, pulse_number(pulse_number)
			, nsamples_integral(nsamples_integral)
			, nsamples_pedestal(nsamples_pedestal)
			, emulated(emulated)
			, integral_emulated(integral_emulated)
			, pedestal_emulated(pedestal_emulated)
			, time_emulated(time_emulated)
			, pulse_peak_emulated(pulse_peak_emulated){}

		
		// from word 1
		uint32_t event_within_block;
		bool     QF_pedestal;
		uint32_t pedestal;
		
		// from word 2
		uint32_t integral;
		bool     QF_NSA_beyond_PTW;
		bool     QF_overflow;
		bool     QF_underflow;
		uint32_t nsamples_over_threshold;
		
		// from word 3
		uint32_t course_time;               //< 4 ns/count
		uint32_t fine_time;                 //< 0.0625 ns/count
		uint32_t pulse_peak;
		bool     QF_vpeak_beyond_NSA;
		bool     QF_vpeak_not_found;
		bool     QF_bad_pedestal;

		// supplemental (not from data stream)
		uint32_t pulse_number;         ///< pulse number for this channel, this event starting from 0		
		uint32_t nsamples_integral;    ///< number of samples used in integral 
		uint32_t nsamples_pedestal;    ///< number of samples used in pedestal
		bool     emulated;             ///< true if made from Window Raw Data
        uint32_t integral_emulated;    ///< Value calculated from raw data (if available)
        uint32_t pedestal_emulated;    ///< Value calculated from raw data (if available)
        uint32_t time_emulated;        ///< Value calculated from raw data (if available)
        uint32_t course_time_emulated;        ///< Value calculated from raw data (if available) - debug
        uint32_t fine_time_emulated;        ///< Value calculated from raw data (if available) - debug
        uint32_t pulse_peak_emulated;  ///< Value calculated from raw data (if available)
        uint32_t QF_emulated;

		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
		
			uint32_t QF = 0; // make single quality factor number for compactness
			if( QF_pedestal         ) QF |= (1<<0);
			if( QF_NSA_beyond_PTW   ) QF |= (1<<1);
			if( QF_overflow         ) QF |= (1<<2);
			if( QF_underflow        ) QF |= (1<<3);
			if( QF_vpeak_beyond_NSA ) QF |= (1<<4);
			if( QF_vpeak_not_found  ) QF |= (1<<5);
			if( QF_bad_pedestal     ) QF |= (1<<6);
			
			uint32_t emulated_all = 0; // make single emulated flags code
			if( emulated          ) emulated_all |= (1<<0);
			if( integral_emulated ) emulated_all |= (1<<1);
			if( pedestal_emulated ) emulated_all |= (1<<2);
		
			DDAQAddress::toStrings(items);
			//AddString(items, "event_within_block"      , "%d", event_within_block      );
			AddString(items, "integral"                , "%d", integral                );
			AddString(items, "pedestal"                , "%d", pedestal                );
			AddString(items, "course_time"             , "%d", course_time             );
			AddString(items, "fine_time"               , "%d", fine_time               );
			AddString(items, "pulse_peak"              , "%d", pulse_peak              );
			AddString(items, "pulse_number"            , "%d", pulse_number            );
			AddString(items, "nsamples_integral"       , "%d", nsamples_integral       );
			AddString(items, "nsamples_pedestal"       , "%d", nsamples_pedestal       );
			AddString(items, "nsamples_over_threshold" , "%d", nsamples_over_threshold );
			AddString(items, "QF"                      , "%x", QF                      );
			AddString(items, "emulated"                , "%x", emulated_all            );

			AddString(items, "integral_emulated"       , "%d", integral_emulated       );
			AddString(items, "pedestal_emulated"       , "%d", pedestal_emulated       );
			AddString(items, "course_time_emulated"           , "%d", course_time_emulated           );
			AddString(items, "fine_time_emulated"           , "%d", fine_time_emulated           );
			AddString(items, "pulse_peak_emulated"     , "%d", pulse_peak_emulated     );
			AddString(items, "QF_emulated"     , "%x", QF_emulated     ); 
		}
};

#endif // _Df250PulseData_

