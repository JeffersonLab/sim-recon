// $Id$
// $HeadURL$
//
//    File: Df125FDCPulse.h
// Created: Fri Nov  13 16:16:00 EDT 2015
// Creator: davidl (on Darwin harriet.jlab.org 13.4.0 x86_64)
//

#ifndef _Df125FDCPulse_
#define _Df125FDCPulse_

#include <DAQ/DDAQAddress.h>

class Df125FDCPulse:public DDAQAddress{
	
	/// Holds pulse integral data for one identified
	/// pulse in one event in one channel of a single
	/// f125 Flash ADC module.
	
	public:
		JOBJECT_PUBLIC(Df125FDCPulse);

		Df125FDCPulse(uint32_t rocid=0, uint32_t slot=0, uint32_t channel=0, uint32_t itrigger=0
		                , uint32_t NPK=0
						, uint32_t le_time=0
						, uint32_t time_quality_bit=0
						, uint32_t overflow_count=0
						, uint32_t pedestal=0
						, uint32_t integral=0
						, uint32_t peak_amp=0
						, uint32_t peak_time=0
						, uint32_t word1=0
						, uint32_t word2=0
						, uint32_t nsamples_pedestal=1
						, uint32_t nsamples_integral=1
						, bool emulated=false )
						  :DDAQAddress(rocid, slot, channel, itrigger)
						  , NPK(NPK)
						  , le_time(le_time)
						  , time_quality_bit(time_quality_bit)
						  , overflow_count(overflow_count)
						  , pedestal(pedestal)
						  , integral(integral)
						  , peak_amp(peak_amp)
						  , peak_time(peak_time)
						  , word1(word1)
						  , word2(word2)
						  , nsamples_pedestal(nsamples_pedestal)
						  , nsamples_integral(nsamples_integral)
						  , emulated(emulated){}

		uint32_t NPK;                  ///< from first word
		uint32_t le_time;              ///< from first word
		uint32_t time_quality_bit;     ///< from first word
		uint32_t overflow_count;       ///< from first word
		uint32_t pedestal;             ///< from second word 
		uint32_t integral;             ///< from second word (type 6)
		uint32_t peak_amp;             ///< from second word (type 9)
		uint32_t peak_time;            ///< from second word
		uint32_t word1;                ///< first word
		uint32_t word2;                ///< second word
		uint32_t nsamples_pedestal;    ///< number of samples used in integral 
		uint32_t nsamples_integral;    ///< number of samples used in pedestal
		bool     emulated;             ///< true if made from Window Raw Data
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			DDAQAddress::toStrings(items);
			AddString(items, "le_time",           "%d", le_time);
			AddString(items, "integral",          "%d", integral);
			AddString(items, "pedestal",          "%d", pedestal);
			AddString(items, "NPK",               "%d", NPK);
			AddString(items, "time_quality_bit",  "%d", time_quality_bit);
			AddString(items, "overflow_count",    "%d", overflow_count);
			AddString(items, "peak_amp",          "%d", peak_amp);
			AddString(items, "peak_time",         "%d", peak_time);
			AddString(items, "nsamples_integral", "%d", nsamples_integral);
			AddString(items, "nsamples_pedestal", "%d", nsamples_pedestal);
			AddString(items, "emulated",          "%d", emulated);
		}
};

#endif // _Df125FDCPulse_

