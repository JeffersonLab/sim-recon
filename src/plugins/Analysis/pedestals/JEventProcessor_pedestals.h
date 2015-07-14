// $Id$
//
//    File: JEventProcessor_pedestals.h
// Created: Fri Jun 20 22:13:58 EDT 2014
// Creator: davidl (on Darwin harriet.local 13.2.0 i386)
//

#ifndef _JEventProcessor_pedestals_
#define _JEventProcessor_pedestals_

#include <TH2.h>

#include <JANA/JEventProcessor.h>
#include <DAQ/DDAQAddress.h>

class JEventProcessor_pedestals:public jana::JEventProcessor{
	public:
		JEventProcessor_pedestals();
		~JEventProcessor_pedestals();
		const char* className(void){return "JEventProcessor_pedestals";}

		class csc_t{
			public:
				csc_t(uint32_t rocid, uint32_t slot, uint32_t channel):rocid(rocid), slot(slot), channel(channel){}
				uint32_t rocid;
				uint32_t slot;
				uint32_t channel;
		};
		
		TH2D* GetHist(const DDAQAddress *hit);
		
		map<csc_t, TH2D*> all_hists;
		//pthread_mutex_t mutex;

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

};

// Need to define a less-than operator for csc_t so it can be used as a key in the map
inline bool operator<(const JEventProcessor_pedestals::csc_t &a, const JEventProcessor_pedestals::csc_t &b){
	if(a.rocid < b.rocid) return true;
	if(a.rocid > b.rocid) return false;
	if(a.slot < b.slot) return true;
	if(a.slot > b.slot) return false;
	if(a.channel < b.channel) return true;
	return false;
}


#endif // _JEventProcessor_pedestals_

