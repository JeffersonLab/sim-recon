// $Id$
//
//    File: JEventProcessor_dumpthrowns.h
// Created: Tue Feb  4 09:29:35 EST 2014
// Creator: davidl (on Linux ifarm1102 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_dumpthrowns_
#define _JEventProcessor_dumpthrowns_

#include <map>
#include <fstream>
using namespace std;

#include <JANA/JEventProcessor.h>
#include <HDGEOMETRY/DGeometry.h>

class JEventProcessor_dumpthrowns:public jana::JEventProcessor{
	public:
		JEventProcessor_dumpthrowns();
		~JEventProcessor_dumpthrowns();
		const char* className(void){return "JEventProcessor_dumpthrowns";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.


		DGeometry *dgeom;
		map<unsigned long, int> wireID;
		ofstream *ofs;
		unsigned int MAX_CANDIDATE_FILTER;
		unsigned long events_written;
		unsigned long events_discarded;
		
		unsigned long GetCDCWireID(const DCDCWire* w){ return w->ring*100 + w->straw;}
		unsigned long GetFDCWireID(const DFDCWire* w){ return 100000 + w->layer*100 + w->wire;}
		int GetWireIndex(unsigned long id){ return wireID[id]; }
};

#endif // _JEventProcessor_dumpthrowns_

