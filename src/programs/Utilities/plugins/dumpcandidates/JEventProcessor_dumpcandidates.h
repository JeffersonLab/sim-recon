// $Id$
//
//    File: JEventProcessor_dumpcandidates.h
// Created: Tue Feb  4 09:29:35 EST 2014
// Creator: davidl (on Linux ifarm1102 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_dumpcandidates_
#define _JEventProcessor_dumpcandidates_

#include <map>
#include <fstream>
using namespace std;

#include <JANA/JEventProcessor.h>
#include <HDGEOMETRY/DGeometry.h>

class JEventProcessor_dumpcandidates:public jana::JEventProcessor{
	public:
		JEventProcessor_dumpcandidates();
		~JEventProcessor_dumpcandidates();
		const char* className(void){return "JEventProcessor_dumpcandidates";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.


		DGeometry *dgeom;
		map<unsigned long, int> wireID;
		ofstream *ofs;
		
		int GetWireIndex(const void *wire){ return wireID[(unsigned long)wire]; }
};

#endif // _JEventProcessor_dumpcandidates_

