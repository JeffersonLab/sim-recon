// $Id$
//
//    File: JEventProcessor_OmegaSkim.h
// Created: Wed May 24 13:46:12 EDT 2017
// Creator: mashephe (on Linux stanley.physics.indiana.edu 2.6.32-642.6.2.el6.x86_64 unknown)
//

#ifndef _JEventProcessor_OmegaSkim_
#define _JEventProcessor_OmegaSkim_

#include <JANA/JEventProcessor.h>

class JEventProcessor_OmegaSkim:public jana::JEventProcessor{
public:
		JEventProcessor_OmegaSkim();
		~JEventProcessor_OmegaSkim();
		const char* className(void){return "JEventProcessor_OmegaSkim";}
  
private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
  
  int WRITE_EVIO_FILE;
  int WRITE_ROOT_TREE;
};

#endif // _JEventProcessor_OmegaSkim_

