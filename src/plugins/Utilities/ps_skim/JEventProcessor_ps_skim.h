// $Id$
//
//    File: JEventProcessor_ps_skim.h
// Created: Mon May 18 09:52:08 EDT 2015
// Creator: nsparks (on Linux cua2.jlab.org 2.6.32-431.5.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_ps_skim_
#define _JEventProcessor_ps_skim_

#include <JANA/JEventProcessor.h>

class JEventProcessor_ps_skim:public jana::JEventProcessor{
public:
    JEventProcessor_ps_skim();
    ~JEventProcessor_ps_skim();
    const char* className(void){return "JEventProcessor_ps_skim";}

private:
    jerror_t init(void); ///< Called once at program start.
    jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber); ///< Called everytime a new run number is detected.
    jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber); ///< Called every event.
    jerror_t erun(void); ///< Called everytime run number changes, provided brun has been called.
    jerror_t fini(void); ///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_ps_skim_

