// $Id$
//
//    File: JEventProcessor_TAGH_doubles.h
// Created: Fri Apr 29 15:19:27 EDT 2016
// Creator: nsparks (on Linux cua2.jlab.org 3.10.0-327.13.1.el7.x86_64 x86_64)
//

#ifndef _JEventProcessor_TAGH_doubles_
#define _JEventProcessor_TAGH_doubles_

#include <JANA/JEventProcessor.h>

class JEventProcessor_TAGH_doubles:public jana::JEventProcessor{
    public:
        JEventProcessor_TAGH_doubles();
        ~JEventProcessor_TAGH_doubles();
        const char* className(void){return "JEventProcessor_TAGH_doubles";}

    private:
        jerror_t init(void);						///< Called once at program start.
        jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
        jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
        jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
        jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_TAGH_doubles_
