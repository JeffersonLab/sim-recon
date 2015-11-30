// $Id$
//
//    File: JEventProcessor_PS_timing.h
// Created: Sat Nov 21 17:21:28 EST 2015
// Creator: nsparks (on Linux cua2.jlab.org 3.10.0-327.el7.x86_64 x86_64)
//

#ifndef _JEventProcessor_PS_timing_
#define _JEventProcessor_PS_timing_

#include <JANA/JEventProcessor.h>

class JEventProcessor_PS_timing:public jana::JEventProcessor{
public:
    JEventProcessor_PS_timing();
    ~JEventProcessor_PS_timing();
    const char* className(void){return "JEventProcessor_PS_timing";}

private:
    jerror_t init(void);                                        ///< Called once at program start.
    jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);  ///< Called everytime a new run number is detected.
    jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);///< Called every event.
    jerror_t erun(void);                                        ///< Called everytime run number changes, provided brun has been called.
    jerror_t fini(void);                                        ///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_PS_timing_

