// $Id$
//
//    File: JEventProcessor_TAGH_timewalk.h
// Created: Fri Nov 13 10:13:02 EST 2015
// Creator: nsparks (on Linux cua2.jlab.org 3.10.0-229.20.1.el7.x86_64 x86_64)
//

#ifndef _JEventProcessor_TAGH_timewalk_
#define _JEventProcessor_TAGH_timewalk_

#include <JANA/JEventProcessor.h>
#include <RF/DRFTime_factory.h>

class JEventProcessor_TAGH_timewalk:public jana::JEventProcessor{
public:
    JEventProcessor_TAGH_timewalk();
    ~JEventProcessor_TAGH_timewalk();
    const char* className(void){return "JEventProcessor_TAGH_timewalk";}

private:
    DRFTime_factory *dRFTimeFactory;
    jerror_t init(void);                                        ///< Called once at program start.
    jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);  ///< Called everytime a new run number is detected.
    jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);///< Called every event.
    jerror_t erun(void);                                        ///< Called everytime run number changes, provided brun has been called.
    jerror_t fini(void);                                        ///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_TAGH_timewalk_

