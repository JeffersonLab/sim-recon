// $Id$
//
//    File: JEventProcessor_es_test.h
// Created: Mon May 18 09:52:08 EDT 2015
//

#ifndef _JEventProcessor_es_test_
#define _JEventProcessor_es_test_

#include <JANA/JEventProcessor.h>

class JEventProcessor_es_test:public jana::JEventProcessor{
public:
    JEventProcessor_es_test();
    ~JEventProcessor_es_test();
    const char* className(void){return "JEventProcessor_es_test";}

private:
    jerror_t init(void); ///< Called once at program start.
    jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber); ///< Called everytime a new run number is detected.
    jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber); ///< Called every event.
    jerror_t erun(void); ///< Called everytime run number changes, provided brun has been called.
    jerror_t fini(void); ///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_es_test_

