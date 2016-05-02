// $Id$
//
//    File: JEventProcessor_es_test.cc
//

#include "JEventProcessor_es_test.h"
using namespace jana;

#include "../eventstore/DESSkimData.h"

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
    void InitPlugin(JApplication *app){
        InitJANAPlugin(app);
        app->AddProcessor(new JEventProcessor_es_test());
    }
} // "C"


//------------------
// JEventProcessor_es_test (Constructor)
//------------------
JEventProcessor_es_test::JEventProcessor_es_test()
{

}

//------------------
// ~JEventProcessor_es_test (Destructor)
//------------------
JEventProcessor_es_test::~JEventProcessor_es_test()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_es_test::init(void)
{

    return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_es_test::brun(JEventLoop *eventLoop, int32_t runnumber)
{
    // This is called whenever the run number changes
    return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_es_test::evnt(JEventLoop *loop, uint64_t eventnumber)
{
  //const DESSkimData *es_data = NULL;
  //loop->GetSingle(es_data);
  vector<const DESSkimData *> es_data;
  loop->Get(es_data);

  cout << "Event " << eventnumber << endl;
  es_data[0]->Print();

  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_es_test::erun(void)
{
    // This is called whenever the run number changes, before it is
    // changed to give you a chance to clean up before processing
    // events from the next run number.
    return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_es_test::fini(void)
{
    // Called before program exit after event processing is finished.
    return NOERROR;
}

