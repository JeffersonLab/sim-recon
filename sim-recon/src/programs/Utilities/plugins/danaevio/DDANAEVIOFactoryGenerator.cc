// Author: Elliott Wolin 19-Mar-2010


#include <JANA/JApplication.h>

#include "DDANAEVIOFactoryGenerator.h"
#include "DDANAEVIO_factory.h"
#include "JEventProcessor_danaevio.h"


jerror_t DDANAEVIOFactoryGenerator::GenerateFactories(JEventLoop *loop) {
  loop->AddFactory(new DDANAEVIO_factory());
  return NOERROR;
}



//----------------------------------------------------------------------------


// for initializing plugins
extern "C" {

  void InitPlugin(JApplication *app) {

    // initialize plugin system
    InitJANAPlugin(app);


    // add DANAEVIO factory
    app->AddFactoryGenerator(new DDANAEVIOFactoryGenerator());


    // Add DANAEVIO event processor to write out events
    // default for writing is true, the normal case if factory is included
    // to turn off specify -PEVIO::WRITEOUT=0
    bool evioWriteOut =true;
    gPARMS->SetDefaultParameter("EVIO:WRITEOUT",evioWriteOut);
    if(evioWriteOut) {
      app->AddProcessor(new JEventProcessor_danaevio(),true);
    } else {
      jout << endl << endl << "    *** No EVIO output file will be generated ***" << endl << endl << endl;
    }

  }

} // "extern C"


//----------------------------------------------------------------------------
