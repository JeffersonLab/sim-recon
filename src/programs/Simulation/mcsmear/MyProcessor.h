// Author: David Lawrence Sat Jan 29 09:37:37 EST 2011
//
//
// MyProcessor.h
//
/// Processor for mcsmear
///

#include <string>

#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>
using namespace jana;

#include <fstream>
#include <HDDM/hddm_s.hpp>



class MyProcessor:public JEventProcessor
{
   public:
      jerror_t init(void);                              ///< Called once at program start.
      jerror_t brun(JEventLoop *loop, int runnumber) {  ///< Called everytime a new run number is detected.
         return NOERROR;
      }
      jerror_t evnt(JEventLoop *loop, int eventnumber); ///< Called every event.
      jerror_t erun(void) {                             ///< Called everytime run number changes, provided brun has been called.
         return NOERROR;
      }
      jerror_t fini(void);                              ///< Called after last event of last event source has been processed.

      ofstream *ofs;
      hddm_s::ostream *fout; 
      unsigned long Nevents_written;

   private:
      bool HDDM_USE_COMPRESSION;
      bool HDDM_USE_INTEGRITY_CHECKS;
};
