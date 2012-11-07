//
//    File: JEventProcessor_danarest.h
// Created: Mon Jul 1 09:08:37 EDT 2012
// Creator: Richard Jones
//

#ifndef _JEventProcessor_danarest_
#define _JEventProcessor_danarest_


#include <string>
using namespace std;

#include <HDDM/hddm_r.hpp>

#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>

class JEventProcessor_danarest : public jana::JEventProcessor
{
 public:
   JEventProcessor_danarest();
   ~JEventProcessor_danarest();

   jerror_t init(void); ///< Called once at program start.
   jerror_t brun(JEventLoop *loop, int runnumber); ///< Called everytime a new run number is detected.
   jerror_t evnt(JEventLoop *loop, int eventnumber); ///< Called every event.
   jerror_t erun(void); ///< Called everytime run number changes, provided brun has been called.
   jerror_t fini(void);	///< Called after last event of last event source has been processed.

 private:
   ofstream *ofs;
   hddm_r::ostream *fout;
   unsigned long int Nevents_written;
};


#endif // _JEventProcessor_danarest_
