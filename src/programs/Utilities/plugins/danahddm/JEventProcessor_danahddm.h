// $Id$
//
//    File: JEventProcessor_danahddm.h
// Created: Mon Mar 15 09:08:37 EDT 2010
// Creator: wolin (on Linux stan.jlab.org 2.6.18-164.el5 x86_64)
//

#ifndef _JEventProcessor_danahddm_
#define _JEventProcessor_danahddm_

#include <string>
#include <fstream>
using namespace std;

#include <HDDM/hddm_s.hpp>

#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>
using namespace jana;



class JEventProcessor_danahddm : public JEventProcessor {

   public:

      JEventProcessor_danahddm();
      ~JEventProcessor_danahddm();

      jerror_t init(void);                                 ///< Called once at program start.
      jerror_t brun(JEventLoop *loop, int runnumber);      ///< Called everytime a new run number is detected.
      jerror_t evnt(JEventLoop *loop, int eventnumber);    ///< Called every event.
      jerror_t erun(void);                                 ///< Called everytime run number changes, provided brun has been called.
      jerror_t fini(void);                                 ///< Called after last event of last event source has been processed.


   private:

      std::ofstream *file;
      hddm_s::ostream *fout;
      unsigned long Nevents_written;

      bool HDDM_USE_COMPRESSION;
      bool HDDM_USE_INTEGRITY_CHECKS;

      void Add_DTrackTimeBased(JEventLoop *loop, 
                               hddm_s::ReconViewList::iterator riter);
      
      string DMatrixDSymToString(const DMatrixDSym &mat);
};


#endif // _JEventProcessor_danahddm_
