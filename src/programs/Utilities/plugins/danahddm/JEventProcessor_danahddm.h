// $Id$
//
//    File: JEventProcessor_danahddm.h
// Created: Mon Mar 15 09:08:37 EDT 2010
// Creator: wolin (on Linux stan.jlab.org 2.6.18-164.el5 x86_64)
//

#ifndef _JEventProcessor_danahddm_
#define _JEventProcessor_danahddm_


#include <string>
using namespace std;

#include <HDDM/hddm_s.h>

#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>
using namespace jana;



class JEventProcessor_danahddm : public JEventProcessor {

	public:

		JEventProcessor_danahddm();
		~JEventProcessor_danahddm();

		jerror_t init(void);											///< Called once at program start.
		jerror_t brun(JEventLoop *loop, int runnumber);		///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Called every event.
		jerror_t erun(void);											///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);											///< Called after last event of last event source has been processed.


	private:

		s_iostream_t *file;
		unsigned long Nevents_written;

		void JEventProcessor_danahddm::Add_DTrackTimeBased(JEventLoop *loop, s_ReconView_t *recon);
};


#endif // _JEventProcessor_danahddm_
