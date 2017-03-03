// $Id$
//
//    File: JEventProcessor_EVIO_TO_HDDM.h
// Created: Mon Feb 13 14:40:43 EST 2017
// Creator: tbritton (on Linux ifarm1101 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_EVIO_TO_HDDM_
#define _JEventProcessor_EVIO_TO_HDDM_

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <JANA/JEventProcessor.h>
#include <HDDM/hddm_s.hpp>
//#include "../hddm_s.hpp"

#include <CDC/DCDCHit.h>
#include <TOF/DTOFHit.h>
#include <FCAL/DFCALHit.h>
#include <BCAL/DBCALDigiHit.h>
#include <BCAL/DBCALTDCDigiHit.h>
#include <START_COUNTER/DSCHit.h>
#include <FDC/DFDCHit.h>
#include <PAIR_SPECTROMETER/DPSHit.h>
#include <PAIR_SPECTROMETER/DPSCHit.h>
#include <TAGGER/DTAGHHit.h>
#include <TAGGER/DTAGMHit.h>
#include <TPOL/DTPOLHit.h>

class JEventProcessor_EVIO_TO_HDDM:public jana::JEventProcessor{
	public:
		JEventProcessor_EVIO_TO_HDDM();
		~JEventProcessor_EVIO_TO_HDDM();
		const char* className(void){return "JEventProcessor_EVIO_TO_HDDM";}

		hddm_s::ostream *fout;
		std::ofstream* ofs;
		int RunNumber;

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_EVIO_TO_HDDM_

