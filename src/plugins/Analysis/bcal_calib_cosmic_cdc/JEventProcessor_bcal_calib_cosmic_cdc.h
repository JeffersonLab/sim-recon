// $Id$
//
//    File: JEventProcessor_bcal_calib_cosmic_cdc.h
// Created: Tue Jul  1 13:11:51 EDT 2014
// Creator: dalton (on Linux gluon05.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_bcal_calib_cosmic_cdc_
#define _JEventProcessor_bcal_calib_cosmic_cdc_

#include <JANA/JEventProcessor.h>
#include "TTree.h"


// Doxygen documentation
/** 
This class is used to project stright lines from cosmic rays in the CDC lines to the BCAL for the purposes of calibrating the BCAL.


*/

class JEventProcessor_bcal_calib_cosmic_cdc:public jana::JEventProcessor{
	public:
		JEventProcessor_bcal_calib_cosmic_cdc();
		~JEventProcessor_bcal_calib_cosmic_cdc();
		const char* className(void){return "JEventProcessor_bcal_calib_cosmic_cdc";}

		TTree *bcal_calib_cosmic_cdc_tree;
		int eventnum;
		int cell; 
		int tlayer;
		int tmodule;
		int tsector;
		int tglobalsect;
		int numcells; ///< Number of BCAL cells intersected by the track
		float tdist;  
		float use;  
		float dse;  
		float track_m;
  		float track_c;
		float chisq;
		int Ndof;

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		/// Command Line Parameters
		int VERBOSE;
};

#endif // _JEventProcessor_bcal_calib_cosmic_cdc_

