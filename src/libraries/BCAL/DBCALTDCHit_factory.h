// $Id$
//
//    File: DBCALTDCHit_factory.h
// Created: Tue Aug  6 11:04:11 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _DBCALTDCHit_factory_
#define _DBCALTDCHit_factory_

#include <JANA/JFactory.h>
#include "DBCALTDCHit.h"

typedef pair<double,double> cell_calib_t;

class DBCALTDCHit_factory:public jana::JFactory<DBCALTDCHit>{
	public:
		DBCALTDCHit_factory(){};
		~DBCALTDCHit_factory(){};

		// overall scale factors
		double t_scale;

		map<int,cell_calib_t> time_offsets;

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		double GetConstant( map<int,cell_calib_t> &the_table, 
				    const DBCALTDCDigiHit *the_digihit);
		void FillCalibTable( map<int,cell_calib_t> &table, 
				     const vector<double> &raw_table);
};

#endif // _DBCALTDCHit_factory_

