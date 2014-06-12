// $Id$
//
//    File: DBCALHit_factory.h
// Created: Tue Aug  6 09:26:13 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _DBCALHit_factory_
#define _DBCALHit_factory_

#include <vector>
#include <map>
#include <utility>
using namespace std;

#include <JANA/JFactory.h>
#include "DBCALHit.h"

typedef pair<double,double> cell_calib_t;

class DBCALHit_factory:public jana::JFactory<DBCALHit>{
	public:
		DBCALHit_factory(){};
		~DBCALHit_factory(){};

		// overall scale factors
		double a_scale;
		double t_scale;
		
		// constants indexed by cell number (module/layer/sector)
		// each cell has an up- and down-stream channel
		// it might be better to index these via a class?
		map<int,cell_calib_t> gains;
		map<int,cell_calib_t> pedestals;
		map<int,cell_calib_t> time_offsets;

		
	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		double GetConstant( map<int,cell_calib_t> &the_table, 
				    const DBCALDigiHit *the_digihit);
		void FillCalibTable( map<int,cell_calib_t> &table, 
				     const vector<double> &raw_table);
		    
};

#endif // _DBCALHit_factory_

