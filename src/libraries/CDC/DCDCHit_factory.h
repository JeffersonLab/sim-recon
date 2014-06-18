// $Id$
//
//    File: DCDCHit_factory.h
// Created: Tue Aug  6 11:29:56 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _DCDCHit_factory_
#define _DCDCHit_factory_

#include <vector>
using namespace std;

#include <JANA/JFactory.h>
#include "HDGEOMETRY/DGeometry.h"
#include "DCDCHit.h"

// store constants indexed by ring/straw number
typedef  vector< vector<double> >  cdc_digi_constants_t;


class DCDCHit_factory:public jana::JFactory<DCDCHit>{
	public:
		DCDCHit_factory(){};
		~DCDCHit_factory(){};

		// overall scale factors.
		double a_scale;
		double t_scale;


		cdc_digi_constants_t gains;
		cdc_digi_constants_t pedestals;
		cdc_digi_constants_t time_offsets;
		
	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		void CalcNstraws(jana::JEventLoop *eventLoop, int runnumber, vector<int> &Nstraws);
		void FillCalibTable(vector< vector<double> > &table, vector<double> &raw_table, 
				    vector<int> &Nstraws);

};

#endif // _DCDCHit_factory_

