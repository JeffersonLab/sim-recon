// $Id$
//
//    File: DFCALHit_factory.h
// Created: Tue Aug  6 12:23:43 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _DFCALHit_factory_
#define _DFCALHit_factory_

#include <vector>
using namespace std;

#include <JANA/JFactory.h>
#include "TTab/DTranslationTable.h"
#include "DFCALDigiHit.h"
#include "DFCALHit.h"

// store constants so that they can be accessed by row/column number
typedef  vector< vector<double> >  fcal_digi_constants_t;

class DFCALHit_factory:public jana::JFactory<DFCALHit>{
	public:
		DFCALHit_factory(){};
		~DFCALHit_factory(){};

		static const int FCAL_MAX_CHANNELS = 2800;  // number of active channels expected

		enum fcal_quality_state {
		    GOOD,
		    BAD,
		    NOISY
		};

		// overall scale factors
		double a_scale;
		double t_scale;
                double t_min;

		// calibration constants stored in row, column format
		fcal_digi_constants_t gains;
		fcal_digi_constants_t pedestals;
		fcal_digi_constants_t time_offsets;
		fcal_digi_constants_t block_qualities;

		const double GetConstant(const fcal_digi_constants_t &the_table,
					 const int in_row, const int in_column) const;
		const double GetConstant(const fcal_digi_constants_t &the_table,
					 const DFCALDigiHit *the_digihit) const;
		const double GetConstant(const fcal_digi_constants_t &the_table,
					 const DFCALHit *the_hit) const;
		//const double GetConstant(const fcal_digi_constants_t &the_table,
		//			 const DTranslationTable *ttab,
		//			 const int in_rocid, const int in_slot, const int in_channel) const;
	private:
		jerror_t init(void);						///< Called once at program start.2
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		void FillCalibTable( fcal_digi_constants_t &table, 
				     const vector<double> &raw_table, 
				     const DFCALGeometry &fcalGeom);

};

#endif // _DFCALHit_factory_

