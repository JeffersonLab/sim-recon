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

// store up/down-stream values for each detector cell
typedef pair<double,double> cell_calib_t;  
// store constants indexed by cell ID (module/layer/sector)
typedef  map<int,cell_calib_t>  bcal_digi_constants_t;    // switch to a hashed map?

class DBCALTDCHit_factory:public jana::JFactory<DBCALTDCHit>{
	public:
		DBCALTDCHit_factory(){};
		~DBCALTDCHit_factory(){};

		// shortcut geometry factors
		// these should really be taken from
		// DBCALGeometry/DGeometry objects
		static const int BCAL_NUM_MODULES  = 48;
		static const int BCAL_NUM_LAYERS   =  4;
		static const int BCAL_NUM_SECTORS  =  4;
		static const int BCAL_MAX_CHANNELS =  1536;

		// overall scale factors
		double t_scale;

		bcal_digi_constants_t time_offsets;

		const double GetConstant( const bcal_digi_constants_t &the_table,
					  const int in_module, const int in_layer,
					  const int in_sector, const int in_end) const;
		const double GetConstant( const bcal_digi_constants_t &the_table,
					  const DBCALTDCDigiHit *the_digihit) const;
		const double GetConstant( const bcal_digi_constants_t &the_table,
					  const DBCALTDCHit *the_hit) const;
		//const double GetConstant( const bcal_digi_constants_t &the_table,
		//			  const DTranslationTable *ttab,
		//			  const int in_rocid, const int in_slot, const int in_channel) const;


	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		void FillCalibTable( map<int,cell_calib_t> &table, 
				     const vector<double> &raw_table);
};

#endif // _DBCALTDCHit_factory_

