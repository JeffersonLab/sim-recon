// $Id$
//
//    File: DBCALTDCHit_factory.h
// Created: Tue Aug  6 11:04:11 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _DBCALTDCHit_factory_
#define _DBCALTDCHit_factory_

#include <JANA/JFactory.h>
#include "TTAB/DTTabUtilities.h"
#include "DBCALTDCHit.h"

// store up/down-stream values for each detector cell
typedef pair<double,double> cell_calib_t;  
typedef vector<cell_calib_t>  bcal_digi_constants_t; 

class DBCALTDCHit_factory:public jana::JFactory<DBCALTDCHit>{
	public:
		DBCALTDCHit_factory(){};
		~DBCALTDCHit_factory(){};

		// shortcut geometry factors
		// these should really be taken from
		// DBCALGeometry/DGeometry objects
		static const int BCAL_NUM_MODULES  = 48;
		static const int BCAL_NUM_TDC_LAYERS =  3;
		static const int BCAL_NUM_SECTORS    =  4;
		static const int BCAL_MAX_TDC_CHANNELS =  1152;

		// overall scale factors
		double t_scale;
                double t_base;
		int t_rollover;

		bcal_digi_constants_t time_offsets;
        bcal_digi_constants_t channel_global_offset;
        bcal_digi_constants_t tdiff_u_d;

		const int GetCalibIndex( int module, int layer, int sector ) const {
			return BCAL_NUM_TDC_LAYERS*BCAL_NUM_SECTORS*(module-1) + BCAL_NUM_SECTORS*(layer-1) + (sector-1);
		}

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

		void FillCalibTable( bcal_digi_constants_t &table, 
				     const vector<double> &raw_table);
        void FillCalibTableShort( bcal_digi_constants_t &table,
                     const vector<double> &raw_table);
};

#endif // _DBCALTDCHit_factory_

