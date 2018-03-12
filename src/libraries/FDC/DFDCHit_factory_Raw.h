// $Id$
//
//    File: DFDCHit_factory_Raw.h
// Created: Wed Aug  7 11:55:02 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _DFDCHit_factory_Raw_
#define _DFDCHit_factory_Raw_

#include <vector>
using namespace std;

#include <JANA/JFactory.h>
#include "TTAB/DTTabUtilities.h"

#include "DFDCHit.h"
#include "DFDCWireDigiHit.h"
#include "DFDCCathodeDigiHit.h"

// store constants indexed by gPlane/element number
typedef  vector< vector<double> >  fdc_digi_constants_t;

class DFDCHit_factory_Raw:public jana::JFactory<DFDCHit>{
	public:
		DFDCHit_factory_Raw(){};
		~DFDCHit_factory_Raw(){};
        const char* Tag(void){return "Raw";}

		// overall scale factors
		double a_scale;
		double t_scale;
		double t_base;
		double fadc_t_base;

		// calibration constant tables
		fdc_digi_constants_t a_gains;
		fdc_digi_constants_t a_pedestals;
		fdc_digi_constants_t timing_offsets;

		const double GetConstant(const fdc_digi_constants_t &the_table,
					 const int in_gPlane, const int in_element) const;
		const double GetConstant(const fdc_digi_constants_t &the_table,
					 const DFDCCathodeDigiHit *the_digihit) const;
		const double GetConstant(const fdc_digi_constants_t &the_table,
					 const DFDCWireDigiHit *the_digihit) const;
		const double GetConstant(const fdc_digi_constants_t &the_table,
					 const DFDCHit *the_hit) const;

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		void LoadPackageCalibTables(jana::JEventLoop *eventLoop, string ccdb_prefix);
};

#endif // _DFDCHit_factory_Raw_

