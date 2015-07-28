#ifndef _DTPOLHit_factory_
#define _DTPOLHit_factory_

#include <JANA/JFactory.h>
#include "TTAB/DTranslationTable.h"
#include "DTPOLHit.h"


class DTPOLHit_factory:public jana::JFactory<DTPOLHit>{
	public:
		DTPOLHit_factory(){};
		~DTPOLHit_factory(){};

		// Geometric information		
		static const int NSECTORS   = 32;
		static const int NRINGS     = 24;

		static const double SECTOR_DIVISION;  
		static const double INNER_RADIUS;      // From "ACTIVE INNER DIAMETER" in catalog
		static const double OUTER_RADIUS ;     // From "ACTIVE OUTER DIAMETER" in catalog
		static const double RING_DIVISION ;    
		// (OUTER_RADIUS - INNER_RADIUS) / DTPOLRingDigiHit::NRINGS;
		// 1mm, agrees with "JUNCTION ELEMENT SEPARATION" in catalog


		// overall scale factors
		double a_scale;
		double t_scale;
                double t_base;

		// calibration constants stored by channel
		vector<double>  a_gains;
		vector<double>  a_pedestals;
		vector<double>  adc_time_offsets;

		double HIT_TIME_WINDOW;
		double ADC_THRESHOLD;

		DTPOLHit* FindMatch(int sector, double T);

		const double GetConstant(const vector<double>  &the_table,
					 const int in_sector) const;
		const double GetConstant(const vector<double>  &the_table,
					 const DSCDigiHit *the_digihit) const;
		const double GetConstant(const vector<double>  &the_table,
					 const DTPOLHit *the_hit) const;

	private:
		jerror_t init(void);						// Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	// Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	// Called every event.
		jerror_t erun(void);						// Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						// Called after last event of last event source has been processed.


};

#endif // _DTPOLHit_factory_
