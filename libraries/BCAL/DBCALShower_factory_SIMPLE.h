// $Id$
//
//    File: DBCALShower_factory_SIMPLE.h
// Created: Tue Mar 20 14:44:12 EDT 2007
// Creator: davidl (on Darwin swire-b241.jlab.org 8.9.0 powerpc)
//

#ifndef _DBCALShower_factory_SIMPLE_
#define _DBCALShower_factory_SIMPLE_

#include <TH1.h>
#include <TH3.h>

#include <JANA/JFactory.h>
#include "DBCALShower.h"
#include "DBCALHit.h"

class DBCALShower_factory_SIMPLE:public JFactory<DBCALShower>{
	public:
		DBCALShower_factory_SIMPLE();
		~DBCALShower_factory_SIMPLE(){};
		const string toString(void);
		const char* Tag(void){return "SIMPLE";}
		
		class pshower_t{
			public:
				const DBCALHit* upstream_hit;
				const DBCALHit* downstream_hit;
				double E;
				double t;
				double z;
				double phi;
				double R;
				bool used;
		};
		
		void ModuleLayerSectorToPhiR(int module, int layer, int sector, double &phi, double &R);
		double FindZFromSingleHit(const DBCALHit *hit);

	private:
		//jerror_t init(void);						///< Called once at program start.
		jerror_t brun(JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		//jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		//jerror_t fini(void);						///< Called after last event of last event source has been processed.
		
		double UP_DOWN_COINCIDENCE_WINDOW;
		double ENERGY_SCALE_FACTOR;
		double SIGNAL_VELOCITY;
		double ATTENUATION_LENGTH;
		double Z_CENTER;
		double MIN_CLUSTER_SPACING;
		double MIN_SHOWER_ENERGY;
		bool DEBUG_HISTS;
		
		TH1F *hit_element_dist;
		TH1F *hit_element_dist_z;
		TH1F *hit_element_dist_phi;
		TH1F *hit_element_dist_r;
		TH3F *x_vs_y_vs_z;
		TH1F *r_shower, *phi_shower, *z_shower;
		TH1F *r_element, *phi_element, *z_element;
		TH1F *E_element, *E_shower;
		TH2F *E_upstream_vs_downstream;
		TH1F *E_clustersum, *E_elementsum;
};

#endif // _DBCALShower_factory_SIMPLE_

