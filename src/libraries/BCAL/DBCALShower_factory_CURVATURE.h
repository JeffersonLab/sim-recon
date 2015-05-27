// $Id$
//
//    File: DBCALShower_factory_CURVATURE.h
// Created: Fri Mar 27 10:57:45 CST 2015
// Creator: beattite (on Linux eos.phys.uregina.ca 2.6.32-504.12.2.el6.x86_64 x86_64)
//

#ifndef _DBCALShower_factory_CURVATURE_
#define _DBCALShower_factory_CURVATURE_

#include <JANA/JFactory.h>
#include "DBCALShower.h"

class DBCALShower_factory_CURVATURE:public jana::JFactory<DBCALShower>{
	public:
		DBCALShower_factory_CURVATURE(){};
		~DBCALShower_factory_CURVATURE(){};
		const char* Tag(void){return "CURVATURE";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

  double m_zTarget;

// energy calibration parameters
  float m_scaleZ_p0;
  float m_scaleZ_p1;
  float m_scaleZ_p2;
  float m_scaleZ_p3;
  float m_nonlinZ_p0;
  float m_nonlinZ_p1;
  float m_nonlinZ_p2;
  float m_nonlinZ_p3;

  double position[2][4][12][32], sigma[2][4][12][32], temptheta, tempenergy, PHITHRESHOLD, ZTHRESHOLD, TTHRESHOLD, ETHRESHOLD;
  double posoffset1, posoffset2, sigoffset1, sigoffset2, zbinposition[2][4], zbinwidth[2][4], dataposition, datasigma;
  int i, j, k, l, bin, k2, l2, layer, angle, energy;
  vector<int> overlap;
  vector<double> recon_showers_phi;
  vector<double> recon_showers_theta;
  vector<double> recon_showers_E;
  vector<double> recon_showers_t;

};

#endif // _DBCALShower_factory_CURVATURE_

