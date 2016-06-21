// $Id: DFCALShower_factory.h 1899 2006-07-13 16:29:56Z davidl $
//
//    File: DFCALShower_factory.h
// Created: Tue May 17 11:57:50 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
// Edited: B. Schaefer 3/23/2012 removed radiation hard insert functionality

#ifndef _DFCALShower_factory_
#define _DFCALShower_factory_

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>
#include <FCAL/DFCALShower.h>
#include <FCAL/DFCALCluster.h>
#include <DANA/DApplication.h>

class DFCALShower_factory:public JFactory<DFCALShower>{
	public:
		DFCALShower_factory();
		~DFCALShower_factory(){};
	
	private:
		jerror_t evnt(JEventLoop *eventLoop, uint64_t eventnumber);	///< Invoked via JEventProcessor virtual method
		jerror_t brun(JEventLoop *loop, int32_t runnumber);

		void GetCorrectedEnergyAndPosition(const DFCALCluster* cluster, double &Ecorrected, DVector3 &pos_corrected, double &errZ, const DVector3 *aVertex);

                double m_zTarget,m_FCALfront;

		double LOAD_CCDB_CONSTANTS;
		double SHOWER_ENERGY_THRESHOLD;
		double energy_cutoff;
		double linfit_slope;
		double linfit_intercept;
		double expfit_param1;
		double expfit_param2;
		double expfit_param3;


		double FCAL_RADIATION_LENGTH;
		double FCAL_CRITICAL_ENERGY;
		double FCAL_SHOWER_OFFSET;
		double FCAL_C_EFFECTIVE;
};


#endif // _DFCALShower_factory_

