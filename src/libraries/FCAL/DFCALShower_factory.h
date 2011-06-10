// $Id: DFCALShower_factory.h 1899 2006-07-13 16:29:56Z davidl $
//
//    File: DFCALShower_factory.h
// Created: Tue May 17 11:57:50 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#ifndef _DFCALShower_factory_
#define _DFCALShower_factory_

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>
#include <FCAL/DFCALShower.h>
#include <FCAL/DFCALCluster.h>
#include <PID/DVertex.h>



class DFCALShower_factory:public JFactory<DFCALShower>{
	public:
		DFCALShower_factory();
		~DFCALShower_factory(){};
	
	private:
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Invoked via JEventProcessor virtual method
	//	jerror_t brun(JEventLoop *loop, int runnumber);

		DFCALShower* makeFcalShower( const DFCALCluster* cluster );
		void GetCorrectedEnergyAndPosition(const DFCALCluster* cluster, double &Ecorrected, DVector3 &pos_corrected, double &errZ, const DVector3 *aVertex);

                double m_zTarget;
		double SHOWER_ENERGY_THRESHOLD;

		double NON_LIN_COEF_A1;
		double NON_LIN_COEF_B1;
		double NON_LIN_COEF_C1;
		double NON_LIN_COEF_alfa1;

		
		double NON_LIN_COEF_A2;
		double NON_LIN_COEF_B2;
		double NON_LIN_COEF_C2;
		double NON_LIN_COEF_alfa2;

		double BUFFER_RADIUS;
		double RHG_RADIUS;

		double FCAL_RADIATION_LENGTH;
		double FCAL_CRITICAL_ENERGY;
		double FCAL_SHOWER_OFFSET;

		// Calibration constants
		// merging clusters, if any, should be done after matching sharged tracks
		//double MIN_CLUSTER_SEPARATION; // minimum separation between 2 clusters for them NOT to be merged
};


#endif // _DFCALShower_factory_

