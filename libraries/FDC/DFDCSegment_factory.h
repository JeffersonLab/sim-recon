//******************************************************************
// DFDCPoint_factory.h: class definition for a factory creating
// space points from pseudopoints
//******************************************************************
#ifndef DFACTORY_DFDCPOINT_H
#define DFACTORY_DFDCPOINT_H

#include "JANA/JFactory.h"
#include "JANA/JException.h"
#include "JANA/JStreamLog.h"

#include "DFDCPoint.h"
#include "DFDCPseudo.h"
#include "DFDCHit.h"
#include "DFDCGeometry.h"

#include <TDecompLU.h>

#include <algorithm>
#include <map>
#include <cmath>

/* The folowing are for interpreting grid of Lorentz deflection data */
#define PACKAGE_Z_POINTS 7
#define LORENTZ_X_POINTS 21
#define LORENTZ_Z_POINTS 4*PACKAGE_Z_POINTS


///
/// class DFDCPoint_factory: definition for a JFactory that
/// produces space points from pseudopoints.
/// 
class DFDCPoint_factory : public JFactory<DFDCPoint> {
	public:
		
		///
		/// DFDCPoint_factory::DFDCPoint_factory():
		/// default constructor -- initializes log file
		///
		DFDCPoint_factory();
		
		///
		/// DFDCPoint_factory::~DFDCPoint_factory():
		/// default destructor -- closes log file
		///
		~DFDCPoint_factory();	

		jerror_t KalmanFilter(vector <DFDCPseudo*>points);
		jerror_t CorrectPointY(float z0,DMatrix S,DFDCPseudo *point);
							
		const string toString(void);

	protected:
		///
		/// DFDCPoint_factory::brun():
		///
		jerror_t brun(JEventLoop *eventLoop, int eventNo);

		///
		/// DFDCPoint_factory::evnt():
		/// this is the place that finds track segments and  
		/// converts pseudopoints into space points.
		///
		jerror_t evnt(JEventLoop *eventLoop, int eventNo);

	private:
		DFDCGeometry _geo;
		JStreamLog* _log;
		ofstream* logFile;
		vector<DMatrix>segments;

		// Variables for implementing lorentz effect
		// due to the magnetic field).
		float lorentz_x[LORENTZ_X_POINTS];
		float lorentz_z[LORENTZ_Z_POINTS];
		float lorentz_nx[LORENTZ_X_POINTS][LORENTZ_Z_POINTS];
		float lorentz_nz[LORENTZ_X_POINTS][LORENTZ_Z_POINTS];

};

#endif // DFACTORY_DFDCPOINT_H

