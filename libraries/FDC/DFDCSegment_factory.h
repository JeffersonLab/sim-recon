//******************************************************************
// DFDCSegment_factory.h: class definition for a factory creating
// track segments from pseudopoints
//******************************************************************
#ifndef DFACTORY_DFDCSEGMENT_H
#define DFACTORY_DFDCSEGMENT_H

#include "JANA/JFactory.h"
#include "JANA/JException.h"
#include "JANA/JStreamLog.h"

#include "DFDCSegment.h"
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
/// class DFDCSegment_factory: definition for a JFactory that
/// produces space points from pseudopoints.
/// 
class DFDCSegment_factory : public JFactory<DFDCSegment> {
	public:
		
		///
		/// DFDCSegment_factory::DFDCSegment_factory():
		/// default constructor -- initializes log file
		///
		DFDCSegment_factory();
		
		///
		/// DFDCSegment_factory::~DFDCSegment_factory():
		/// default destructor -- closes log file
		///
		~DFDCSegment_factory();	

		jerror_t KalmanFilter(vector <DFDCPseudo*>points);
		jerror_t CorrectPointY(DMatrix S,DFDCPseudo *point);
		jerror_t GetProcessNoiseCovariance(DMatrix S, 
	                  vector<DFDCPseudo*>points,DMatrix &Q);
	        jerror_t GetHelicalTrackPosition(double z, DMatrix S,
	                  double &xpos,double &ypos);
		jerror_t GetTrackProjectionMatrix(double z,DMatrix S,
						  DMatrix &H);
		jerror_t GetStateTransportMatrix(double oldx, double oldy,
	                  double x,double y, DMatrix S, DMatrix &F);
		jerror_t GetStateVector(double oldx, double oldy,
	                  double old_z,double x, double y,double z,
			  DMatrix S,DMatrix &S1);

		const string toString(void);

	protected:
		///
		/// DFDCSegment_factory::brun():
		///
		jerror_t brun(JEventLoop *eventLoop, int eventNo);

		///
		/// DFDCSegment_factory::evnt():
		/// this is the place that finds track segments and  
		/// converts pseudopoints into space points.
		///
		jerror_t evnt(JEventLoop *eventLoop, int eventNo);

	private:
		DFDCGeometry _geo;
		JStreamLog* _log;
		ofstream* logFile;

	        // Average magnetic field in each package
	        double BField[4];

		// Variables for implementing lorentz effect
		// due to the magnetic field).
		double lorentz_x[LORENTZ_X_POINTS];
		double lorentz_z[LORENTZ_Z_POINTS];
		double lorentz_nx[LORENTZ_X_POINTS][LORENTZ_Z_POINTS];
		double lorentz_nz[LORENTZ_X_POINTS][LORENTZ_Z_POINTS];
	
};

#endif // DFACTORY_DFDCSEGMENT_H

