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

#include <DMatrix.h>
#include <TDecompLU.h>

#include <algorithm>
#include <map>
#include <cmath>

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

		jerror_t KalmanFilter(vector <const DFDCPoint*>points);
		jerror_t CorrectPointY(float z0,DMatrix S,
				       const DFDCPoint *point);
							
		const string toString(void);

	protected:
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
};

#endif // DFACTORY_DFDCPOINT_H

