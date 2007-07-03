//*****************************************************************************
// DTrackLinker_factory.h: class definition for a factory linking segments
//*****************************************************************************
#ifndef DFACTORY_DTRACKLINKER_H
#define DFACTORY_DTRACKLINKER_H

#include "JANA/JFactory.h"
#include "JANA/JException.h"
#include "JANA/JStreamLog.h"
#include <DMatrix.h>

#include "DTrackLinker.h"
#include "FDC/DFDCSegment.h"

#include "TRACKING/DMagneticFieldMap.h"
#include <TDecompLU.h>

#include <algorithm>
#include <map>
#include <cmath>

///
/// class DTrackLinker_factory: definition for a JFactory that
/// produces space points from pseudopoints.
/// 
class DTrackLinker_factory : public JFactory<DTrackLinker> {
	public:
		
		///
		/// DTrackLinker_factory::DTrackLinker_factory():
		/// default constructor -- initializes log file
		///
		DTrackLinker_factory();
		
		///
		/// DTrackLinker_factory::~DTrackLinker_factory():
		/// default destructor -- closes log file
		///
		~DTrackLinker_factory();	

		jerror_t GetPositionAndMomentum(DFDCSegment *segment,
						DVector3 &pos, DVector3 &mom);

		const string toString(void);

	protected:
		///
		/// DTrackLinker_factory::brun():
		///
		jerror_t brun(JEventLoop *eventLoop, int eventNo);

		///
		/// DTrackLinker_factory::evnt():
		/// this is the place that links track segments into tracks
		///
		jerror_t evnt(JEventLoop *eventLoop, int eventNo);

	private:
		 JStreamLog* _log;
		 ofstream* logFile;

		 const DMagneticFieldMap *bfield;
};

#endif // DFACTORY_DTRACKLINKER_H

