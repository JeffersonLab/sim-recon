// $Id$
//
//    File: DTrackCandidate_factory_CDCCOSMIC.h
// Created: Sat Jun 28 16:50:07 EDT 2014
// Creator: davidl (on Darwin harriet.local 13.2.0 i386)
//

#ifndef _DTrackCandidate_factory_CDCCOSMIC_
#define _DTrackCandidate_factory_CDCCOSMIC_

#include <JANA/JFactory.h>
#include <TRACKING/DTrackCandidate.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <HDGEOMETRY/DMagneticFieldMapNoField.h>

class DTrackCandidate_factory_CDCCOSMIC:public jana::JFactory<DTrackCandidate>{
	public:
		DTrackCandidate_factory_CDCCOSMIC():rt(NULL),bfield(NULL){};
		~DTrackCandidate_factory_CDCCOSMIC(){};
		const char* Tag(void){return "CDCCOSMIC";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		DReferenceTrajectory *rt;
		DMagneticFieldMapNoField *bfield;
		
		void CalcChisq(DTrackCandidate *can, vector<const DCDCTrackHit*> &axial_hits, vector<const DCDCTrackHit*> &stereo_hits);
};

#endif // _DTrackCandidate_factory_CDCCOSMIC_

