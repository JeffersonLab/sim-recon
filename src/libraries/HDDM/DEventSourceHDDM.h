// Author: David Lawrence  June 24, 2004
//
//
// DEventSourceHDDM
//
/// Implements JEventSource for HDDM files
//
// Changes:	Oct 3, 2012 Yi Qiang: add classes for Cerenkov detector
//			OCt 10, 2012 Yi Qiang: modifed Cerenkov classes with general Cere hits
//			Oct 8, 2013 Yi Qiang: added dedicated object for RICH Truth Hit
//

#ifndef _JEVENT_SOURCEHDDM_H_
#define _JEVENT_SOURCEHDDM_H_

#include <vector>
#include <string>
using namespace std;

#include <pthread.h>

#include <JANA/JEventSource.h>
#include <JANA/jerror.h>
#include <JANA/JCalibration.h>

#include "hddm_s.h"
#include "hddm_s.hpp"

#include "TRACKING/DMCTrackHit.h"
#include "TRACKING/DMCThrown.h"
#include "TRACKING/DMCTrajectoryPoint.h"
#include "BCAL/DBCALSiPMHit.h"
#include "BCAL/DBCALSiPMSpectrum.h"
#include "BCAL/DBCALHit.h"
#include "BCAL/DBCALIncidentParticle.h"
#include "BCAL/DBCALTDCHit.h"
#include "BCAL/DBCALTruthShower.h"
#include "BCAL/DBCALTruthCell.h"
#include "CDC/DCDCHit.h"
#include "FDC/DFDCHit.h"
#include "FCAL/DFCALTruthShower.h"
#include "FCAL/DFCALHit.h"
#include "CCAL/DCCALTruthShower.h"
#include "CCAL/DCCALHit.h"
#include "TOF/DTOFTruth.h"
#include "TOF/DTOFHit.h"
#include "TOF/DTOFHitMC.h"
#include "START_COUNTER/DSCTruthHit.h"
#include "START_COUNTER/DSCHit.h"
#include <PID/DMCReaction.h>
#include <PID/DBeamPhoton.h>
#include <TRACKING/DTrackTimeBased.h>
#include <TAGGER/DTagger.h>
// load CERE headers, yqiang Oct 3, 2012
// modified by yqiang, Oct 10 2012
// added RichTruthHit object, yqiang, Oct 7, 2013
#include <CERE/DCereHit.h>
#include <RICH/DRichHit.h>
#include <RICH/DRichTruthHit.h>

class DEventSourceHDDM:public JEventSource
{
	public:
		DEventSourceHDDM(const char* source_name);
		virtual ~DEventSourceHDDM();		
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DEventSourceHDDM";}

		jerror_t GetEvent(JEvent &event);
		void FreeEvent(JEvent &event);
		jerror_t GetObjects(JEvent &event, JFactory_base *factory);
		
		jerror_t Extract_DMCTrackHit(s_HDDM_t *hddm_s, JFactory<DMCTrackHit> *factory);
		jerror_t GetCDCTruthHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data);
		jerror_t GetFDCTruthHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data);
		jerror_t GetBCALTruthHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data);
		jerror_t GetTOFTruthHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data);
		jerror_t GetCherenkovTruthHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data);
		jerror_t GetFCALTruthHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data);
		jerror_t GetCCALTruthHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data);
		jerror_t GetSCTruthHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data);

		jerror_t Extract_DBCALTruthShower(s_HDDM_t *hddm_s, JFactory<DBCALTruthShower> *factory);
		jerror_t Extract_DBCALTruthCell(s_HDDM_t *hddm_s, JFactory<DBCALTruthCell> *factory);
		jerror_t Extract_DBCALSiPMHit(s_HDDM_t *hddm_s, JFactory<DBCALSiPMHit> *factory);
		jerror_t Extract_DBCALSiPMSpectrum(s_HDDM_t *hddm_s, JFactory<DBCALSiPMSpectrum>* factory);
		jerror_t Extract_DBCALHit(s_HDDM_t *hddm_s, JFactory<DBCALHit> *factory);
		jerror_t Extract_DBCALIncidentParticle(s_HDDM_t *hddm_s, JFactory<DBCALIncidentParticle> *factory);
		jerror_t Extract_DBCALTDCHit(s_HDDM_t *hddm_s, JFactory<DBCALTDCHit> *factory);
		jerror_t Extract_DMCReaction(s_HDDM_t *hddm_s, JFactory<DMCReaction> *factory, JEventLoop *loop);
		jerror_t Extract_DBeamPhoton(s_HDDM_t *hddm_s, JFactory<DBeamPhoton> *factory, JEventLoop *loop);
		jerror_t Extract_DMCThrown(s_HDDM_t *hddm_s, JFactory<DMCThrown> *factory);
		jerror_t Extract_DCDCHit(s_HDDM_t *hddm_s, JFactory<DCDCHit> *factory, string tag);
		jerror_t Extract_DFDCHit(s_HDDM_t *hddm_s, JFactory<DFDCHit> *factory, string tag);
		jerror_t Extract_DFCALTruthShower(s_HDDM_t *hddm_s, JFactory<DFCALTruthShower> *factory);
		jerror_t Extract_DCCALHit(s_HDDM_t *hddm_s, JFactory<DCCALHit> *factory, JEventLoop* eventLoop, string tag);
		jerror_t Extract_DCCALTruthShower(s_HDDM_t *hddm_s, JFactory<DCCALTruthShower> *factory);
		jerror_t Extract_DFCALHit(s_HDDM_t *hddm_s, JFactory<DFCALHit> *factory, JEventLoop* eventLoop, string tag);
		jerror_t Extract_DMCTrajectoryPoint(s_HDDM_t *hddm_s, JFactory<DMCTrajectoryPoint> *factory);
		jerror_t Extract_DTOFTruth(s_HDDM_t *hddm_s,  JFactory<DTOFTruth> *factory);

		jerror_t Extract_DTOFHit( s_HDDM_t *hddm_s,  JFactory<DTOFHit>* factory, JFactory<DTOFHitMC>* factoryMC,string tag);
		jerror_t Extract_DTOFHitMC( s_HDDM_t *hddm_s,  JFactory<DTOFHitMC>* factoryMC, JFactory<DTOFHit>* factory,string tag);

		jerror_t Extract_DSCHit(s_HDDM_t *hddm_s,  JFactory<DSCHit> *factory);
		jerror_t Extract_DSCTruthHit(s_HDDM_t *hddm_s,  JFactory<DSCTruthHit> *factory);

		jerror_t Extract_DTrackTimeBased(s_HDDM_t *hddm_s,  JFactory<DTrackTimeBased> *factory, int runnumber);
		string StringToDMatrixDSym(string &str_vals, DMatrixDSym &mat, int &Nrows, int Ncols);

		jerror_t Extract_DTagger( s_HDDM_t *hddm_s,  JFactory<DTagger>* factory);

		Particle_t IDTrack(float locCharge, float locMass) const;

		// add RICH hit and Truth, yqiang Oct 3, 2012
		// modifed by yqiang, Oct 10 2012 now include both truth hits in DMCThrown
		// Oct 8, 2013, added dedicated object for RICH truth hit
		jerror_t GetRichTruthHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data);
		jerror_t Extract_DCereHit(s_HDDM_t *hddm_s, JFactory<DCereHit> *factory);
		jerror_t Extract_DRichHit(s_HDDM_t *hddm_s, JFactory<DRichHit> *factory);
		jerror_t Extract_DRichTruthHit(s_HDDM_t *hddm_s, JFactory<DRichTruthHit> *factory);

		s_iostream_t *fin;
		s_HDDM_t *hddm_s;
		bool flush_on_free;
		DApplication *dapp;
		const DMagneticFieldMap *bfield;
		const DGeometry *geom;
		
	private:
		bool initialized;
	
		pthread_mutex_t rt_mutex;
		map<s_HDDM_t*, vector<DReferenceTrajectory*> > rt_by_event;
		list<DReferenceTrajectory*> rt_pool;

		JCalibration *jcalib;
		float uscale[192],vscale[192];

};

#endif //_JEVENT_SOURCEHDDM_H_
