// $Id$
//
//    File: DFactory_DMCTrackHit.cc
// Created: Mon Apr  4 08:18:07 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include "DFactory_DMCTrackHit.h"
#include "DEventLoop.h"
#include "DTrackHit.h"
#include "GlueX.h"

//------------------------------------------------------------------
// Binary predicate used to sort hits
//------------------------------------------------------------------
class MCTrackHitSort{
	public:
		bool operator()(DMCTrackHit* const &thit1, DMCTrackHit* const &thit2) const {
			return thit1->z < thit2->z;
		}
};

bool MCTrackHitSort_C(DMCTrackHit* const &thit1, DMCTrackHit* const &thit2) {
	return thit1->z < thit2->z;
}


//------------------
// evnt
//------------------
derror_t DFactory_DMCTrackHit::evnt(DEventLoop *eventLoop, int eventnumber)
{
	/// This doesn't do anything. All of the work is done in  Extract_HDDM()
	/// and the GetXXXHits() methods.

	return NOERROR;
}

//------------------
// Extract_HDDM
//------------------
derror_t DFactory_DMCTrackHit::Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v)
{
	/// Copies the data from the given hddm_s structure. This is called
	/// from DEventSourceHDDM::GetObjects.
	
	v.clear();

	// These will put the data in the correctly typed _data vector. The
	// pointers are copied to the v vector below. Note this means that 
	// _data will be overwritten by the contents of v later by DEvent::Get()
	// but the contents will be identical. This somewhat convoluted way
	// of doing things is needed to implement a generic API for event sources.
	GetCDCHits(hddm_s);
	GetFDCHits(hddm_s);
	GetBCALHits(hddm_s);
	GetTOFHits(hddm_s);
	GetCherenkovHits(hddm_s);
	GetFCALHits(hddm_s);
	GetUPVHits(hddm_s);
	
	// sort hits by z
	sort(_data.begin(), _data.end(), MCTrackHitSort_C);
	
	// Some systems will use negative phis. Force them all to
	// be in the 0 to 2pi range
	for(unsigned int i=0;i<_data.size();i++){
		DMCTrackHit *mctrackhit = _data[i];
		if(mctrackhit->phi<0.0)mctrackhit->phi += 2.0*M_PI;
	}
	
	// Copy into v
	for(unsigned int i=0; i<_data.size(); i++)v.push_back(_data[i]);

	return NOERROR;
}
//-------------------
// GetCDCHits
//-------------------
derror_t DFactory_DMCTrackHit::GetCDCHits(s_HDDM_t *hddm_s)
{
	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL ||
			hits->centralDC == HDDM_NULL ||
			hits->centralDC->cdcTruthPoints == HDDM_NULL)continue;
		
		s_CdcTruthPoints_t *cdctruthpoints = hits->centralDC->cdcTruthPoints;
		s_CdcTruthPoint_t *cdctruthpoint = cdctruthpoints->in;
		for(unsigned int j=0; j<cdctruthpoints->mult; j++, cdctruthpoint++){
			DMCTrackHit *mctrackhit = new DMCTrackHit;
			mctrackhit->r			= cdctruthpoint->r;
			mctrackhit->phi		= cdctruthpoint->phi;
			mctrackhit->z			= cdctruthpoint->z;
			mctrackhit->track		= cdctruthpoint->track;
			mctrackhit->primary	= cdctruthpoint->primary;
			mctrackhit->system	= SYS_CDC;
			_data.push_back(mctrackhit);
		}
	}
		
	return NOERROR;
}

//-------------------
// GetFDCHits
//-------------------
derror_t DFactory_DMCTrackHit::GetFDCHits(s_HDDM_t *hddm_s)
{
	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL ||
			hits->forwardDC == HDDM_NULL ||
			hits->forwardDC->fdcChambers == HDDM_NULL)continue;

		s_FdcChambers_t* fdcChambers = hits->forwardDC->fdcChambers;
		s_FdcChamber_t *fdcChamber = fdcChambers->in;
		for(unsigned int j=0; j<fdcChambers->mult; j++, fdcChamber++){
			s_FdcTruthPoints_t *fdcTruthPoints = fdcChamber->fdcTruthPoints;
			if(fdcTruthPoints == HDDM_NULL)continue;
			
			s_FdcTruthPoint_t *truth = fdcTruthPoints->in;
			for(unsigned int k=0; k<fdcTruthPoints->mult; k++, truth++){
				float x = truth->x;
				float y = truth->y;
				DMCTrackHit *mctrackhit = new DMCTrackHit;
				mctrackhit->r			= sqrt(x*x + y*y);
				mctrackhit->phi		= atan2(y,x);
				mctrackhit->z			= truth->z;
				mctrackhit->track		= truth->track;
				mctrackhit->primary	= truth->primary;
				mctrackhit->system	= SYS_FDC;
				_data.push_back(mctrackhit);
			}
		}
	}

	return NOERROR;
}

//-------------------
// GetBCALHits
//-------------------
derror_t DFactory_DMCTrackHit::GetBCALHits(s_HDDM_t *hddm_s)
{
	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL ||
			hits->barrelEMcal == HDDM_NULL ||
			hits->barrelEMcal->bcalTruthShowers == HDDM_NULL)continue;
		
		s_BcalTruthShowers_t *bcalTruthShowers = hits->barrelEMcal->bcalTruthShowers;
		s_BcalTruthShower_t *bcalTruthShower = bcalTruthShowers->in;
		for(unsigned int j=0; j<bcalTruthShowers->mult; j++, bcalTruthShower++){
			DMCTrackHit *mctrackhit = new DMCTrackHit;
			mctrackhit->r			= bcalTruthShower->r;
			mctrackhit->phi		= bcalTruthShower->phi;
			mctrackhit->z			= bcalTruthShower->z;
			mctrackhit->track		= bcalTruthShower->track;
			mctrackhit->primary	= bcalTruthShower->primary;
			mctrackhit->system	= SYS_BCAL;
			_data.push_back(mctrackhit);
		}
	}

	return NOERROR;
}

//-------------------
// GetTOFHits
//-------------------
derror_t DFactory_DMCTrackHit::GetTOFHits(s_HDDM_t *hddm_s)
{
	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL ||
			hits->forwardTOF == HDDM_NULL ||
			hits->forwardTOF->ftofTruthPoints == HDDM_NULL)continue;
		
		s_FtofTruthPoints_t *ftoftruthpoints = hits->forwardTOF->ftofTruthPoints;
		s_FtofTruthPoint_t *ftoftruthpoint = ftoftruthpoints->in;
		for(unsigned int j=0; j<ftoftruthpoints->mult; j++, ftoftruthpoint++){
			float x = ftoftruthpoint->x;
			float y = ftoftruthpoint->y;
			DMCTrackHit *mctrackhit = new DMCTrackHit;
			mctrackhit->r			= sqrt(x*x + y*y);
			mctrackhit->phi		= atan2(y,x);
			mctrackhit->z			= ftoftruthpoint->z;
			mctrackhit->track		= ftoftruthpoint->track;
			mctrackhit->primary	= ftoftruthpoint->primary;
			mctrackhit->system	= SYS_TOF;
			_data.push_back(mctrackhit);
		}
	}
		
	return NOERROR;
}

//-------------------
// GetCherenkovHits
//-------------------
derror_t DFactory_DMCTrackHit::GetCherenkovHits(s_HDDM_t *hddm_s)
{

	return NOERROR;
}

//-------------------
// GetFCALHits
//-------------------
derror_t DFactory_DMCTrackHit::GetFCALHits(s_HDDM_t *hddm_s)
{
	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL ||
			hits->forwardEMcal == HDDM_NULL ||
			hits->forwardEMcal->fcalTruthShowers == HDDM_NULL)continue;
		
		s_FcalTruthShowers_t *fcalTruthShowers = hits->forwardEMcal->fcalTruthShowers;
		s_FcalTruthShower_t *fcalTruthShower = fcalTruthShowers->in;
		for(unsigned int j=0; j<fcalTruthShowers->mult; j++, fcalTruthShowers++){
			float x = fcalTruthShower->x;
			float y = fcalTruthShower->y;
			DMCTrackHit *mctrackhit = new DMCTrackHit;
			mctrackhit->r			= sqrt(x*x + y*y);
			mctrackhit->phi		= atan2(y,x);
			mctrackhit->z			= fcalTruthShower->z;
			mctrackhit->track		= fcalTruthShower->track;
			mctrackhit->primary	= fcalTruthShower->primary;
			mctrackhit->system	= SYS_TOF;
			_data.push_back(mctrackhit);
		}
	}

	return NOERROR;
}

//-------------------
// GetUPVHits
//-------------------
derror_t DFactory_DMCTrackHit::GetUPVHits(s_HDDM_t *hddm_s)
{

	return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DMCTrackHit::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	GetNrows();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("id:   r(cm): phi(rad):  z(cm): track: primary:    system:");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DMCTrackHit *mctrackhit = _data[i];

		printnewrow();
		
		printcol("%x", mctrackhit->id);
		printcol("%3.1f", mctrackhit->r);
		printcol("%1.3f", mctrackhit->phi);
		printcol("%3.1f", mctrackhit->z);
		printcol("%d", mctrackhit->track);
		printcol(mctrackhit->primary ? "Y":"N");
		printcol(SystemName(mctrackhit->system));
		printrow();
	}

	return _table;
}
