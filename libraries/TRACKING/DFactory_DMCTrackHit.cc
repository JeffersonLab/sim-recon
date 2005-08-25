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
	sort(_data.begin(), _data.end(), MCTrackHitSort());
	
	// Set id values of all hits to be unique. At the same time ...
	// Some systems will use negative phis. Force them all to
	// be in the 0 to 2pi range
	identifier_t idcntr = 1;
	for(unsigned int i=0;i<_data.size();i++){
		DMCTrackHit *mctrackhit = _data[i];
		if(mctrackhit->phi<0.0)mctrackhit->phi += 2.0*M_PI;
		mctrackhit->id = idcntr++;
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
		// ------------ CdcPoints, Hits --------------
		s_Rings_t *rings=NULL;
		if(PE->in[i].hitView)
			if(PE->in[i].hitView->centralDC)
				rings = PE->in[i].hitView->centralDC->rings;
		if(rings){
			for(unsigned int j=0;j<rings->mult;j++){
				//float radius = rings->in[j].radius;
				s_Straws_t *straws = rings->in[j].straws;
				if(straws){
					for(unsigned int k=0;k<straws->mult;k++){
						//float phim = straws->in[k].phim;
						s_CdcPoints_t *cdcpoints = straws->in[k].cdcPoints;
						if(cdcpoints){
							for(unsigned int m=0;m<cdcpoints->mult;m++){
								DMCTrackHit *mctrackhit = new DMCTrackHit;
								mctrackhit->r			= cdcpoints->in[m].r;
								mctrackhit->phi		= cdcpoints->in[m].phi;
								mctrackhit->z			= cdcpoints->in[m].z;
								mctrackhit->track		= cdcpoints->in[m].track;
								mctrackhit->primary	= cdcpoints->in[m].primary;
								mctrackhit->system	= SYS_CDC;
								_data.push_back(mctrackhit);
							}
						}
					}
				}
			}
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
		s_Chambers_t *chambers = NULL;
		if(PE->in[i].hitView)
			if(PE->in[i].hitView->forwardDC)
				chambers = PE->in[i].hitView->forwardDC->chambers;
		if(!chambers)continue;
		
		for(unsigned int j=0;j<chambers->mult;j++){

			s_AnodePlanes_t *anodeplanes = chambers->in[j].anodePlanes;
			if(anodeplanes){
			
				for(unsigned int k=0;k<anodeplanes->mult;k++){
					s_Wires_t *wires = anodeplanes->in[k].wires;
					if(!wires)continue;
				
					for(unsigned int m=0;m<wires->mult;m++){
						s_FdcPoints_t *fdcPoints = wires->in[m].fdcPoints;
						if(!fdcPoints)continue;
						for(unsigned int n=0;n<fdcPoints->mult;n++){
							DMCTrackHit *mctrackhit = new DMCTrackHit;
							float x = fdcPoints->in[n].x;
							float y = fdcPoints->in[n].y;
							mctrackhit->r			= sqrt(x*x + y*y);
							mctrackhit->phi		= atan2(y,x);
							mctrackhit->z			= fdcPoints->in[n].z;
							mctrackhit->track		= fdcPoints->in[n].track;
							mctrackhit->primary	= fdcPoints->in[n].primary;
							mctrackhit->system	= SYS_FDC;
							_data.push_back(mctrackhit);
						}
					}
				}
			}

			s_CathodePlanes_t *cathodeplanes = chambers->in[j].cathodePlanes;
			if(cathodeplanes){
			
				for(unsigned int k=0;k<cathodeplanes->mult;k++){
					//float tau = cathodeplanes->in[k].tau;
					//float z = cathodeplanes->in[k].z;
					s_Strips_t *strips = cathodeplanes->in[k].strips;
					if(!strips)continue;
				
					for(unsigned int m=0;m<strips->mult;m++){
						// Just skip cathode hits
					}
				}
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
		s_BarrelShowers_t *barrelShowers = NULL;
		if(PE->in[i].hitView)
			if(PE->in[i].hitView->barrelEMcal)
				barrelShowers = PE->in[i].hitView->barrelEMcal->barrelShowers;
		if(!barrelShowers)continue;
		
		for(unsigned int j=0;j<barrelShowers->mult;j++){
			DMCTrackHit *mctrackhit = new DMCTrackHit;
			mctrackhit->r			= barrelShowers->in[j].r;
			mctrackhit->phi		= barrelShowers->in[j].phi;
			mctrackhit->z			= barrelShowers->in[j].z;
			mctrackhit->track		= barrelShowers->in[j].track;
			mctrackhit->primary	= barrelShowers->in[j].primary;
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
		s_TofPoints_t *tofPoints = NULL;
		if(PE->in[i].hitView)
			if(PE->in[i].hitView->forwardTOF)
				tofPoints = PE->in[i].hitView->forwardTOF->tofPoints;
		if(!tofPoints)continue;
		
		for(unsigned int j=0;j<tofPoints->mult;j++){
			DMCTrackHit *mctrackhit = new DMCTrackHit;
			float x = tofPoints->in[j].x;
			float y = tofPoints->in[j].y;
			mctrackhit->r			= sqrt(x*x + y*y);
			mctrackhit->phi		= atan2(y,x);
			if(mctrackhit->phi<0.0)mctrackhit->phi += 2.0*M_PI;
			mctrackhit->z			= tofPoints->in[j].z;
			mctrackhit->track		= tofPoints->in[j].track;
			mctrackhit->primary	= tofPoints->in[j].primary;
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
		s_ForwardShowers_t *forwardShowers = NULL;
		if(PE->in[i].hitView)
			if(PE->in[i].hitView->forwardEMcal)
				forwardShowers = PE->in[i].hitView->forwardEMcal->forwardShowers;
		if(!forwardShowers)continue;
		
		for(unsigned int j=0;j<forwardShowers->mult;j++){
			DMCTrackHit *mctrackhit = new DMCTrackHit;
			float x = forwardShowers->in[j].x;
			float y = forwardShowers->in[j].y;
			mctrackhit->r			= sqrt(x*x + y*y);
			mctrackhit->phi		= atan2(y,x);
			if(mctrackhit->phi<0.0)mctrackhit->phi += 2.0*M_PI;
			mctrackhit->z			= forwardShowers->in[j].z;
			mctrackhit->track		= forwardShowers->in[j].track;
			mctrackhit->primary	= forwardShowers->in[j].primary;
			mctrackhit->system	= SYS_FCAL;
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
		
		printcol("%d", mctrackhit->id);
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
