// $Id$
//
//    File: DFactory_DMCCheatHit.cc
// Created: Mon Apr  4 08:18:07 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include "DFactory_DMCCheatHit.h"
#include "DEventLoop.h"

//------------------------------------------------------------------
// Binary predicate used to sort cheat hits
//------------------------------------------------------------------
template<class T>
class CheatHitSort{
	public:
		bool operator()(const T &chit1,const T &chit2) const {
			if(chit1->track != chit2->track){
				return chit1->track < chit2->track;
			}
			return chit1->z < chit2->z;
		}
};

//------------------
// evnt
//------------------
derror_t DFactory_DMCCheatHit::evnt(DEventLoop *eventLoop, int eventnumber)
{
	/// This doesn't do anything. All of the work is done in  Extract_HDDM()
	/// and the GetXXXHits() methods.

	return NOERROR;
}

//------------------
// Extract_HDDM
//------------------
derror_t DFactory_DMCCheatHit::Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v)
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
	
	// Some systems will use negative phis. Force them all to
	// be in the 0 to 2pi range
	for(unsigned int i=0;i<_data.size();i++){
		DMCCheatHit *mccheathit = _data[i];
		if(mccheathit->phi<0.0)mccheathit->phi += 2.0*M_PI;
	}
	
	// sort hits by track, then z
	sort(_data.begin(), _data.end(), CheatHitSort<DMCCheatHit*>());
	
	// Copy into v
	for(unsigned int i=0; i<_data.size(); i++)v.push_back(_data[i]);

	return NOERROR;
}
//-------------------
// GetCDCHits
//-------------------
derror_t DFactory_DMCCheatHit::GetCDCHits(s_HDDM_t *hddm_s)
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
								DMCCheatHit *mccheathit = new DMCCheatHit;
								mccheathit->r			= cdcpoints->in[m].r;
								mccheathit->phi		= cdcpoints->in[m].phi;
								mccheathit->z			= cdcpoints->in[m].z;
								mccheathit->track		= cdcpoints->in[m].track;
								mccheathit->primary	= cdcpoints->in[m].primary;
								mccheathit->system	= 1;
								_data.push_back(mccheathit);
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
derror_t DFactory_DMCCheatHit::GetFDCHits(s_HDDM_t *hddm_s)
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
							DMCCheatHit *mccheathit = new DMCCheatHit;
							float x = fdcPoints->in[n].x;
							float y = fdcPoints->in[n].y;
							mccheathit->r			= sqrt(x*x + y*y);
							mccheathit->phi		= atan2(y,x);
							mccheathit->z			= fdcPoints->in[n].z;
							mccheathit->track		= fdcPoints->in[n].track;
							mccheathit->primary	= fdcPoints->in[n].primary;
							mccheathit->system	= 2;
							_data.push_back(mccheathit);
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
derror_t DFactory_DMCCheatHit::GetBCALHits(s_HDDM_t *hddm_s)
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
			DMCCheatHit *mccheathit = new DMCCheatHit;
			mccheathit->r			= barrelShowers->in[j].r;
			mccheathit->phi		= barrelShowers->in[j].phi;
			mccheathit->z			= barrelShowers->in[j].z;
			mccheathit->track		= barrelShowers->in[j].track;
			mccheathit->primary	= barrelShowers->in[j].primary;
			mccheathit->system	= 3;
			_data.push_back(mccheathit);
		}
	}

	return NOERROR;
}

//-------------------
// GetTOFHits
//-------------------
derror_t DFactory_DMCCheatHit::GetTOFHits(s_HDDM_t *hddm_s)
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
			DMCCheatHit *mccheathit = new DMCCheatHit;
			float x = tofPoints->in[j].x;
			float y = tofPoints->in[j].y;
			mccheathit->r			= sqrt(x*x + y*y);
			mccheathit->phi		= atan2(y,x);
			if(mccheathit->phi<0.0)mccheathit->phi += 2.0*M_PI;
			mccheathit->z			= tofPoints->in[j].z;
			mccheathit->track		= tofPoints->in[j].track;
			mccheathit->primary	= tofPoints->in[j].primary;
			mccheathit->system	= 4;
			_data.push_back(mccheathit);
		}
	}

	return NOERROR;
}

//-------------------
// GetCherenkovHits
//-------------------
derror_t DFactory_DMCCheatHit::GetCherenkovHits(s_HDDM_t *hddm_s)
{

	return NOERROR;
}

//-------------------
// GetFCALHits
//-------------------
derror_t DFactory_DMCCheatHit::GetFCALHits(s_HDDM_t *hddm_s)
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
			DMCCheatHit *mccheathit = new DMCCheatHit;
			float x = forwardShowers->in[j].x;
			float y = forwardShowers->in[j].y;
			mccheathit->r			= sqrt(x*x + y*y);
			mccheathit->phi		= atan2(y,x);
			if(mccheathit->phi<0.0)mccheathit->phi += 2.0*M_PI;
			mccheathit->z			= forwardShowers->in[j].z;
			mccheathit->track		= forwardShowers->in[j].track;
			mccheathit->primary	= forwardShowers->in[j].primary;
			mccheathit->system	= 6;
			_data.push_back(mccheathit);
		}
	}

	return NOERROR;
}

//-------------------
// GetUPVHits
//-------------------
derror_t DFactory_DMCCheatHit::GetUPVHits(s_HDDM_t *hddm_s)
{

	return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DMCCheatHit::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	GetNrows();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("row:   r(cm): phi(rad):  z(cm): track: primary:    system:");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DMCCheatHit *mccheathit = _data[i];

		printnewrow();
		
		printcol("%d", i);
		printcol("%3.1f", mccheathit->r);
		printcol("%1.3f", mccheathit->phi);
		printcol("%3.1f", mccheathit->z);
		printcol("%d", mccheathit->track);
		printcol(mccheathit->primary ? "Y":"N");

		char *system = "<unknown>";
		switch(mccheathit->system){
			case 1: system = "CDC";				break;
			case 2: system = "FDC";				break;
			case 3: system = "BCAL";			break;
			case 4: system = "TOF";				break;
			case 5: system = "Cherenkov";		break;
			case 6: system = "FCAL";			break;
			case 7: system = "UPV";				break;
		}
		printcol(system);
		printrow();
	}

	return _table;
}
