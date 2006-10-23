// $Id$
//
// Author: David Lawrence  June 24, 2004
//
//
// DEventSourceHDDM methods
//

#include <iostream>
#include <iomanip>
using namespace std;

#include <JANA/JFactory_base.h>
#include <JANA/JEventLoop.h>
#include <JANA/JEvent.h>

#include "DEventSourceHDDM.h"
#include "FDC/DFDCGeometry.h"

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


//----------------
// Constructor
//----------------
DEventSourceHDDM::DEventSourceHDDM(const char* source_name):JEventSource(source_name)
{
	/// Constructor for DEventSourceHDDM object
	fin = open_s_HDDM((char*)source_name);
	hddm_s = NULL;
	
	if(fin)source_is_open = 1;
	flush_on_free = true;
}

//----------------
// Destructor
//----------------
DEventSourceHDDM::~DEventSourceHDDM()
{
	// Close file and set pointers to NULL
	if(fin)close_s_HDDM(fin);
	hddm_s = NULL;
	fin = NULL;
}

//----------------
// GetEvent
//----------------
jerror_t DEventSourceHDDM::GetEvent(JEvent &event)
{
	/// Implementation of JEventSource virtual function
		
	if(!fin){
		return EVENT_SOURCE_NOT_OPEN;
	}
	hddm_s = read_s_HDDM(fin);
	++Nevents_read;
	
	// each open HDDM file takes up about 1M of memeory so it's
	// worthwhile to close it as soon as we can
	if(!hddm_s){
		if(fin)close_s_HDDM(fin);
		fin = NULL;
	}

	int event_number = -1;
	int run_number = -1;

	// Copy the reference info into the JEvent object
	event.SetJEventSource(this);
	event.SetEventNumber(event_number);
	event.SetRunNumber(run_number);
	event.SetRef(hddm_s);

	return hddm_s==NULL ? NO_MORE_EVENTS_IN_SOURCE:NOERROR;
}

//----------------
// FreeEvent
//----------------
void DEventSourceHDDM::FreeEvent(JEvent &event)
{
	s_HDDM_t *my_hddm_s = (s_HDDM_t*)event.GetRef();
	if(flush_on_free && my_hddm_s!=NULL)flush_s_HDDM(my_hddm_s, 0);
}

//----------------
// GetObjects
//----------------
jerror_t DEventSourceHDDM::GetObjects(JEvent &event, JFactory_base *factory)
{
	/// This gets called through the virtual method of the
	/// JEventSource base class. It creates the objects of the type
	/// on which factory is based. It uses the s_HDDM_t object
	/// kept in the ref field of the JEvent object passed.

	// We must have a factory to hold the data
	if(!factory)throw RESOURCE_UNAVAILABLE;
	
	// HDDM doesn't support tagged factories
	if(strcmp(factory->Tag(), ""))return OBJECT_NOT_AVAILABLE;
	
	// The ref field of the JEvent is just the s_HDDM_t pointer.
	s_HDDM_t *my_hddm_s = (s_HDDM_t*)event.GetRef();
	if(!my_hddm_s)throw RESOURCE_UNAVAILABLE;
	
	// Get name of data class we're trying to extract
	string dataClassName = factory->dataClassName();
	
	if(dataClassName =="DMCTrackHit")
		return Extract_DMCTrackHit(my_hddm_s, dynamic_cast<JFactory<DMCTrackHit>*>(factory));

	if(dataClassName =="DMCThrown")
		return Extract_DMCThrown(my_hddm_s, dynamic_cast<JFactory<DMCThrown>*>(factory));

	if(dataClassName =="DBCALHit")
		return Extract_DBCALHit(my_hddm_s, dynamic_cast<JFactory<DBCALHit>*>(factory));

	if(dataClassName =="DCDCHit")
		return Extract_DCDCHit(my_hddm_s, dynamic_cast<JFactory<DCDCHit>*>(factory));

	if(dataClassName =="DFDCHit")
		return Extract_DFDCHit(my_hddm_s, dynamic_cast<JFactory<DFDCHit>*>(factory));

	
	return OBJECT_NOT_AVAILABLE;
}


//------------------
// Extract_DMCTrackHit
//------------------
jerror_t DEventSourceHDDM::Extract_DMCTrackHit(s_HDDM_t *hddm_s, JFactory<DMCTrackHit> *factory)
{
	/// Copies the data from the given hddm_s structure. This is called
	/// from JEventSourceHDDM::GetObjects. If factory is NULL, this
	/// returns OBJECT_NOT_AVAILABLE immediately.
	
	if(factory==NULL)return OBJECT_NOT_AVAILABLE;
	
	// The following routines will create DMCTrackHit objects and add them
	// to data.
	vector<DMCTrackHit*> data;
	GetCDCTruthHits(hddm_s, data);
	GetFDCTruthHits(hddm_s, data);
	GetBCALTruthHits(hddm_s, data);
	GetTOFTruthHits(hddm_s, data);
	GetCherenkovTruthHits(hddm_s, data);
	GetFCALTruthHits(hddm_s, data);
	GetUPVTruthHits(hddm_s, data);
	
	// sort hits by z
	sort(data.begin(), data.end(), MCTrackHitSort_C);
	
	// Some systems will use negative phis. Force them all to
	// be in the 0 to 2pi range
	for(unsigned int i=0;i<data.size();i++){
		DMCTrackHit *mctrackhit = data[i];
		if(mctrackhit->phi<0.0)mctrackhit->phi += 2.0*M_PI;
	}
	
	// Copy into factory
	factory->CopyTo(data);

	return NOERROR;
}
//-------------------
// GetCDCTruthHits
//-------------------
jerror_t DEventSourceHDDM::GetCDCTruthHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data)
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
			data.push_back(mctrackhit);
		}
	}
		
	return NOERROR;
}

//-------------------
// GetFDCTruthHits
//-------------------
jerror_t DEventSourceHDDM::GetFDCTruthHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data)
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
				data.push_back(mctrackhit);
			}
		}
	}

	return NOERROR;
}

//-------------------
// GetBCALTruthHits
//-------------------
jerror_t DEventSourceHDDM::GetBCALTruthHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data)
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
			data.push_back(mctrackhit);
		}
	}

	return NOERROR;
}

//-------------------
// GetTOFTruthHits
//-------------------
jerror_t DEventSourceHDDM::GetTOFTruthHits(s_HDDM_t *hddm_s,  vector<DMCTrackHit*>& data)
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
			data.push_back(mctrackhit);
		}
	}
		
	return NOERROR;
}

//-------------------
// GetCherenkovTruthHits
//-------------------
jerror_t DEventSourceHDDM::GetCherenkovTruthHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data)
{

	return NOERROR;
}

//-------------------
// GetFCALTruthHits
//-------------------
jerror_t DEventSourceHDDM::GetFCALTruthHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data)
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
		for(unsigned int j=0; j<fcalTruthShowers->mult; j++, fcalTruthShower++){
			float x = fcalTruthShower->x;
			float y = fcalTruthShower->y;
			DMCTrackHit *mctrackhit = new DMCTrackHit;
			mctrackhit->r			= sqrt(x*x + y*y);
			mctrackhit->phi		= atan2(y,x);
			mctrackhit->z			= fcalTruthShower->z;
			mctrackhit->track		= fcalTruthShower->track;
			mctrackhit->primary	= fcalTruthShower->primary;
			mctrackhit->system	= SYS_FCAL;
			data.push_back(mctrackhit);
		}
	}

	return NOERROR;
}

//-------------------
// GetUPVTruthHits
//-------------------
jerror_t DEventSourceHDDM::GetUPVTruthHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data)
{

	return NOERROR;
}

//------------------
// Extract_DMCThrown
//------------------
jerror_t DEventSourceHDDM::Extract_DBCALHit(s_HDDM_t *hddm_s, JFactory<DBCALHit> *factory)
{
	/// Copies the data from the given hddm_s structure. This is called
	/// from JEventSourceHDDM::GetObjects. If factory is NULL, this
	/// returns OBJECT_NOT_AVAILABLE immediately.
	
	if(factory==NULL)return OBJECT_NOT_AVAILABLE;
	
	vector<DBCALHit*> data;

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL ||
			hits->barrelEMcal == HDDM_NULL ||
			hits->barrelEMcal->bcalCells == HDDM_NULL)continue;
		
		// Loop over BCAL cells
		s_BcalCells_t *cells = hits->barrelEMcal->bcalCells;
		for(unsigned int j=0;j<cells->mult;j++){
			s_BcalCell_t *cell = &cells->in[j];
			if(cell->bcalUpstreamHits != HDDM_NULL){
				for(unsigned int k=0; k<cell->bcalUpstreamHits->mult; k++){
					s_BcalUpstreamHit_t *upstreamhit = &cell->bcalUpstreamHits->in[k];
					
					DBCALHit *bcalhit = new DBCALHit();
					bcalhit->module = cell->module;
					bcalhit->layer = cell->layer;
					bcalhit->sector = cell->sector;
					bcalhit->end = DBCALHit::UPSTREAM;
					bcalhit->E = upstreamhit->E;
					bcalhit->t = upstreamhit->t;
					data.push_back(bcalhit);
				}
			}

			if(cell->bcalDownstreamHits != HDDM_NULL){
				for(unsigned int k=0; k<cell->bcalDownstreamHits->mult; k++){
					s_BcalDownstreamHit_t *downstreamhit = &cell->bcalDownstreamHits->in[k];
					
					DBCALHit *bcalhit = new DBCALHit();
					bcalhit->module = cell->module;
					bcalhit->layer = cell->layer;
					bcalhit->sector = cell->sector;
					bcalhit->end = DBCALHit::DOWNSTREAM;
					bcalhit->E = downstreamhit->E;
					bcalhit->t = downstreamhit->t;
					data.push_back(bcalhit);
				}
			}
		} // j   (cells)
	} // i   (physicsEvents)
	
	// Copy into factory
	factory->CopyTo(data);

	return NOERROR;
}

//------------------
// Extract_DMCThrown
//------------------
jerror_t DEventSourceHDDM::Extract_DMCThrown(s_HDDM_t *hddm_s,  JFactory<DMCThrown> *factory)
{
	/// Copies the data from the given hddm_s structure. This is called
	/// from JEventSourceHDDM::GetObjects. If factory is NULL, this
	/// returns OBJECT_NOT_AVAILABLE immediately.
	
	if(factory==NULL)return OBJECT_NOT_AVAILABLE;
	
	vector<DMCThrown*> data;

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	for(unsigned int i=0; i<PE->mult; i++){
		// ------------ Reactions --------------
		s_Reactions_t *reactions=PE->in[i].reactions;
		if(!reactions)continue;

		for(unsigned int j=0; j<reactions->mult; j++){
			s_Vertices_t *vertices = reactions->in[j].vertices;
			if(vertices){
				for(unsigned int k=0; k<vertices->mult; k++){
					s_Origin_t *origin = vertices->in[k].origin;
					s_Products_t *products = vertices->in[k].products;
					if(products && origin){
						for(unsigned int m=0;m<products->mult;m++){
							s_Product_t *product = &products->in[m];
							
							DMCThrown *mcthrown = new DMCThrown;
							mcthrown->x = origin->vx;
							mcthrown->y = origin->vy;
							mcthrown->z = origin->vz;
							mcthrown->type = product->type;
							mcthrown->q = (float)ParticleCharge(product->type);
							mcthrown->mass = ParticleMass(product->type);
							mcthrown->E = product->momentum->E;
							double px = product->momentum->px;
							double py = product->momentum->py;
							double pz = product->momentum->pz;
							mcthrown->p = sqrt(px*px + py*py + pz*pz);
							mcthrown->phi = atan2(py, px);
							if(mcthrown->phi<0.0)mcthrown->phi += 2.0*M_PI;
							mcthrown->theta = acos(pz/mcthrown->p);
							data.push_back(mcthrown);
						}
					}
				}
			}
		}
	}
	
	// Copy into factory
	factory->CopyTo(data);

	return NOERROR;
}

//------------------
// Extract_DCDCHit
//------------------
jerror_t DEventSourceHDDM::Extract_DCDCHit(s_HDDM_t *hddm_s,  JFactory<DCDCHit> *factory)
{
	/// Copies the data from the given hddm_s structure. This is called
	/// from JEventSourceHDDM::GetObjects. If factory is NULL, this
	/// returns OBJECT_NOT_AVAILABLE immediately.
	
	if(factory==NULL)return OBJECT_NOT_AVAILABLE;
	
	vector<DCDCHit*> data;

	// Acquire the pointer to the physics events
	s_PhysicsEvents_t* allEvents = hddm_s->physicsEvents;
	if(!allEvents)return NOERROR;
       
	for (unsigned int m=0; m < allEvents->mult; m++) {
		// Acquire the pointer to the overall hits section of the data
		s_HitView_t *hits = allEvents->in[m].hitView;
		
		if (hits == HDDM_NULL)return NOERROR;
		if (hits->centralDC == HDDM_NULL)return NOERROR;
		if (hits->centralDC->cdcStraws == HDDM_NULL)return NOERROR;
		for(unsigned int k=0; k<hits->centralDC->cdcStraws->mult; k++){
			s_CdcStraw_t *cdcstraw = &hits->centralDC->cdcStraws->in[k];
			for(unsigned int j=0; j<cdcstraw->cdcStrawHits->mult; j++){
				s_CdcStrawHit_t *strawhit = &cdcstraw->cdcStrawHits->in[j];
				
				DCDCHit *hit = new DCDCHit;
				hit->ring = cdcstraw->ring;
				hit->straw = cdcstraw->straw;
				hit->dE = strawhit->dE;
				hit->t = strawhit->t;

				data.push_back(hit);
			}
		}
	}
	
	// Copy into factory
	factory->CopyTo(data);

	return NOERROR;
}


//------------------
// Extract_DFDCHit
//------------------
jerror_t DEventSourceHDDM::Extract_DFDCHit(s_HDDM_t *hddm_s,  JFactory<DFDCHit> *factory)
{
	/// Copies the data from the given hddm_s structure. This is called
	/// from JEventSourceHDDM::GetObjects. If factory is NULL, this
	/// returns OBJECT_NOT_AVAILABLE immediately.
	
	if(factory==NULL)return OBJECT_NOT_AVAILABLE;
	
	vector<DFDCHit*> data;

	// Acquire the pointer to the physics events
	s_PhysicsEvents_t* allEvents = hddm_s->physicsEvents;
	if(!allEvents) {
	  //throw JException("Attempt to get physics events from HDDM source failed.");
		return NOERROR;
	}
       
	for (unsigned int m=0; m < allEvents->mult; m++) {
	
		// Acquire the pointer to the overall hits section of the data
		s_HitView_t *hits = allEvents->in[m].hitView;
		
		if (hits == HDDM_NULL) {
		  //throw JException("HDDM source has no hits.");
			return NOERROR;
		}

		if (hits->forwardDC == HDDM_NULL) {
		  //throw JException("HDDM source has no forwardDC information.");
			return NOERROR;
		}

		if (hits->forwardDC->fdcChambers == HDDM_NULL) {
		  // throw JException("HDDM source has no hits in the FDC.");		
			return NOERROR;
		}

		// Acquire the pointer to the beginning of the FDC hit tree
		s_FdcChambers_t* fdcChamberSet = hits->forwardDC->fdcChambers;

		for (unsigned int i=0; i < fdcChamberSet->mult; i++) {
			// Each chamber in the ChamberSet has a wire set and a strip set
			s_FdcChamber_t &fdcChamber 		= fdcChamberSet->in[i];		
			s_FdcAnodeWires_t* wireSet 		= fdcChamber.fdcAnodeWires;
			s_FdcCathodeStrips_t* stripSet 	= fdcChamber.fdcCathodeStrips;
		
			// Each set of wires has (obviously) wires inside of it, and each wire
			// may have one or more hits on it. Make a DFDCHit object for each one
			// of these hits.
			for (unsigned int j=0; j < wireSet->mult; j++) {
				s_FdcAnodeWire_t anodeWire		= wireSet->in[j];
				s_FdcAnodeHits_t* wireHitSet	= anodeWire.fdcAnodeHits;
				for (unsigned int k=0; k < wireHitSet->mult; k++) {
					s_FdcAnodeHit_t wireHit		= wireHitSet->in[k];
					DFDCHit* newHit				= new DFDCHit();
					newHit->layer		 		= fdcChamber.layer;
					newHit->module		 		= fdcChamber.module;
					newHit->element				= anodeWire.wire;
					newHit->dE					= wireHit.dE;
					newHit->t					= wireHit.t;
					newHit->plane				= 2;
					newHit->type				= 0;
					newHit->gPlane				= DFDCGeometry::gPlane(newHit);
					newHit->gLayer				= DFDCGeometry::gLayer(newHit);
					newHit->r					= DFDCGeometry::getWireR(newHit); 
					data.push_back(newHit);
				}
			}
		
			// Ditto for the cathodes.
			for (unsigned int j=0; j < stripSet->mult; j++) {
				s_FdcCathodeStrip_t cathodeStrip = stripSet->in[j];
				s_FdcCathodeHits_t* stripHitSet = cathodeStrip.fdcCathodeHits;
				for (unsigned int k=0; k < stripHitSet->mult; k++) {
					s_FdcCathodeHit_t stripHit	= stripHitSet->in[k];
					DFDCHit* newHit				= new DFDCHit();
					newHit->layer				= fdcChamber.layer;
					newHit->module				= fdcChamber.module;
					newHit->element				= cathodeStrip.strip;
					newHit->plane				= cathodeStrip.plane;
					newHit->dE					= stripHit.dE;
					newHit->t					= stripHit.t;
					newHit->type				= 1;
					newHit->gPlane				= DFDCGeometry::gPlane(newHit);	 
					newHit->gLayer				= DFDCGeometry::gLayer(newHit);
					newHit->r					= DFDCGeometry::getStripR(newHit);
					data.push_back(newHit);
				}
			}	
		}
	}
	
	// Copy into factory
	factory->CopyTo(data);

	return NOERROR;
}

