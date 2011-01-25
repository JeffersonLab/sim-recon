// $Id$
//
// Author: David Lawrence  June 24, 2004
//
// changes: Wed Jun 20 17:08:13 EDT 2007 B. Zihlmann
//          modify TOF section to add several new variables incuding the 
//          GEANT particle type to the Truth hits and the hit and track-hit list.
//
// DEventSourceHDDM methods
//

#include <iostream>
#include <iomanip>
using namespace std;

#include <JANA/JFactory_base.h>
#include <JANA/JEventLoop.h>
#include <JANA/JEvent.h>

#include "BCAL/DBCALGeometry.h"

#include <DVector2.h>
#include <DEventSourceHDDM.h>
#include <FDC/DFDCGeometry.h>
#include <FCAL/DFCALGeometry.h>
#include <FCAL/DFCALHit.h>
#include <CCAL/DCCALGeometry.h>
#include <CCAL/DCCALHit.h>

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
	
	pthread_mutex_init(&rt_mutex, NULL);
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

	// each open HDDM file takes up about 1M of memory so it's
	// worthwhile to close it as soon as we can
	if(!hddm_s){
		if(fin)close_s_HDDM(fin);
		fin = NULL;
		return NO_MORE_EVENTS_IN_SOURCE;
	}
	
	++Nevents_read;
	
	int event_number = -1;
	int run_number = -1;
	
	// Get event/run numbers from HDDM
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(PE && PE->mult>0){
		event_number = PE->in[0].eventNo;
		run_number = PE->in[0].runNo;
	}

	// Copy the reference info into the JEvent object
	event.SetJEventSource(this);
	event.SetEventNumber(event_number);
	event.SetRunNumber(run_number);
	event.SetRef(hddm_s);
	
	return NOERROR;
}

//----------------
// FreeEvent
//----------------
void DEventSourceHDDM::FreeEvent(JEvent &event)
{
	s_HDDM_t *my_hddm_s = (s_HDDM_t*)event.GetRef();
	if(flush_on_free && my_hddm_s!=NULL)flush_s_HDDM(my_hddm_s, 0);

	// Check for DReferenceTrajectory objects we need to delete
	pthread_mutex_lock(&rt_mutex);
	map<s_HDDM_t*, vector<DReferenceTrajectory*> >::iterator iter = rt_by_event.find(my_hddm_s);
	if(iter != rt_by_event.end()){
		vector<DReferenceTrajectory*> &rts = iter->second;
		for(unsigned int i=0; i<rts.size(); i++)rt_pool.push_back(rts[i]);
		rt_by_event.erase(iter);
	}
	pthread_mutex_unlock(&rt_mutex);
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
	
	// HDDM doesn't exactly support tagged factories, but the tag
	// can be used to direct filling of the correct factory.
	string tag = factory->Tag()==NULL ? "":factory->Tag();
	
	// The ref field of the JEvent is just the s_HDDM_t pointer.
	s_HDDM_t *my_hddm_s = (s_HDDM_t*)event.GetRef();
	if(!my_hddm_s)throw RESOURCE_UNAVAILABLE;

	// Get pointer to the B-field object and Geometry object
	JEventLoop *loop = event.GetJEventLoop();
	if(loop){
		DApplication *dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
		if(dapp){
			bfield = dapp->GetBfield();
			geom = dapp->GetDGeometry(event.GetRunNumber());
		}
	}

	// Get name of data class we're trying to extract
	string dataClassName = factory->GetDataClassName();
	
	if(dataClassName =="DTagger" && tag=="")
	  return Extract_DTagger(my_hddm_s, dynamic_cast<JFactory<DTagger>*>(factory));
	
	if(dataClassName =="DMCTrackHit" && tag=="")
	  return Extract_DMCTrackHit(my_hddm_s, dynamic_cast<JFactory<DMCTrackHit>*>(factory));
	
	if(dataClassName =="DBeamPhoton" && tag=="")
	  return Extract_DBeamPhoton(my_hddm_s, dynamic_cast<JFactory<DBeamPhoton>*>(factory));
	
	if(dataClassName =="DMCThrown" && tag=="")
	  return Extract_DMCThrown(my_hddm_s, dynamic_cast<JFactory<DMCThrown>*>(factory));
	
	if(dataClassName == "DBCALTruthShower" && tag=="")
	  return Extract_DBCALTruthShower(my_hddm_s, dynamic_cast<JFactory<DBCALTruthShower>*>(factory));
	
	if(dataClassName =="DBCALHit" && tag=="")
	  return Extract_DBCALHit(my_hddm_s, dynamic_cast<JFactory<DBCALHit>*>(factory));
	
	if(dataClassName =="DCDCHit" && (tag=="" || tag=="TRUTH") )
	  return Extract_DCDCHit(my_hddm_s, dynamic_cast<JFactory<DCDCHit>*>(factory) , tag );
	
	if(dataClassName =="DFDCHit" && (tag=="" || tag=="TRUTH") )
	  return Extract_DFDCHit(my_hddm_s, dynamic_cast<JFactory<DFDCHit>*>(factory), tag);
	
	if(dataClassName == "DFCALTruthShower" && tag=="")
	  return Extract_DFCALTruthShower(my_hddm_s, dynamic_cast<JFactory<DFCALTruthShower>*>(factory));
	
	if(dataClassName == "DFCALHit" && (tag=="" || tag=="TRUTH"))
	  return Extract_DFCALHit(my_hddm_s, dynamic_cast<JFactory<DFCALHit>*>(factory), event.GetJEventLoop(), tag);
	
	if(dataClassName == "DCCALTruthShower" && tag=="")
	  return Extract_DCCALTruthShower(my_hddm_s, dynamic_cast<JFactory<DCCALTruthShower>*>(factory));
	
	if(dataClassName == "DCCALHit" && (tag=="" || tag=="TRUTH"))
	  return Extract_DCCALHit(my_hddm_s, dynamic_cast<JFactory<DCCALHit>*>(factory), event.GetJEventLoop(), tag);
	
	if(dataClassName =="DMCTrajectoryPoint" && tag=="")
	  return Extract_DMCTrajectoryPoint(my_hddm_s, dynamic_cast<JFactory<DMCTrajectoryPoint>*>(factory));
	
	if(dataClassName =="DTOFTruth" && tag=="")
	  return Extract_DTOFTruth(my_hddm_s, dynamic_cast<JFactory<DTOFTruth>*>(factory));
	
	// TOF is a special case: TWO factories are needed at the same time
	// DTOFRawHit and DTOFRawHitMC
	if(dataClassName == "DTOFRawHit" && (tag=="" || tag=="TRUTH")) {
	  JFactory_base* factory2 = loop->GetFactory("DTOFRawHitMC", tag.c_str()); 
	  return Extract_DTOFRawHit(my_hddm_s, dynamic_cast<JFactory<DTOFRawHit>*>(factory),  dynamic_cast<JFactory<DTOFRawHitMC>*>(factory2), tag);
	}
	if(dataClassName == "DTOFRawHitMC" && (tag=="" || tag=="TRUTH")) {
	  JFactory_base* factory2 = loop->GetFactory("DTOFRawHit", tag.c_str()); 
	  return Extract_DTOFRawHit(my_hddm_s, dynamic_cast<JFactory<DTOFRawHit>*>(factory2), dynamic_cast<JFactory<DTOFRawHitMC>*>(factory),  tag);
	}

	if(dataClassName =="DSCHit" && tag=="")
	  return Extract_DSCHit(my_hddm_s, dynamic_cast<JFactory<DSCHit>*>(factory));

	if(dataClassName =="DSCTruthHit" && tag=="")
	  return Extract_DSCTruthHit(my_hddm_s, dynamic_cast<JFactory<DSCTruthHit>*>(factory));

	if(dataClassName =="DTrackTimeBased" && tag=="")
	  return Extract_DTrackTimeBased(my_hddm_s, dynamic_cast<JFactory<DTrackTimeBased>*>(factory));

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
	GetSCTruthHits(hddm_s, data);

	// It has happened that some CDC hits have "nan" for the drift time
	// in a peculiar event Alex Somov came across. This ultimately caused
	// a seg. fault in MCTrackHitSort_C. I hate doing this since it
	// is treating the symptom rather than the cause, but nonetheless,
	// it patches up the problem for now until there is time to revisit
	// it later.
	for(unsigned int i=0;i<data.size(); i++)if(!finite(data[i]->z))data[i]->z=-1000.0;
	
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
		if (hits == HDDM_NULL)continue;
		if(hits->centralDC == HDDM_NULL)continue;
		if(hits->centralDC->cdcTruthPoints == HDDM_NULL)continue;
		
		s_CdcTruthPoints_t *cdctruthpoints = hits->centralDC->cdcTruthPoints;
		s_CdcTruthPoint_t *cdctruthpoint = cdctruthpoints->in;
		for(unsigned int j=0; j<cdctruthpoints->mult; j++, cdctruthpoint++){
			DMCTrackHit *mctrackhit = new DMCTrackHit;
			mctrackhit->r			= cdctruthpoint->r;
			mctrackhit->phi		= cdctruthpoint->phi;
			mctrackhit->z			= cdctruthpoint->z;
			mctrackhit->track		= cdctruthpoint->track;
			mctrackhit->primary	= cdctruthpoint->primary;
			mctrackhit->ptype		= cdctruthpoint->ptype;
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
		if(hits == HDDM_NULL)continue;
		if(hits->forwardDC == HDDM_NULL)continue;
		if(hits->forwardDC->fdcChambers == HDDM_NULL)continue;

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
				mctrackhit->ptype		= truth->ptype;
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
			mctrackhit->ptype		= bcalTruthShower->ptype;
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
			mctrackhit->ptype    = ftoftruthpoint->ptype; // save GEANT particle type 
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
			mctrackhit->ptype		= fcalTruthShower->ptype;
			mctrackhit->system	= SYS_FCAL;
			data.push_back(mctrackhit);
		}
	}

	return NOERROR;
}

//-------------------
// GetCCALTruthHits
//-------------------
jerror_t DEventSourceHDDM::GetCCALTruthHits(s_HDDM_t *hddm_s, vector<DMCTrackHit*>& data)
{
	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
#if 0	
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
			mctrackhit->ptype		= fcalTruthShower->ptype;
			mctrackhit->system	= SYS_FCAL;
			data.push_back(mctrackhit);
		}
	}
#endif

	return NOERROR;
}


//-------------------
// GetSCTruthHits
//-------------------
jerror_t DEventSourceHDDM::GetSCTruthHits(s_HDDM_t *hddm_s,  vector<DMCTrackHit*>& data)
{
	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		
		if (hits == HDDM_NULL ||
			hits->startCntr == HDDM_NULL ||
			hits->startCntr->stcTruthPoints == HDDM_NULL)continue;
		
		s_StcTruthPoints_t *stctruthpoints = hits->startCntr->stcTruthPoints;
		s_StcTruthPoint_t *stctruthpoint   = stctruthpoints->in;

		for(unsigned int j=0; j<stctruthpoints->mult; j++, stctruthpoint++){
			DMCTrackHit *mctrackhit = new DMCTrackHit;
			mctrackhit->r			=   stctruthpoint->r;
			mctrackhit->phi		        =   stctruthpoint->phi;
			mctrackhit->z			=   stctruthpoint->z;
			mctrackhit->track		=   stctruthpoint->track;
			mctrackhit->primary	        =   stctruthpoint->primary;
			mctrackhit->ptype               =   stctruthpoint->ptype;    // save GEANT particle type 
			mctrackhit->system 	        =   SYS_START;
			data.push_back(mctrackhit);
		}
	}
		
	return NOERROR;
}


//------------------
// Extract_DBCALHit
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
		s_BcalfADCCells_t *cells = hits->barrelEMcal->bcalfADCCells;
		for(unsigned int j=0;j<cells->mult;j++){
			s_BcalfADCCell_t *cell = &cells->in[j];
			int cellIdent = DBCALGeometry::cellId( cell->module, cell->layer, 
				  cell->sector );
			if(cell->bcalfADCUpHits != HDDM_NULL || cell->bcalfADCDownHits != HDDM_NULL){
				for(unsigned int k=0; k<cell->bcalfADCUpHits->mult; k++){
			     
					s_BcalfADCUpHit_t *uphit = &cell->bcalfADCUpHits->in[k];

					DBCALHit *response = new DBCALHit;
					
					response->module =cell->module;
					response->layer = cell->layer;
					response->sector = cell->sector;
					response->E = uphit->E;
					response->t = uphit->t;
					response->end = DBCALGeometry::kUpstream;
					response->cellId = cellIdent;

					data.push_back(response);
				}

				for(unsigned int k=0; k<cell->bcalfADCDownHits->mult; k++){

					s_BcalfADCDownHit_t *downhit = &cell->bcalfADCDownHits->in[k];
					DBCALHit *response = new DBCALHit;
	    
					response->module =cell->module;
					response->layer = cell->layer;
					response->sector = cell->sector;
					response->E = downhit->E;
					response->t = downhit->t;
					response->end = DBCALGeometry::kDownstream;

					response->cellId = cellIdent;

					data.push_back(response);

					//Old BCALHit code
					//DHDDMBCALHit *bcalhit = new DHDDMBCALHit();
					//bcalhit->module = cell->module;
					//bcalhit->layer = cell->layer;
					//bcalhit->sector = cell->sector;
					//bcalhit->E = hit->E;
					//bcalhit->t = hit->t;
					//bcalhit->zLocal = hit->zLocal;
					//data.push_back(bcalhit);
				}
			}
		} // j   (cells)
	} // i   (physicsEvents)
	
	// Copy into factory
	factory->CopyTo(data);

	return NOERROR;
}

//------------------
// Extract_DBeamPhoton
//------------------
jerror_t DEventSourceHDDM::Extract_DBeamPhoton(s_HDDM_t *hddm_s,  JFactory<DBeamPhoton> *factory)
{
	/// Copies the data from the given hddm_s structure. This is called
	/// from JEventSourceHDDM::GetObjects. If factory is NULL, this
	/// returns OBJECT_NOT_AVAILABLE immediately.
	
	if(factory==NULL)return OBJECT_NOT_AVAILABLE;
	
	vector<DBeamPhoton*> data;

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	for(unsigned int i=0; i<PE->mult; i++){
		// ------------ Reactions --------------
		s_Reactions_t *reactions=PE->in[i].reactions;
		if(!reactions)continue;

		for(unsigned int j=0; j<reactions->mult; j++){
			s_Beam_t *beam = reactions->in[j].beam;
			if(beam!=HDDM_NULL){
				s_Momentum_t *momentum = beam->momentum;
				
				DVector3 pos(0.0, 0.0, 65.0);
				DVector3 mom(momentum->px, momentum->py, momentum->pz);
				DBeamPhoton *photon = new DBeamPhoton;
				photon->setPosition(pos);
				photon->setMomentum(mom);
				photon->setMass(0.0);
				photon->setCharge(0.0);
				photon->clearErrorMatrix();
				photon->t = 0.0;
		
				data.push_back(photon);
			}
		}
	}
	
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

							double E  = product->momentum->E;
							double px = product->momentum->px;
							double py = product->momentum->py;
							double pz = product->momentum->pz;
							double mass = sqrt(E*E - (px*px + py*py + pz*pz));
							if(!finite(mass))mass = 0.0;
							
							DMCThrown *mcthrown = new DMCThrown;
							mcthrown->type = product->type;
							mcthrown->myid = product->id;
							mcthrown->parentid = product->parentid;
							mcthrown->mech = product->mech;
							mcthrown->pdgtype = product->pdgtype;
							mcthrown->setMass(mass);
							mcthrown->setMomentum(DVector3(px, py, pz));
							mcthrown->setPosition(DVector3(origin->vx, origin->vy, origin->vz));
							mcthrown->setCharge(ParticleCharge(product->type));
							
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
jerror_t DEventSourceHDDM::Extract_DCDCHit(s_HDDM_t *hddm_s,  JFactory<DCDCHit> *factory, string tag)
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
		
		if (hits == HDDM_NULL)continue;
		if (hits->centralDC == HDDM_NULL)continue;
		if (hits->centralDC->cdcStraws == HDDM_NULL)continue;

		
		for(unsigned int k=0; k<hits->centralDC->cdcStraws->mult; k++){
		  s_CdcStraw_t *cdcstraw = &hits->centralDC->cdcStraws->in[k];

		  if (tag=="") {
		    for(unsigned int j=0; j<cdcstraw->cdcStrawHits->mult; j++){
		      s_CdcStrawHit_t *strawhit = &cdcstraw->cdcStrawHits->in[j];
		      
		      DCDCHit *hit = new DCDCHit;
		      hit->ring    = cdcstraw->ring;
		      hit->straw   = cdcstraw->straw;
		      hit->dE      = strawhit->dE;
		      hit->t       = strawhit->t;
		      hit->d       = 0.; // initialize to zero to avoid any NaN
		      hit->itrack  = strawhit->itrack;
		      hit->ptype   = strawhit->ptype;
		      
		      data.push_back(hit);
		    }
		  } // end of tag==""

		  if (tag=="TRUTH") {
		    for(unsigned int j=0; j<cdcstraw->cdcStrawTruthHits->mult; j++){
		      s_CdcStrawTruthHit_t *strawhit = &cdcstraw->cdcStrawTruthHits->in[j];
		      
		      DCDCHit *hit = new DCDCHit;
		      hit->ring    = cdcstraw->ring;
		      hit->straw   = cdcstraw->straw;
		      hit->dE      = strawhit->dE;
		      hit->t       = strawhit->t;
		      hit->d       = strawhit->d;
		      hit->itrack  = strawhit->itrack;
		      hit->ptype   = strawhit->ptype;
		      
		      data.push_back(hit);
		    }
		  } // end of tag=="TRUTH"
		}
	}
	
	// Copy into factory
	factory->CopyTo(data);

	return NOERROR;
}


//------------------
// Extract_DFDCHit
//------------------
jerror_t DEventSourceHDDM::Extract_DFDCHit(s_HDDM_t *hddm_s,  JFactory<DFDCHit> *factory, string tag)
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
			continue;
		}

		if (hits->forwardDC == HDDM_NULL) {
		  //throw JException("HDDM source has no forwardDC information.");
			continue;
		}

		if (hits->forwardDC->fdcChambers == HDDM_NULL) {
		  // throw JException("HDDM source has no hits in the FDC.");		
			continue;
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

		  if (tag==""){
		    for (unsigned int j=0; j < wireSet->mult; j++) {

		      s_FdcAnodeWire_t anodeWire		= wireSet->in[j];
		      s_FdcAnodeHits_t* wireHitSet	        = anodeWire.fdcAnodeHits;

		      for (unsigned int k=0; k < wireHitSet->mult; k++) {

			s_FdcAnodeHit_t wireHit		        = wireHitSet->in[k];
			DFDCHit* newHit				= new DFDCHit();
			newHit->layer		 		= fdcChamber.layer;
			newHit->module		 		= fdcChamber.module;
			newHit->element				= anodeWire.wire;
			newHit->q				= wireHit.dE;
			newHit->t				= wireHit.t;
 		        newHit->d                               = 0.; // initialize to zero to avoid any NaN
			newHit->itrack                          = wireHit.itrack;
			newHit->ptype                           = wireHit.ptype;
			newHit->plane				= 2;
			newHit->type				= 0;
			newHit->gPlane				= DFDCGeometry::gPlane(newHit);
			newHit->gLayer				= DFDCGeometry::gLayer(newHit);
			newHit->r				= DFDCGeometry::getWireR(newHit);
			
			data.push_back(newHit);
		      }
		    }
		    
		    
		    // Ditto for the cathodes.
		    for (unsigned int j=0; j < stripSet->mult; j++) {

		      s_FdcCathodeStrip_t cathodeStrip          = stripSet->in[j];
		      s_FdcCathodeHits_t* stripHitSet           = cathodeStrip.fdcCathodeHits;

		      for (unsigned int k=0; k < stripHitSet->mult; k++) {

			s_FdcCathodeHit_t stripHit	        = stripHitSet->in[k];
			DFDCHit* newHit				= new DFDCHit();
			newHit->layer				= fdcChamber.layer;
			newHit->module				= fdcChamber.module;
			newHit->element				= cathodeStrip.strip;
			if (newHit->element>1000) newHit->element-=1000;
			
			newHit->plane				= cathodeStrip.plane;
			newHit->q				= stripHit.q;
			newHit->t				= stripHit.t;
		        newHit->d                               = 0.; // initialize to zero to avoid any NaN
			newHit->itrack                          = stripHit.itrack;
			newHit->ptype                           = stripHit.ptype;
			newHit->type				= 1;
			newHit->gPlane				= DFDCGeometry::gPlane(newHit);	 
			newHit->gLayer				= DFDCGeometry::gLayer(newHit);
			newHit->r				= DFDCGeometry::getStripR(newHit);
			
			data.push_back(newHit);
		      }
		    }

		  } else if (tag=="TRUTH"){
		    
		    for (unsigned int j=0; j < wireSet->mult; j++) {

		      s_FdcAnodeWire_t anodeWire		= wireSet->in[j];
		      s_FdcAnodeTruthHits_t* wireHitSet	        = anodeWire.fdcAnodeTruthHits;

		      for (unsigned int k=0; k < wireHitSet->mult; k++) {

			s_FdcAnodeTruthHit_t wireHit		= wireHitSet->in[k];
			DFDCHit* newHit				= new DFDCHit();
			newHit->layer		 		= fdcChamber.layer;
			newHit->module		 		= fdcChamber.module;
			newHit->element				= anodeWire.wire;
			newHit->q				= wireHit.dE;
			newHit->t				= wireHit.t;
			newHit->d				= wireHit.d;
			newHit->itrack                          = wireHit.itrack;
			newHit->ptype                           = wireHit.ptype;
			newHit->plane				= 2;
			newHit->type				= 0;
			newHit->gPlane				= DFDCGeometry::gPlane(newHit);
			newHit->gLayer				= DFDCGeometry::gLayer(newHit);
			newHit->r				= DFDCGeometry::getWireR(newHit);
			
			data.push_back(newHit);
		      }
		    }
		    
		    
		    // Ditto for the cathodes.
		    for (unsigned int j=0; j < stripSet->mult; j++) {

		      s_FdcCathodeStrip_t cathodeStrip          = stripSet->in[j];
		      s_FdcCathodeTruthHits_t* stripHitSet      = cathodeStrip.fdcCathodeTruthHits;

		      for (unsigned int k=0; k < stripHitSet->mult; k++) {

			s_FdcCathodeTruthHit_t stripHit	= stripHitSet->in[k];
			DFDCHit* newHit				= new DFDCHit();
			newHit->layer				= fdcChamber.layer;
			newHit->module				= fdcChamber.module;
			newHit->element				= cathodeStrip.strip;
			if (newHit->element>1000) newHit->element-=1000;
			
			newHit->plane				= cathodeStrip.plane;
			newHit->q				= stripHit.q;
			newHit->t				= stripHit.t;
	   	        newHit->d                               = 0.; // initialize to zero to avoid any NaN
			newHit->itrack                          = stripHit.itrack;
			newHit->ptype                           = stripHit.ptype;
			newHit->type				= 1;
			newHit->gPlane				= DFDCGeometry::gPlane(newHit);	 
			newHit->gLayer				= DFDCGeometry::gLayer(newHit);
			newHit->r				= DFDCGeometry::getStripR(newHit);
			
			data.push_back(newHit);
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
// Extract_DBCALTruthShower
//------------------
jerror_t DEventSourceHDDM::Extract_DBCALTruthShower(s_HDDM_t *hddm_s,  JFactory<DBCALTruthShower> *factory)
{
	/// Copies the data from the given hddm_s structure. This is called
	/// from JEventSourceHDDM::GetObjects. If factory is NULL, this
	/// returns OBJECT_NOT_AVAILABLE immediately.

	if(factory==NULL)return OBJECT_NOT_AVAILABLE;

	vector<DBCALTruthShower*> data;

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	for(unsigned int i = 0; i < PE->mult; i++) {
		s_BcalTruthShowers_t* bcalTruthShowers = NULL;
		if(PE->in[i].hitView == HDDM_NULL)continue;
		if(PE->in[i].hitView->barrelEMcal == HDDM_NULL)continue;
		if((bcalTruthShowers = PE->in[i].hitView->barrelEMcal->bcalTruthShowers) == HDDM_NULL)continue;
       
		for(unsigned int j = 0; j < bcalTruthShowers->mult; j++) {
			DBCALTruthShower *bcaltruth = new DBCALTruthShower;
			bcaltruth->track = bcalTruthShowers->in[j].track;
			bcaltruth->primary = bcalTruthShowers->in[j].primary ? 1 : 0;
			bcaltruth->phi = bcalTruthShowers->in[j].phi;
			bcaltruth->r = bcalTruthShowers->in[j].r;
			bcaltruth->z = bcalTruthShowers->in[j].z;
			bcaltruth->t = bcalTruthShowers->in[j].t;
			bcaltruth->E = bcalTruthShowers->in[j].E;
			data.push_back(bcaltruth);
		}
	}

	// Copy into factory
	factory->CopyTo(data);

	return NOERROR;
}

//------------------
// Extract_DFCALTruthShower
//------------------
jerror_t DEventSourceHDDM::Extract_DFCALTruthShower(s_HDDM_t *hddm_s,  JFactory<DFCALTruthShower> *factory)
{
  	/// Copies the data from the given hddm_s structure. This is called
	/// from JEventSourceHDDM::GetObjects. If factory is NULL, this
	/// returns OBJECT_NOT_AVAILABLE immediately.
	
	if(factory==NULL)return OBJECT_NOT_AVAILABLE;
	
	vector<DFCALTruthShower*> data;

		// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	JObject::oid_t id=1;
	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL ||
			hits->forwardEMcal == HDDM_NULL ||
			hits->forwardEMcal->fcalTruthShowers == HDDM_NULL)continue;

		s_FcalTruthShowers_t *showers = hits->forwardEMcal->fcalTruthShowers;
		for(unsigned int j=0; j<showers->mult; j++){
			s_FcalTruthShower_t *shower = &showers->in[j];
			
			DFCALTruthShower *dfcaltruthshower = new DFCALTruthShower(
				id++,
				shower->x,
				shower->y,
				shower->z,
				shower->px,
				shower->py,
				shower->pz,
				shower->E,
				shower->t,
				shower->primary,
				shower->track,
				shower->ptype
				);
			
			data.push_back(dfcaltruthshower);
		}

	} // i  (physicsEvents)

	// Copy into factory
	factory->CopyTo(data);

	return NOERROR;
}

//------------------
// Extract_DFCALHit
//------------------
jerror_t DEventSourceHDDM::Extract_DFCALHit(s_HDDM_t *hddm_s,  JFactory<DFCALHit> *factory, JEventLoop* eventLoop, string tag)
{
  /// Copies the data from the given hddm_s structure. This is called
  /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
  /// returs OBJECT_NOT_AVAILABLE immediately.
	
  if(factory==NULL)return OBJECT_NOT_AVAILABLE;

	// extract the FCAL Geometry (for isBlockActive() and positionOnFace())
	vector<const DFCALGeometry*> fcalGeomVect;
	eventLoop->Get( fcalGeomVect );
	if(fcalGeomVect.size()<1)return OBJECT_NOT_AVAILABLE;
	const DFCALGeometry& fcalGeom = *(fcalGeomVect[0]);

	vector<DFCALHit*> data;

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;

	int hitId = 0;
	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL ||
			 hits->forwardEMcal == HDDM_NULL ||
			 hits->forwardEMcal->fcalBlocks == HDDM_NULL)continue;
		
		s_FcalBlocks_t *blocks = hits->forwardEMcal->fcalBlocks;		
		for(unsigned int j=0; j<blocks->mult; j++){
			 s_FcalBlock_t *block = &blocks->in[j];
			 
			 // Filter out non-physical blocks here
			 if(!fcalGeom.isBlockActive(block->row, block->column))continue;
			 
			 // Get position of blocks on front face. (This should really come from
			 // hdgeant directly so the poisitions can be shifted in mcsmear.)
			 DVector2 pos = fcalGeom.positionOnFace(block->row, block->column);
			 
			 // Real hits
			 if(tag==""){
				 for(unsigned int k=0; k<block->fcalHits->mult; k++){
					  s_FcalHit_t *fcalhit = &block->fcalHits->in[k];
					  
					  DFCALHit *mchit = new DFCALHit();
					  mchit->row = block->row;
					  mchit->column = block->column;
					  mchit->x = pos.X();
					  mchit->y = pos.Y();
					  mchit->E = fcalhit->E;
					  mchit->t = fcalhit->t;
					  mchit->id = hitId++;
					  
					  data.push_back(mchit);
					  
				 } // k  (fcalhits)
			}else if(tag=="TRUTH"){
				 for(unsigned int k=0; k<block->fcalTruthHits->mult; k++){
					  s_FcalTruthHit_t *fcalhit = &block->fcalTruthHits->in[k];
					  
					  DFCALHit *mchit = new DFCALHit();
					  mchit->row = block->row;
					  mchit->column = block->column;
					  mchit->x = pos.X();
					  mchit->y = pos.Y();
					  mchit->E = fcalhit->E;
					  mchit->t = fcalhit->t;
					  mchit->id = hitId++;
					  
					  data.push_back(mchit);
					  
				 } // k  (fcalhits)
			} // tag
		} // j  (blocks)
	} // i  (physicsEvents)
  
  // Copy into factory
  factory->CopyTo(data);
  
  return NOERROR;
}

//------------------
// Extract_DCCALTruthShower
//------------------
jerror_t DEventSourceHDDM::Extract_DCCALTruthShower(s_HDDM_t *hddm_s,  JFactory<DCCALTruthShower> *factory)
{
  	/// Copies the data from the given hddm_s structure. This is called
	/// from JEventSourceHDDM::GetObjects. If factory is NULL, this
	/// returns OBJECT_NOT_AVAILABLE immediately.
	
	if(factory==NULL)return OBJECT_NOT_AVAILABLE;
	
	vector<DCCALTruthShower*> data;

		// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	JObject::oid_t id=1;
	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL ||
			hits->ComptonEMcal == HDDM_NULL ||
			hits->ComptonEMcal->ccalTruthShowers == HDDM_NULL)continue;

		s_CcalTruthShowers_t *showers = hits->ComptonEMcal->ccalTruthShowers;
		for(unsigned int j=0; j<showers->mult; j++){
			s_CcalTruthShower_t *shower = &showers->in[j];
			
			DCCALTruthShower *dccaltruthshower = new DCCALTruthShower(
				id++,
				shower->x,
				shower->y,
				shower->z,
				shower->px,
				shower->py,
				shower->pz,
				shower->E,
				shower->t,
				shower->primary,
				shower->track,
				shower->ptype
				);
			
			data.push_back(dccaltruthshower);
		}

	} // i  (physicsEvents)

	// Copy into factory
	factory->CopyTo(data);

	return NOERROR;
}

//------------------
// Extract_DCCALHit
//------------------
jerror_t DEventSourceHDDM::Extract_DCCALHit(s_HDDM_t *hddm_s,  JFactory<DCCALHit> *factory, JEventLoop* eventLoop, string tag)
{
  /// Copies the data from the given hddm_s structure. This is called
  /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
  /// returs OBJECT_NOT_AVAILABLE immediately.
	
  if(factory==NULL)return OBJECT_NOT_AVAILABLE;

	// extract the CCAL Geometry (for isBlockActive() and positionOnFace())
	vector<const DCCALGeometry*> ccalGeomVect;
	eventLoop->Get( ccalGeomVect );
	if(ccalGeomVect.size()<1)return OBJECT_NOT_AVAILABLE;
	const DCCALGeometry& ccalGeom = *(ccalGeomVect[0]);

	vector<DCCALHit*> data;

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;

	int hitId = 0;
	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL ||
			 hits->ComptonEMcal == HDDM_NULL ||
			 hits->ComptonEMcal->ccalBlocks == HDDM_NULL)continue;
		
		s_CcalBlocks_t *blocks = hits->ComptonEMcal->ccalBlocks;		
		for(unsigned int j=0; j<blocks->mult; j++){
			 s_CcalBlock_t *block = &blocks->in[j];
			 
			 // Filter out non-physical blocks here
			 if(!ccalGeom.isBlockActive(block->row, block->column))continue;
			 
			 // Get position of blocks on front face. (This should really come from
			 // hdgeant directly so the poisitions can be shifted in mcsmear.)
			 DVector2 pos = ccalGeom.positionOnFace(block->row, block->column);
			 
			 // Real hits
			 if(tag==""){
				 for(unsigned int k=0; k<block->ccalHits->mult; k++){
					  s_CcalHit_t *ccalhit = &block->ccalHits->in[k];
					  
					  DCCALHit *mchit = new DCCALHit();
					  mchit->row = block->row;
					  mchit->column = block->column;
					  mchit->x = pos.X();
					  mchit->y = pos.Y();
					  mchit->E = ccalhit->E;
					  mchit->t = ccalhit->t;
					  mchit->id = hitId++;
					  
					  data.push_back(mchit);
					  
				 } // k  (ccalhits)
			}else if(tag=="TRUTH"){
				 for(unsigned int k=0; k<block->ccalTruthHits->mult; k++){
					  s_CcalTruthHit_t *ccalhit = &block->ccalTruthHits->in[k];
					  
					  DCCALHit *mchit = new DCCALHit();
					  mchit->row = block->row;
					  mchit->column = block->column;
					  mchit->x = pos.X();
					  mchit->y = pos.Y();
					  mchit->E = ccalhit->E;
					  mchit->t = ccalhit->t;
					  mchit->id = hitId++;
					  
					  data.push_back(mchit);
					  
				 } // k  (ccalhits)
			} // tag
		} // j  (blocks)
	} // i  (physicsEvents)
  
  // Copy into factory
  factory->CopyTo(data);
  
  return NOERROR;
}

//------------------
// Extract_DMCTrajectoryPoint
//------------------
jerror_t DEventSourceHDDM::Extract_DMCTrajectoryPoint(s_HDDM_t *hddm_s,  JFactory<DMCTrajectoryPoint> *factory)
{
	/// Copies the data from the given hddm_s structure. This is called
	/// from JEventSourceHDDM::GetObjects. If factory is NULL, this
	/// returns OBJECT_NOT_AVAILABLE immediately.
	
	if(factory==NULL)return OBJECT_NOT_AVAILABLE;
	
	vector<DMCTrajectoryPoint*> data;

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if(hits == HDDM_NULL || hits==NULL)continue;
		if(hits->mcTrajectory == HDDM_NULL || hits->mcTrajectory==NULL)continue;
		if(hits->mcTrajectory->mcTrajectoryPoints == HDDM_NULL || hits->mcTrajectory->mcTrajectoryPoints==NULL)continue;

		s_McTrajectoryPoints_t *points = hits->mcTrajectory->mcTrajectoryPoints;
		for(unsigned int i=0; i<points->mult; i++){
			DMCTrajectoryPoint *p = new DMCTrajectoryPoint;
			
			p->x = points->in[i].x;
			p->y = points->in[i].y;
			p->z = points->in[i].z;
			p->t = points->in[i].t;
			p->px = points->in[i].px;
			p->py = points->in[i].py;
			p->pz = points->in[i].pz;
			p->E = points->in[i].E;

			p->dE = points->in[i].dE;
			p->primary_track = points->in[i].primary_track;
			p->track = points->in[i].track;
			p->part = points->in[i].part;

			p->radlen = points->in[i].radlen;
			p->step = points->in[i].step;
			p->mech = points->in[i].mech;
			
			data.push_back(p);
		}		
	}
	
	// Copy into factory
	factory->CopyTo(data);

	return NOERROR;
}

//------------------
// Extract_DTOFTruth
//------------------
jerror_t DEventSourceHDDM::Extract_DTOFTruth(s_HDDM_t *hddm_s,  JFactory<DTOFTruth>* factory)
{
  /// Copies the data from the given hddm_s structure. This is called
  /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
  /// returns OBJECT_NOT_AVAILABLE immediately.
	
  if(factory==NULL)return OBJECT_NOT_AVAILABLE;
  
  vector<DTOFTruth*> data;

  // Loop over Physics Events
  s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
  if(!PE) return NOERROR;
	
  for(unsigned int i=0; i<PE->mult; i++){
    s_HitView_t *hits = PE->in[i].hitView;
    if (hits == HDDM_NULL ||
	hits->forwardTOF == HDDM_NULL ||
	hits->forwardTOF->ftofTruthPoints == HDDM_NULL)continue;

    s_FtofTruthPoints_t* ftofTruthPoints = hits->forwardTOF->ftofTruthPoints;

    // Loop truth hits
    s_FtofTruthPoint_t *ftofTruthPoint = ftofTruthPoints->in;
    for(unsigned int j=0;j<ftofTruthPoints->mult; j++, ftofTruthPoint++){
      DTOFTruth *toftruth = new DTOFTruth;
		
      toftruth->primary     = ftofTruthPoint->primary;
      toftruth->track       = ftofTruthPoint->track;
      toftruth->x           = ftofTruthPoint->x;
      toftruth->y           = ftofTruthPoint->y;
      toftruth->z           = ftofTruthPoint->z;
      toftruth->t           = ftofTruthPoint->t;
      toftruth->px          = ftofTruthPoint->px;
      toftruth->py          = ftofTruthPoint->py;
      toftruth->pz          = ftofTruthPoint->pz;
      toftruth->E           = ftofTruthPoint->E;
      toftruth->ptype       = ftofTruthPoint->ptype;

      data.push_back(toftruth);
    }
  }

  // Copy into factory
  factory->CopyTo(data);

  return NOERROR;
}

//------------------
// Extract_DTOFRawHit
//------------------
jerror_t DEventSourceHDDM::Extract_DTOFRawHit( s_HDDM_t *hddm_s,  JFactory<DTOFRawHit>* factory, JFactory<DTOFRawHitMC> *factoryMC, string tag)
{
  /// Copies the data from the given hddm_s structure. This is called
  /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
  /// returns OBJECT_NOT_AVAILABLE immediately.
	
  if(factory==NULL)return OBJECT_NOT_AVAILABLE;

  vector<DTOFRawHit*> data;
  vector<DTOFRawHitMC*> dataMC;

  s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
  if(!PE) return NOERROR;
  
  for(unsigned int i=0; i<PE->mult; i++){
    s_HitView_t *hits = PE->in[i].hitView;
    if (hits == HDDM_NULL ||
	hits->forwardTOF == HDDM_NULL ||
	hits->forwardTOF->ftofCounters == HDDM_NULL)continue;
		
    s_FtofCounters_t* ftofCounters = hits->forwardTOF->ftofCounters;
		
    // Loop over counters
    s_FtofCounter_t *ftofCounter = ftofCounters->in;
    for(unsigned int j=0;j<ftofCounters->mult; j++, ftofCounter++){
			 

      if (tag==""){
	
	// Loop over north AND south hits
	s_FtofNorthHits_t *ftofNorthHits = ftofCounter->ftofNorthHits;
	s_FtofNorthHit_t  *ftofNorthHit =  ftofNorthHits->in;
	
	for(unsigned int k=0;k<ftofNorthHits->mult; k++, ftofNorthHit++){
	  DTOFRawHit *tofhit = new DTOFRawHit;
	  tofhit->bar	         = ftofCounter->bar;
	  tofhit->plane	         = ftofCounter->plane;
	  tofhit->lr             = 0;
	  tofhit->dE     	 = ftofNorthHit->dE;
	  tofhit->t     	 = ftofNorthHit->t;
	  data.push_back(tofhit);

	  s_FtofMCHits_t *MCHits = ftofCounter->ftofNorthTruthHits->in[k].ftofMCHits;
	  for (unsigned int j=0; j<MCHits->mult; j++){
	    DTOFRawHitMC *tofmchit = new DTOFRawHitMC;
	    tofmchit->bar      = tofhit->bar;
	    tofmchit->plane    = tofhit->plane;
	    tofmchit->lr       = tofhit->lr;
	    tofmchit->itrack   = MCHits->in[j].itrack;
	    tofmchit->ptype    = MCHits->in[j].ptype;
	    tofmchit->dist     = MCHits->in[j].dist;
	    tofmchit->x        = MCHits->in[j].x;
	    tofmchit->y        = MCHits->in[j].y;
	    tofmchit->z        = MCHits->in[j].z;
	    tofmchit->px       = MCHits->in[j].px;
	    tofmchit->py       = MCHits->in[j].py;
	    tofmchit->pz       = MCHits->in[j].pz;
	    tofmchit->E        = MCHits->in[j].E;
	    dataMC.push_back(tofmchit);
	    tofhit->AddAssociatedObject(tofmchit);
	  }
	}
	
	s_FtofSouthHits_t *ftofSouthHits = ftofCounter->ftofSouthHits;
	s_FtofSouthHit_t  *ftofSouthHit =  ftofSouthHits->in;
	
	for(unsigned int k=0;k<ftofSouthHits->mult; k++, ftofSouthHit++){
	  DTOFRawHit *tofhit = new DTOFRawHit;
	  tofhit->bar	         = ftofCounter->bar;
	  tofhit->plane	         = ftofCounter->plane;
	  tofhit->lr             = 1;
	  tofhit->dE     	 = ftofSouthHit->dE;
	  tofhit->t     	 = ftofSouthHit->t;
	  data.push_back(tofhit);

	  s_FtofMCHits_t *MCHits = ftofCounter->ftofSouthTruthHits->in[k].ftofMCHits;
	  for (unsigned int j=0; j<MCHits->mult; j++){
	    DTOFRawHitMC *tofmchit = new DTOFRawHitMC;
	    tofmchit->bar      = tofhit->bar;
	    tofmchit->plane    = tofhit->plane;
	    tofmchit->lr       = tofhit->lr;
	    tofmchit->itrack   = MCHits->in[j].itrack;
	    tofmchit->ptype    = MCHits->in[j].ptype;
	    tofmchit->dist     = MCHits->in[j].dist;
	    tofmchit->x        = MCHits->in[j].x;
	    tofmchit->y        = MCHits->in[j].y;
	    tofmchit->z        = MCHits->in[j].z;
	    tofmchit->px       = MCHits->in[j].px;
	    tofmchit->py       = MCHits->in[j].py;
	    tofmchit->pz       = MCHits->in[j].pz;
	    tofmchit->E        = MCHits->in[j].E;
	    dataMC.push_back(tofmchit);
	    tofhit->AddAssociatedObject(tofmchit);
	  }
	}
      } else if (tag=="TRUTH"){
	
	// Loop over north AND south hits
	s_FtofNorthTruthHits_t *ftofNorthTruthHits = ftofCounter->ftofNorthTruthHits;
	s_FtofNorthTruthHit_t *ftofNorthTruthHit = ftofNorthTruthHits->in;
	
	for(unsigned int k=0;k<ftofNorthTruthHits->mult; k++, ftofNorthTruthHit++){
	  DTOFRawHit *tofhit = new DTOFRawHit;
	  tofhit->bar	         = ftofCounter->bar;
	  tofhit->plane	         = ftofCounter->plane;
	  tofhit->lr	         = 0;
	  tofhit->dE    	 = ftofNorthTruthHit->dE;
	  tofhit->t     	 = ftofNorthTruthHit->t;
	  data.push_back(tofhit);

	  s_FtofMCHits_t *MCHits = ftofNorthTruthHit->ftofMCHits;
	  for (unsigned int j=0; j<MCHits->mult; j++){
	    DTOFRawHitMC *tofmchit = new DTOFRawHitMC;
	    tofmchit->bar      = tofhit->bar;
	    tofmchit->plane    = tofhit->plane;
	    tofmchit->lr       = tofhit->lr;
	    tofmchit->itrack   = MCHits->in[j].itrack;
	    tofmchit->ptype    = MCHits->in[j].ptype;
	    tofmchit->dist     = MCHits->in[j].dist;
	    tofmchit->x        = MCHits->in[j].x;
	    tofmchit->y        = MCHits->in[j].y;
	    tofmchit->z        = MCHits->in[j].z;
	    tofmchit->px       = MCHits->in[j].px;
	    tofmchit->py       = MCHits->in[j].py;
	    tofmchit->pz       = MCHits->in[j].pz;
	    tofmchit->E        = MCHits->in[j].E;
	    dataMC.push_back(tofmchit);
	    tofhit->AddAssociatedObject(tofmchit);
	  }

	}

	s_FtofSouthTruthHits_t *ftofSouthTruthHits = ftofCounter->ftofSouthTruthHits;
	s_FtofSouthTruthHit_t *ftofSouthTruthHit = ftofSouthTruthHits->in;

	for(unsigned int k=0;k<ftofSouthTruthHits->mult; k++, ftofSouthTruthHit++){
	  DTOFRawHit *tofhit = new DTOFRawHit;
	  tofhit->bar	         = ftofCounter->bar;
	  tofhit->plane	         = ftofCounter->plane;
	  tofhit->lr	         = 1;
	  tofhit->dE    	 = ftofSouthTruthHit->dE;
	  tofhit->t     	 = ftofSouthTruthHit->t;
	  data.push_back(tofhit);

	  s_FtofMCHits_t *MCHits = ftofSouthTruthHit->ftofMCHits;
	  for (unsigned int j=0; j<MCHits->mult; j++){
	    DTOFRawHitMC *tofmchit = new DTOFRawHitMC;
	    tofmchit->bar      = tofhit->bar;
	    tofmchit->plane    = tofhit->plane;
	    tofmchit->lr       = tofhit->lr;
	    tofmchit->itrack   = MCHits->in[j].itrack;
	    tofmchit->ptype    = MCHits->in[j].ptype;
	    tofmchit->dist     = MCHits->in[j].dist;
	    tofmchit->x        = MCHits->in[j].x;
	    tofmchit->y        = MCHits->in[j].y;
	    tofmchit->z        = MCHits->in[j].z;
	    tofmchit->px       = MCHits->in[j].px;
	    tofmchit->py       = MCHits->in[j].py;
	    tofmchit->pz       = MCHits->in[j].pz;
	    tofmchit->E        = MCHits->in[j].E;
	    dataMC.push_back(tofmchit);
	    tofhit->AddAssociatedObject(tofmchit);
	  }
	}
      }      
    }
  }

  // Copy into factory
  factory->CopyTo(data);
  factoryMC->CopyTo(dataMC);

  return NOERROR;
}

//------------------
// Extract_DSCHit
//------------------
jerror_t DEventSourceHDDM::Extract_DSCHit( s_HDDM_t *hddm_s,  JFactory<DSCHit>* factory)
{
  /// Copies the data from the given hddm_s structure. This is called
  /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
  /// returns OBJECT_NOT_AVAILABLE immediately.
	
  if(factory==NULL)return OBJECT_NOT_AVAILABLE;

  vector<DSCHit*> data;
  
  s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
  if(!PE) return NOERROR;
  
	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL ||
			hits->startCntr == HDDM_NULL)continue;
		
		s_StcPaddles_t* stcPaddles = hits->startCntr->stcPaddles;
		if(stcPaddles==HDDM_NULL)continue;
		//s_StcTruthPoints_t* stcTruthPoints = hits->startCntr->stcTruthPoints;
		
		// Loop over counters
		s_StcPaddle_t *stcPaddle = stcPaddles->in;
		for(unsigned int j=0;j<stcPaddles->mult; j++, stcPaddle++){
			s_StcHits_t *stcHits = stcPaddle->stcHits;
			if(stcHits==HDDM_NULL)continue;
			
			s_StcHit_t *stchit = stcHits->in;
			for(unsigned int k=0; k<stcHits->mult; k++, stchit++){
				DSCHit *hit = new DSCHit;
				hit->sector = stcPaddle->sector;
				hit->dE = stchit->dE;
				hit->t = stchit->t;
				data.push_back(hit);
			}
		}
	}

  // Copy into factory
  factory->CopyTo(data);

  return NOERROR;
}


//------------------
// Extract_DSCTruthHit
//------------------
jerror_t DEventSourceHDDM::Extract_DSCTruthHit( s_HDDM_t *hddm_s,  JFactory<DSCTruthHit>* factory)
{
  /// Copies the data from the given hddm_s structure. This is called
  /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
  /// returns OBJECT_NOT_AVAILABLE immediately.
	
  if(factory==NULL)return OBJECT_NOT_AVAILABLE;

  vector<DSCTruthHit*> data;
  
  s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
  if(!PE) return NOERROR;
  
	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL ||
			hits->startCntr == HDDM_NULL)continue;
		
		s_StcTruthPoints_t* stcTruthPoints = hits->startCntr->stcTruthPoints;
		if(stcTruthPoints==HDDM_NULL)continue;
		
		// Loop over counters
		s_StcTruthPoint_t *stcTruthPoint = stcTruthPoints->in;
		for(unsigned int j=0;j<stcTruthPoints->mult; j++, stcTruthPoint++){

			DSCTruthHit *hit = new DSCTruthHit;
			hit->dEdx = stcTruthPoint->dEdx;
			hit->phi = stcTruthPoint->phi;
			hit->primary = stcTruthPoint->primary;
			hit->ptype = stcTruthPoint->ptype;
                        hit->r = stcTruthPoint->r;
			hit->t = stcTruthPoint->t;
			hit->z = stcTruthPoint->z;
			hit->track = stcTruthPoint->track;
			hit->sector = stcTruthPoint->sector;
			data.push_back(hit);
		}
	}

  // Copy into factory
  factory->CopyTo(data);

  return NOERROR;
}

//------------------
// Extract_DTrackTimeBased
//------------------
jerror_t DEventSourceHDDM::Extract_DTrackTimeBased(s_HDDM_t *hddm_s,  JFactory<DTrackTimeBased> *factory)
{
	// Note: Since this is a reconstructed factory, we want to generally return OBJECT_NOT_AVAILABLE
	// rather than NOERROR. The reason being that the caller interprets "NOERROR" to mean "yes I
	// usually can provide objects of that type, but this event has none." This will cause it to
	// skip any attempt at reconstruction. On the other hand, a value of "OBJECT_NOT_AVAILABLE" tells
	// it "I cannot provide those type of objects for this event.

  if(factory==NULL)return OBJECT_NOT_AVAILABLE;

  vector<DTrackTimeBased*> data;
  vector<DReferenceTrajectory*> rts;
  
  // Loop over physics events.
  s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
  if(!PE) return OBJECT_NOT_AVAILABLE;
  
	bool event_had_tracktimebaseds = false;
	for(unsigned int i=0; i<PE->mult; i++){
		s_ReconView_t *recon = PE->in[i].reconView;
		if (recon == HDDM_NULL || recon->tracktimebaseds == HDDM_NULL)continue;

		// Get enough DReferenceTrajectory objects for all of the DTrackTimeBased Objects
		// we're about to read in. This seems a little complicated, but that's because it
		// is expensive to allocate these things so we recycle as much as possible.
		list<DReferenceTrajectory*> my_rts;
		pthread_mutex_lock(&rt_mutex);
		while(my_rts.size() < recon->tracktimebaseds->mult){
			if(rt_pool.size()>0){
				my_rts.push_back(rt_pool.back());
				rt_pool.pop_back();
			}else{
				my_rts.push_back(new DReferenceTrajectory(bfield));
			}
		}
		pthread_mutex_unlock(&rt_mutex);
		
		// Loop over timebased tracks
		event_had_tracktimebaseds = true;
		for(unsigned int i=0; i< recon->tracktimebaseds->mult; i++){
			s_Tracktimebased_t *tbt = &(recon->tracktimebaseds->in[i]);
			
			DVector3 mom(tbt->momentum->px, tbt->momentum->py, tbt->momentum->pz);
			DVector3 pos(tbt->origin->vx, tbt->origin->vy, tbt->origin->vz);
			
			DTrackTimeBased *track = new DTrackTimeBased();

			track->setMomentum(mom);
			track->setPosition(pos);
			track->setCharge(tbt->properties->charge);
			track->setMass(tbt->properties->mass);
			track->chisq = tbt->chisq;
			track->Ndof = tbt->Ndof;
			track->FOM = tbt->FOM;
			track->candidateid = tbt->candidateid;
			track->id = tbt->id;
			
			// Reconstitute errorMatrix
			string str_vals = tbt->errorMatrix->vals;
			DMatrixDSym errMatrix;
			StringToDMatrixDSym(str_vals, errMatrix, tbt->errorMatrix->Nrows, tbt->errorMatrix->Ncols);
			track->setErrorMatrix(errMatrix);
			
			// Reconstitute TrackingErrorMatrix
			str_vals = tbt->TrackingErrorMatrix->vals;
			DMatrixDSym TrackingErrorMatrix;
			StringToDMatrixDSym(str_vals, TrackingErrorMatrix, tbt->TrackingErrorMatrix->Nrows, tbt->TrackingErrorMatrix->Ncols);
			track->setTrackingErrorMatrix(TrackingErrorMatrix);
			
			// Use DReferenceTrajectory objects (either recycled or new)
			DReferenceTrajectory *rt = my_rts.back();
			my_rts.pop_back();
			if(rt){
				rt->SetMass(track->mass());
				rt->SetDGeometry(geom);
				rt->Swim(pos, mom, track->charge());
				rts.push_back(rt);
			}
			track->rt = rt;

			data.push_back(track);
		}
	}

	// Copy into factory
	if(event_had_tracktimebaseds){
		factory->CopyTo(data);
		
		// Add DReferenceTrajectory objects to rt_by_event so they can be deleted later.
		// The rt_by_event maintains lists indexed by the hddm_s pointer since multiple
		// threads may be calling us. Note that we first look to see if a list already
		// exists for this event and append to it if it does. This is so we can the
		// same list for all objects that use DReferenceTrajectories.
		pthread_mutex_lock(&rt_mutex);
		map<s_HDDM_t*, vector<DReferenceTrajectory*> >::iterator iter = rt_by_event.find(hddm_s);
		if(iter != rt_by_event.end()){
			vector<DReferenceTrajectory*> &my_rts = iter->second;
			my_rts.insert(my_rts.end(), rts.begin(), rts.end());
		}else{
			rt_by_event[hddm_s] = rts;
		}
		pthread_mutex_unlock(&rt_mutex);

		// If the event had a s_Tracktimebased_t pointer, then report back that
		// we read them in from the file. Otherwise, report OBJECT_NOT_AVAILABLE
		return NOERROR;
	}

	// If we get to here then there was not even a placeholder in the HDDM file.
	// Return OBJECT_NOT_AVAILABLE to indicate reconstruction should be tried.
	return OBJECT_NOT_AVAILABLE;
}


//-------------------------------
// StringToDMatrixDSym
//-------------------------------
string DEventSourceHDDM::StringToDMatrixDSym(string &str_vals, DMatrixDSym &mat, int &Nrows, int Ncols)
{
	/// This is the inverse of the DMatrixDSymToString method in the
	/// danahddm plugin.

	// Convert the given string into a symmetric matrix
	mat.ResizeTo(Nrows, Ncols);
	stringstream ss(str_vals);
	for(int irow=0; irow<mat.GetNrows(); irow++) {
		for(int icol=irow; icol<mat.GetNcols(); icol++) {
			ss >> mat[irow][icol];
			mat[icol][irow] = mat[irow][icol];
		}
	}
	
	return ss.str();
}

//------------------
// Extract_DTagger
//------------------
jerror_t DEventSourceHDDM::Extract_DTagger( s_HDDM_t *hddm_s,  JFactory<DTagger>* factory)
{
  /// Copies the data from the given hddm_s structure. This is called
  /// from JEventSourceHDDM::GetObjects. If factory is NULL, this
  /// returns OBJECT_NOT_AVAILABLE immediately.

  if(factory==NULL)return OBJECT_NOT_AVAILABLE;

  vector<DTagger*> data;
  
  s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
  if(!PE) return NOERROR;
  
  for(unsigned int i=0; i<PE->mult; i++){
    s_HitView_t *hits = PE->in[i].hitView;
    if (hits == HDDM_NULL ||
	hits->tagger == HDDM_NULL)continue;
    s_MicroChannels_t *microChannels=hits->tagger->microChannels;
    if (microChannels==HDDM_NULL) continue;

    s_MicroChannel_t *microChannel=microChannels->in;
    for (unsigned int j=0;j<microChannels->mult;j++,microChannel++){
      s_TaggerHits_t *taggerHits=microChannel->taggerHits;
      if (taggerHits==HDDM_NULL) continue;

      s_TaggerHit_t *taggerHit=taggerHits->in;
      for (unsigned int k=0;k<taggerHits->mult;k++,taggerHit++){
	DTagger *tagger=new DTagger();

	tagger->E=microChannel->E;
	tagger->t=taggerHit->t;
	tagger->row=microChannel->row;
	tagger->column=microChannel->column;

	data.push_back(tagger);
      }
    }
  }
  // Copy into factory
  factory->CopyTo(data);

  return NOERROR;

}
