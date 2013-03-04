// $Id$
//
//    File: DCDCTrackHit_factory.cc
// Created: Mon Oct 16 10:20:07 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.7.0 powerpc)
//

#include <cmath>
#include <pthread.h>
using namespace std;

#include "DCDCTrackHit_factory.h"
#include "DCDCHit.h"
#include <HDGEOMETRY/DGeometry.h>
#include <TRACKING/DTrackHitSelectorTHROWN.h>
#include <TRACKING/DMCTrackHit.h>
 
DCDCTrackHit_factory::~DCDCTrackHit_factory(){
  if (cdcwires.size()){
    for (unsigned int i=0;i<cdcwires.size();i++){
      for (unsigned int j=0;j<cdcwires[i].size();j++){
	delete cdcwires[i][j];
      }
    }    
  }
}

//------------------
// init
//------------------
jerror_t DCDCTrackHit_factory::init(void)
{
   MATCH_TRUTH_HITS=false;
  
   gPARMS->SetDefaultParameter("CDC:MATCH_TRUTH_HITS",MATCH_TRUTH_HITS,"Match truth hits to CDC hits (DEF=false)");

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DCDCTrackHit_factory::brun(JEventLoop *loop, int runnumber)
{
  // Get pointer to DGeometry object
  DApplication* dapp=dynamic_cast<DApplication*>(eventLoop->GetJApplication());
  dgeom  = dapp->GetDGeometry(runnumber);
  
  // Get the CDC wire table from the XML
  //jout<< "Getting map of cdc wires from the XML" <<endl;
  dgeom->GetCDCWires(cdcwires);
  
  // Fill array with the number of straws for each layer
  for (unsigned int i=0;i<cdcwires.size();i++){
    Nstraws[i]=cdcwires[i].size();
  }

  return NOERROR;
}

jerror_t DCDCTrackHit_factory::erun(void){
  if (cdcwires.size()){
    for (unsigned int i=0;i<cdcwires.size();i++){
      for (unsigned int j=0;j<cdcwires[i].size();j++){
	delete cdcwires[i][j];
      }
    }    
  }
  cdcwires.clear();
  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DCDCTrackHit_factory::evnt(JEventLoop *loop, int eventnumber)
{
	/// Convert from ring/straw indexing to x/y position
	/// of wire center and stereo angle.
	vector<const DCDCHit*> cdchits;
	loop->Get(cdchits);
	if (cdchits.size()==0) return NOERROR;
	
	// If this is simulated data then we want to match up the truth hit
	// with this "real" hit. Ideally, this would be done at the
	// DCDCHit object level, but the organization of the data in HDDM
	// makes that difficult. Here we have the full wire definition so
	// we make the connection here.
	vector<const DMCTrackHit*> mctrackhits;
	if (MATCH_TRUTH_HITS)loop->Get(mctrackhits);
	
	for(unsigned int i=0; i<cdchits.size(); i++){
		const DCDCHit* cdchit = cdchits[i];

		if(cdchit->ring>CDC_MAX_RINGS || cdchit->ring<1
			|| cdchit->straw>Nstraws[cdchit->ring-1] || cdchit->straw<1){
			cerr<<__FILE__<<":"<<__LINE__<<" Ring or straw out of range! ring="
				<<cdchit->ring<<" (should be 1-"<<CDC_MAX_RINGS<<")  straw="
				<<cdchit->straw<<" (should be 1-"<<Nstraws[cdchit->ring-1]<<")"<<endl;
			continue;
		}

		DCDCTrackHit *hit = new DCDCTrackHit;
		hit->wire = cdcwires[cdchit->ring-1][cdchit->straw-1];
		hit->is_stereo=((cdchit->ring>4&&cdchit->ring<13)
				||(cdchit->ring>16&&cdchit->ring<25))
				?true:false;
		hit->tdrift = cdchit->t;
		double w_eff=29.5e-9;
		double gas_gain=1e5;
		double electron_charge=1.6022e-4; /* fC */
		hit->dE=cdchit->q*w_eff/(gas_gain*electron_charge);
		hit->dist = hit->tdrift*55.0E-4; // Use number hardwired in simulation for now
		hit->AddAssociatedObject(cdchit);

		if (MATCH_TRUTH_HITS==true&&mctrackhits.size()>0){
		  // Try matching truth hit with this "real" hit.
		  const DMCTrackHit *mctrackhit = DTrackHitSelectorTHROWN::GetMCTrackHit(hit->wire, hit->dist, mctrackhits);
		
		  if(mctrackhit)hit->AddAssociatedObject(mctrackhit);
		}

		_data.push_back(hit);
	}
	
	return NOERROR;
}

