#include "JEventProcessor_FCAL_TimingOffsets.h"
#include <JANA/JApplication.h>
#include "FCAL/DFCALShower.h"
#include "FCAL/DFCALGeometry.h"
#include "FCAL/DFCALHit.h"
#include "FCAL/DFCALDigiHit.h"
#include "FCAL/DFCALCluster.h"
#include "TRACKING/DTrackTimeBased.h"
#include "PID/DVertex.h"
#include "DVector3.h"
#include "ANALYSIS/DAnalysisUtilities.h"
#include <TTree.h>
#include "DVector3.h"
#include "PID/DParticleID.h"
#include "GlueX.h"
#include <vector>
#include <map>
#include <deque>
#include <string>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <TDirectory.h>
#include <TH1I.h>
#include <TH2F.h>
#include <TProfile.h>
#include "ANALYSIS/DTreeInterface.h"
#include <thread>

using namespace jana;

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>

const int nChan = 2800;

// Define Histograms
static TH1I* Offsets[nChan];

extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_FCAL_TimingOffsets());
  }
} // "C"

//------------------
// JEventProcessor_FCAL_TimingOffsets (Constructor)
//------------------
JEventProcessor_FCAL_TimingOffsets::JEventProcessor_FCAL_TimingOffsets()
{

}

//------------------
// ~JEventProcessor_FCAL_TimingOffsets (Destructor)
//------------------
JEventProcessor_FCAL_TimingOffsets::~JEventProcessor_FCAL_TimingOffsets()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_FCAL_TimingOffsets::init(void)
{
  // This is called once at program startup. If you are creating
  // and filling historgrams in this plugin, you should lock the
  // ROOT mutex like this:
  //
	TDirectory *main = gDirectory;
  gDirectory->mkdir("FCAL_TimingOffsets")->cd();
  
  
  for (int i = 0; i < nChan; ++i) {
    Offsets[i] = new TH1I(Form("Offset_%i",i),Form("Timing Offset for Channel %i",i),800,-50,50);
  }

  main->cd();
  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_FCAL_TimingOffsets::brun(JEventLoop *eventLoop, 
					     int32_t runnumber)
{

  // get the FCAL z position from the global geometry interface
  DApplication *dapp = 
    dynamic_cast<DApplication*>(eventLoop->GetJApplication());
  const DGeometry *geom = dapp->GetDGeometry(runnumber);
  if( geom ) {

    geom->GetFCALZ( m_FCALfront );
  }
  else{
      
    cerr << "No geometry accessbile." << endl;
    return RESOURCE_UNAVAILABLE;
  }


  // we need an FCAL Geometry object
  vector< const DFCALGeometry* > geomVec;
  eventLoop->Get( geomVec );

  if( geomVec.size() != 1 ){

    cerr << "No geometry accessbile." << endl;
    return RESOURCE_UNAVAILABLE;
  }

  m_fcalGeom = geomVec[0];

  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_FCAL_TimingOffsets::evnt(JEventLoop *eventLoop, 
					     uint64_t eventnumber)
{
double FCAL_C_EFFECTIVE = 15.0;
DApplication* locApplication = dynamic_cast<DApplication*>(eventLoop->GetJApplication());
    DGeometry* locGeometry = locApplication->GetDGeometry(eventLoop->GetJEvent().GetRunNumber());
    double locTargetZCenter = 0.0;
    locTargetZCenter = locGeometry->GetTargetZ(locTargetZCenter);
    dTargetCenter.SetXYZ(0.0, 0.0, locTargetZCenter);
 
 vector< const DFCALGeometry* > geomVec;
  vector< const DFCALDigiHit*  > digiHits;
  vector< const DFCALHit*      > hits;
 vector<const DEventRFBunch*> locEventRFBunches;
 vector< const DVertex* > kinfitVertex;
 vector< const DTrackTimeBased* > locTrackTimeBased;
 vector < const DFCALShower * > matchedShowers;
  

  eventLoop->Get( geomVec );
  eventLoop->Get( hits );
  
  vector< const DFCALShower* > locFCALShowers;
 
   if( hits.size() < 500 ){  // only form clusters and showers if there aren't too many hits
    eventLoop->Get(locFCALShowers);    
    eventLoop->Get(locEventRFBunches);
    eventLoop->Get(locTrackTimeBased);
   }

  double locStartTime = locEventRFBunches.empty() ? 0.0 : locEventRFBunches[0]->dTime;

 DVector3 norm(0.0,0.0,-1);
  DVector3 pos,mom;
  
 
for (unsigned int i=0; i < locTrackTimeBased.size() ; ++i){
    for (unsigned int j=0; j< locFCALShowers.size(); ++j){

      Double_t x = locFCALShowers[j]->getPosition().X();
      Double_t y = locFCALShowers[j]->getPosition().Y();
      
     
      DVector3 pos_FCAL(0,0,638);

      if (locTrackTimeBased[i]->rt->GetIntersectionWithPlane(pos_FCAL,norm,pos,mom,NULL,NULL,NULL,SYS_FCAL)==NOERROR)
	{
	  Double_t trkmass = locTrackTimeBased[i]->mass();
	  Double_t FOM = TMath::Prob(locTrackTimeBased[i]->chisq, locTrackTimeBased[i]->Ndof);
	  Double_t dRho = sqrt(((pos.X() - x)*(pos.X() - x)) + ((pos.Y() - y)* (pos.Y() - y)));
	  if(trkmass < 0.15 && dRho < 5 && FOM > 0.0001 ) {  
	    matchedShowers.push_back(locFCALShowers[j]);
	  }
	}

    }
  }


   
  

    for(unsigned int k=0; k<locFCALShowers.size(); k++)
      {
     if (find(matchedShowers.begin(), matchedShowers.end(),locFCALShowers[k]) != matchedShowers.end()) continue;
	
const DFCALShower *s1 = locFCALShowers[k];

     
	double pos_corrected_Z = locFCALShowers[k]->getPosition().Z();

        // Get the clusters from the showers
        vector <const DFCALCluster *> clusterVector;
        s1->Get(clusterVector);
        
        // Loop over clusters within the shower
        for (unsigned int iCluster = 0; iCluster < clusterVector.size(); iCluster++){
            // Get the hits
            const vector<DFCALCluster::DFCALClusterHit_t> hitVector = clusterVector[iCluster]->GetHits();

            //Loop over hits
            for (unsigned int iHit = 0; iHit < 1; iHit++){ 
		double hitEnergy = hitVector[iHit].E;
		if (hitEnergy <= 0.4) continue;                
		double hitTime = hitVector[iHit].t;
		double tCorr = ( m_FCALfront + DFCALGeometry::blockLength() - pos_corrected_Z )/FCAL_C_EFFECTIVE;
                hitTime -= tCorr; // Apply the t corection used for the cluster/shower conversion 
                
		int chanx =  hitVector[iHit].x;
		int chany =  hitVector[iHit].y;
		int ChannelNumber = hitVector[iHit].ch;
		
		dFCALblockcenter.SetXYZ(chanx, chany,  pos_corrected_Z);

		double locPathLength = (dFCALblockcenter - dTargetCenter).Mag();
		double locDeltaT = hitTime - locPathLength/29.9792458 - locStartTime;

		Offsets[ChannelNumber]->Fill(locDeltaT);
  
            }
        }
    }
	








return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_FCAL_TimingOffsets::erun(void)
{
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_FCAL_TimingOffsets::fini(void)
{

  // Called before program exit after event processing is finished.   
  return NOERROR;
}


