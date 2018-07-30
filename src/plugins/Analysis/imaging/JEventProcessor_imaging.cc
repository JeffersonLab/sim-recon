// $Id$
//
//    File: JEventProcessor_imaging.cc
// Created: Thu Nov  9 10:49:12 EST 2017
// Creator: staylor (on Linux ifarm1402.jlab.org 3.10.0-327.el7.x86_64 x86_64)
//

#include "JEventProcessor_imaging.h"
using namespace jana;
#include <TDirectory.h>

#include <PID/DChargedTrack.h>

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_imaging());
}
} // "C"


//------------------
// JEventProcessor_imaging (Constructor)
//------------------
JEventProcessor_imaging::JEventProcessor_imaging()
{

}

//------------------
// ~JEventProcessor_imaging (Destructor)
//------------------
JEventProcessor_imaging::~JEventProcessor_imaging()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_imaging::init(void)
{
  gDirectory->mkdir("Vertexes")->cd();

  TwoTrackXYZ= new TH3I("TwoTrackXYZ","z vs y vs x",400,-10,10,
			400,-10,10,140,30,100);
  TwoTrackXYZ->SetXTitle("x (cm)");
  TwoTrackXYZ->SetYTitle("y (cm)");
  TwoTrackXYZ->SetZTitle("z (cm)");
  
  TwoTrackDoca=new TH1F("TwoTrackDoca","#Deltad between tracks",1000,0,20);
  TwoTrackDoca->SetXTitle("doca [cm]"); 

  TwoTrackZ=new TH1F("TwoTrackZ","z for r<0.5 cm",1000,0,200);
  TwoTrackZ->SetXTitle("z [cm]");

  TwoTrackPocaCut=new TH2F("TwoTrackPocaCut","2track POCA,doca cut",4000,0,400,650,0,65);
  TwoTrackPocaCut->SetXTitle("z (cm)");
  TwoTrackPocaCut->SetYTitle("r (cm)"); 

  TwoTrackXY_at_65cm=new TH2F("TwoTrackXY_at_65cm","y vs x near 65 cm",400,-5,5,400,-5,5);
  TwoTrackXY_at_65cm->SetXTitle("x [cm]");
  TwoTrackXY_at_65cm->SetYTitle("y [cm]");
  
  //  TwoTrackRelCosTheta=new TH1F("TwoTrackRelCosTheta","relative direction",100,-1.,1.);
  TwoTrackChi2=new TH1F("TwoTrackChi2","vertex #chi^2",1000,0,1000);
  DocaPull=new TH1F("DocaPull","#deltad/#sigma(#deltad)",100,0.,5);
  TwoTrackProb=new TH1F("TwoTrackProb","vertex probability",100,0,1.);

  TwoTrackZFit=new TH1F("TwoTrackZFit","z for r<0.5 cm",1000,0,200);
  TwoTrackZFit->SetXTitle("z [cm]");

  TwoTrackPocaCutFit=new TH2F("TwoTrackPocaCutFit","2track POCA,doca cut",4000,0,400,650,0,65);
  TwoTrackPocaCutFit->SetXTitle("z (cm)");
  TwoTrackPocaCutFit->SetYTitle("r (cm)"); 

  TwoTrackXYFit_at_65cm=new TH2F("TwoTrackXYFit_at_65cm","y vs x near 65 cm",400,-5,5,400,-5,5);
  TwoTrackXYFit_at_65cm->SetXTitle("x [cm]");
  TwoTrackXYFit_at_65cm->SetYTitle("y [cm]"); 

  TwoTrackXYZFit= new TH3I("TwoTrackXYZFit","z vs y vs x",400,-10,10,
			   400,-10,10,140,30,100);
  TwoTrackXYZFit->SetXTitle("x (cm)");
  TwoTrackXYZFit->SetYTitle("y (cm)");
  TwoTrackXYZFit->SetZTitle("z (cm)");
  
  gDirectory->cd("../");

  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_imaging::brun(JEventLoop *eventLoop, int32_t runnumber)
{
  // This is called whenever the run number changes
  DApplication* dapp=dynamic_cast<DApplication*>(eventLoop->GetJApplication());
  geom = dapp->GetDGeometry(runnumber);
  bfield=dapp->GetBfield(runnumber);
  
  //Pre-allocate memory for DReferenceTrajectory objects early
  //The swim-step objects of these arrays take up a significant amount of memory, and it can be difficult to find enough free contiguous space for them.
  //Therefore, allocate them at the beginning before the available memory becomes randomly populated
  while(rtv.size() < 50)
    rtv.push_back(new DReferenceTrajectory(bfield));

  

  FIT_VERTEX=true;
  gPARMS->SetDefaultParameter("IMAGING:FIT_VERTEX",FIT_VERTEX, "Turn on/off vertex fitting");
  FIT_CL_CUT=0.01;
  gPARMS->SetDefaultParameter("IMAGING:FIT_CL_CUT",FIT_CL_CUT, "CL cut for vertex fit"); 
  TRACK_CL_CUT=1e-4;
  gPARMS->SetDefaultParameter("IMAGING:TRACK_CL_CUT",TRACK_CL_CUT, "CL cut for tracks");
  DOCA_CUT=1.0; 
  gPARMS->SetDefaultParameter("IMAGING:DOCA_CUT",DOCA_CUT, "Maximum doca between tracks");

  
  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_imaging::evnt(JEventLoop *loop, uint64_t eventnumber)
{
  // This is called for every event. Use of common resources like writing
  // to a file or filling a histogram should be mutex protected. Using
  // loop->Get(...) to get reconstructed objects (and thereby activating the
  // reconstruction algorithm) should be done outside of any mutex lock
  // since multiple threads may call this method at the same time.
  // Here's an example:
  //
  // vector<const MyDataClass*> mydataclasses;
  // loop->Get(mydataclasses);
  //
  // japp->RootFillLock(this);
  //  ... fill historgrams or trees ...
  // japp->RootFillUnLock(this);
  
  vector<const DChargedTrack*>tracks; 
  loop->Get(tracks); 

  japp->RootWriteLock();

  // Reset the number of used reference trajectories from the pool
  num_used_rts=0;
 
  for (unsigned int i=0;i<tracks.size();i++){
    const DTrackTimeBased *track1=tracks[i]->Get_BestTrackingFOM()->Get_TrackTimeBased();
    if (TMath::Prob(track1->chisq,track1->Ndof)>TRACK_CL_CUT){
      for (unsigned int j=i+1;j<tracks.size();j++){
	const DTrackTimeBased *track2=tracks[j]->Get_BestTrackingFOM()->Get_TrackTimeBased();

	if (TMath::Prob(track2->chisq,track2->Ndof)>TRACK_CL_CUT){
	  // Make sure there are enough DReferenceTrajectory objects
	  unsigned int locNumInitialReferenceTrajectories = rtv.size();
	  while(rtv.size()<=num_used_rts){
               //printf("Adding %d\n",rtv.size());
	    rtv.push_back(new DReferenceTrajectory(bfield));
	  }
	  DReferenceTrajectory *rt1 = rtv[num_used_rts];
	  if(locNumInitialReferenceTrajectories == rtv.size()) //didn't create a new one
	    rt1->Reset();
	  rt1->SetDGeometry(geom);
	  rt1->SetMass(track1->mass());
	  //rt1->SetStepSize(0.1);
	  rt1->FastSwim(track1->position(),track1->momentum(),track1->charge(),
		    2000.0,0.,370.);
	  num_used_rts++;

	  locNumInitialReferenceTrajectories = rtv.size();
	  while(rtv.size()<=num_used_rts){
               //printf("Adding %d\n",rtv.size());
	    rtv.push_back(new DReferenceTrajectory(bfield));
	  }
	  DReferenceTrajectory *rt2 = rtv[num_used_rts];
	  if(locNumInitialReferenceTrajectories == rtv.size()) //didn't create a new one
	    rt2->Reset();
	  rt2->SetDGeometry(geom);
	  rt2->SetMass(track2->mass());
	  //rt2->SetStepSize(0.1);
	  rt2->FastSwim(track2->position(),track2->momentum(),track2->charge(),
		    2000.0,0.,370.);
	  num_used_rts++;

	  DVector3 pos;
	  double doca,var_doca,vertex_chi2,vertex_prob=1.;
	  DKinematicData kd1=*track1,kd2=*track2;
	  rt1->IntersectTracks(rt2,&kd1,&kd2,pos,doca,var_doca,vertex_chi2,FIT_VERTEX);
	  // rt1->IntersectTracks(rt2,NULL,NULL,pos,doca,var_doca,vertex_chi2);
	  TwoTrackDoca->Fill(doca);
	  DocaPull->Fill(doca/sqrt(var_doca));
	  if (doca<DOCA_CUT){
	    TwoTrackPocaCut->Fill(pos.z(),pos.Perp());
	    TwoTrackXYZ->Fill(pos.x(),pos.y(),pos.z());
	    if (pos.z()>64.5 && pos.z()<65.5){
	      TwoTrackXY_at_65cm->Fill(pos.x(),pos.y());
	    }
	    if (pos.Perp()<0.5){
	      TwoTrackZ->Fill(pos.z());
	    }

	    if (FIT_VERTEX){
	      TwoTrackChi2->Fill(vertex_chi2);
	      vertex_prob=TMath::Prob(vertex_chi2,1);
	      TwoTrackProb->Fill(vertex_prob);
	      if (vertex_prob>FIT_CL_CUT){  
		TwoTrackPocaCutFit->Fill(pos.z(),pos.Perp());
		TwoTrackXYZFit->Fill(pos.x(),pos.y(),pos.z());
		if (pos.z()>64.5 && pos.z()<65.5){
		  TwoTrackXYFit_at_65cm->Fill(pos.x(),pos.y());
		}
		if (pos.Perp()<0.5){
		  TwoTrackZFit->Fill(pos.z());
		}
	      }
	    }
	  }
	}
      }
    }
  }

  japp->RootUnLock();


  
  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_imaging::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_imaging::fini(void)
{
	// Called before program exit after event processing is finished.
	return NOERROR;
}

