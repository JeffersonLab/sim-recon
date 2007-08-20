// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.cc
//

#include <iostream>
#include <vector>
#include <string>
using namespace std;

#include <TLorentzVector.h>
#include <TLorentzRotation.h>
#include <TEllipse.h>
#include <TArc.h>
#include <TBox.h>
#include <TLine.h>
#include <TText.h>
#include <TVector3.h>
#include <TColor.h>

#include "hdview2.h"
#include "hdv_mainframe.h"
#include "MyProcessor.h"
#include "TRACKING/DTrackHit_factory.h"
#include "TRACKING/DQuickFit.h"
#include "TRACKING/DMagneticFieldStepper.h"
#include "TRACKING/DTrackCandidate_factory.h"
#include "TRACKING/DMCTrackHit.h"
#include "TRACKING/DMCThrown.h"
#include "TRACKING/DTrack.h"
#include "TRACKING/DReferenceTrajectory.h"
#include "JANA/JGeometry.h"
#include "TRACKING/DMCTrajectoryPoint.h"
#include "FCAL/DFCALHit.h"
#include "FDC/DFDCGeometry.h"
#include "CDC/DCDCTrackHit.h"
#include "FDC/DFDCPseudo.h"

extern hdv_mainframe *hdvmf;


MyProcessor *gMYPROC=NULL;

//------------------------------------------------------------------
// MyProcessor 
//------------------------------------------------------------------
MyProcessor::MyProcessor()
{
	Bfield = NULL;
	loop = NULL;

	// Tell factory to keep around a few density histos
	//gPARMS->SetParameter("TRKFIND:MAX_DEBUG_BUFFERS",	16);
	
	gMYPROC = this;
}

//------------------------------------------------------------------
// ~MyProcessor 
//------------------------------------------------------------------
MyProcessor::~MyProcessor()
{

}

//------------------------------------------------------------------
// init 
//------------------------------------------------------------------
jerror_t MyProcessor::init(void)
{
	// Make sure detectors have been drawn
	//if(!drew_detectors)DrawDetectors();
	
	vector<JEventLoop*> loops = app->GetJEventLoops();
	if(loops.size()>0){
		vector<string> facnames;
		loops[0]->GetFactoryNames(facnames);
		
		hdvmf = new hdv_mainframe(gClient->GetRoot(), 1000, 600);
		hdvmf->SetTrackFactories(facnames);
		hdvmf->SetCandidateFactories(facnames);
		hdvmf->SetReconstructedFactories(facnames);
	}
	
	return NOERROR;
}

//------------------------------------------------------------------
// brun 
//------------------------------------------------------------------
jerror_t MyProcessor::brun(JEventLoop *eventloop, int runnumber)
{

	// Read in Magnetic field map
	DApplication* dapp = dynamic_cast<DApplication*>(eventloop->GetJApplication());
	Bfield = dapp->GetBfield();

	return NOERROR;
}

//------------------------------------------------------------------
// evnt 
//------------------------------------------------------------------
jerror_t MyProcessor::evnt(JEventLoop *eventLoop, int eventnumber)
{
	if(!eventLoop)return NOERROR;
	loop = eventLoop;
	last_jevent.FreeEvent();
	last_jevent = loop->GetJEvent();
	
	cout<<"----------- New Event "<<eventnumber<<" -------------"<<endl;
	hdvmf->SetEvent(eventnumber);
	hdvmf->DoRedraw();	

	return NOERROR;
}

//------------------------------------------------------------------
// FillGraphics 
//------------------------------------------------------------------
void MyProcessor::FillGraphics(void)
{
	/// Create "graphics" objects for this event given the current GUI settings.
	///
	/// This method will create DGraphicSet objects that represent tracks, hits,
	/// and showers for the event. It creates objects for both hits and
	/// reconstructed entities. The "graphics" objects created here are
	/// really just collections of 3D space points with flags indicating
	/// whether they should be drawn as markers or lines and with what
	/// color and size. The actual graphics objects are created for the
	/// various views of the detector in hdv_mainframe.
		
	graphics.clear();
	
	if(!loop)return;
	
	// FCAL hits
	if(hdvmf->GetCheckButton("fcal")){
		vector<const DFCALHit*> fcalhits;
		loop->Get(fcalhits);
		
		for(unsigned int i=0; i<fcalhits.size(); i++){
			const DFCALHit *hit = fcalhits[i];
			TPolyLine *poly = hdvmf->GetFCALPolyLine(hit->x, hit->y);
			if(!poly)continue;
			
			double a = hit->E/0.25;
			double f = sqrt(a>1.0 ? 1.0:a<0.0 ? 0.0:a);
			double grey = 0.8;
			double s = 1.0 - f;

			float r = s*grey;
			float g = s*grey;
			float b = f*(1.0-grey) + grey;

			poly->SetFillColor(TColor::GetColor(r,g,b));
		}
	}	
	
	// CDC hits
	if(hdvmf->GetCheckButton("cdc")){
		vector<const DCDCTrackHit*> cdctrackhits;
		loop->Get(cdctrackhits);
		
		for(unsigned int i=0; i<cdctrackhits.size(); i++){
			const DCDCWire *wire = cdctrackhits[i]->wire;
			
			DGraphicSet gset(kCyan, kLine, 1.0);
			gset.points.push_back(wire->origin-(wire->L/2.0)*wire->udir);
			gset.points.push_back(wire->origin+(wire->L/2.0)*wire->udir);
			graphics.push_back(gset);
		}
	}	

	// FDC wire
	if(hdvmf->GetCheckButton("fdcwire")){
		vector<const DFDCPseudo*> fdcpseudos;
		loop->Get(fdcpseudos);
		
		for(unsigned int i=0; i<fdcpseudos.size(); i++){
			const DFDCWire *wire = fdcpseudos[i]->wire;
			
			// Wire
			DGraphicSet gset(kCyan, kLine, 1.0);
			gset.points.push_back(wire->origin-(wire->L/2.0)*wire->udir);
			gset.points.push_back(wire->origin+(wire->L/2.0)*wire->udir);
			graphics.push_back(gset);
		}
	}
	
	// FDC psuedo hits
	if(hdvmf->GetCheckButton("fdcpseudo")){
		vector<const DFDCPseudo*> fdcpseudos;
		loop->Get(fdcpseudos);
		DGraphicSet gsetp(38, kMarker, 0.5);
		
		for(unsigned int i=0; i<fdcpseudos.size(); i++){
			const DFDCWire *wire = fdcpseudos[i]->wire;
			
			// Pseudo point
			TVector3 pos(fdcpseudos[i]->x, fdcpseudos[i]->y, wire->origin.Z());
			gsetp.points.push_back(pos);
		}
		graphics.push_back(gsetp);
	}

	// DMCThrown
	if(hdvmf->GetCheckButton("thrown")){
		vector<const DMCThrown*> mcthrown;
		loop->Get(mcthrown);
		for(unsigned int i=0; i<mcthrown.size(); i++){
			AddKinematicDataTrack(mcthrown[i], kGreen, 2.0);
		}
	}
	
	// CDC Truth points
	if(hdvmf->GetCheckButton("cdctruth")){	
		vector<const DMCTrackHit*> mctrackhits;
		loop->Get(mctrackhits);
		DGraphicSet gset(12, kMarker, 0.5);
		for(unsigned int i=0; i<mctrackhits.size(); i++){
			const DMCTrackHit *hit = mctrackhits[i];
			if(hit->system != SYS_CDC)continue;

			TVector3 pos(hit->r*cos(hit->phi), hit->r*sin(hit->phi), hit->z);
			gset.points.push_back(pos);
		}
		graphics.push_back(gset);
	}
	
	// FDC Truth points
	if(hdvmf->GetCheckButton("fdctruth")){	
		vector<const DMCTrackHit*> mctrackhits;
		loop->Get(mctrackhits);
		DGraphicSet gset(12, kMarker, 0.5);
		for(unsigned int i=0; i<mctrackhits.size(); i++){
			const DMCTrackHit *hit = mctrackhits[i];
			if(hit->system != SYS_FDC)continue;

			TVector3 pos(hit->r*cos(hit->phi), hit->r*sin(hit->phi), hit->z);
			gset.points.push_back(pos);
		}
		graphics.push_back(gset);
	}

	// TOF Truth points
	if(hdvmf->GetCheckButton("toftruth")){	
		vector<const DMCTrackHit*> mctrackhits;
		loop->Get(mctrackhits);
		DGraphicSet gset(kBlack, kMarker, 0.5);
		for(unsigned int i=0; i<mctrackhits.size(); i++){
			const DMCTrackHit *hit = mctrackhits[i];
			if(hit->system != SYS_TOF)continue;

			TVector3 pos(hit->r*cos(hit->phi), hit->r*sin(hit->phi), hit->z);
			gset.points.push_back(pos);
		}
		graphics.push_back(gset);
	}

	// BCAL Truth points
	if(hdvmf->GetCheckButton("bcaltruth")){	
		vector<const DMCTrackHit*> mctrackhits;
		loop->Get(mctrackhits);
		DGraphicSet gset(kBlack, kMarker, 1.0);
		for(unsigned int i=0; i<mctrackhits.size(); i++){
			const DMCTrackHit *hit = mctrackhits[i];
			if(hit->system != SYS_BCAL)continue;

			TVector3 pos(hit->r*cos(hit->phi), hit->r*sin(hit->phi), hit->z);
			gset.points.push_back(pos);
		}
		graphics.push_back(gset);
	}

	// FCAL Truth points
	if(hdvmf->GetCheckButton("fcaltruth")){	
		vector<const DMCTrackHit*> mctrackhits;
		loop->Get(mctrackhits);
		DGraphicSet gset(kBlack, kMarker, 1.0);
		for(unsigned int i=0; i<mctrackhits.size(); i++){
			const DMCTrackHit *hit = mctrackhits[i];
			if(hit->system != SYS_FCAL)continue;

			TVector3 pos(hit->r*cos(hit->phi), hit->r*sin(hit->phi), hit->z);
			gset.points.push_back(pos);
		}
		graphics.push_back(gset);
	}

	// DMCTrajectoryPoints
	if(hdvmf->GetCheckButton("trajectories")){
		vector<const DMCTrajectoryPoint*> mctrajectorypoints;
		loop->Get(mctrajectorypoints);
		
		DGraphicSet gset(kBlack, kLine, 1.5);
		for(unsigned int i=0; i<mctrajectorypoints.size(); i++){
			const DMCTrajectoryPoint *pt = mctrajectorypoints[i];
			
			TVector3 v(pt->x, pt->y, pt->z);
			gset.points.push_back(v);
		}
		graphics.push_back(gset);
	}

	// DTrackCandidate
	if(hdvmf->GetCheckButton("candidates")){
		vector<const DTrackCandidate*> trackcandidates;
		loop->Get(trackcandidates, hdvmf->GetFactoryTag("DTrackCandidate"));
		for(unsigned int i=0; i<trackcandidates.size(); i++){
			AddKinematicDataTrack(trackcandidates[i], kMagenta, 1.75);
		}
	}

	// DTrack
	if(hdvmf->GetCheckButton("tracks")){
		vector<const DTrack*> tracks;
		loop->Get(tracks, hdvmf->GetFactoryTag("DTrack"));
		for(unsigned int i=0; i<tracks.size(); i++){
			AddKinematicDataTrack(tracks[i], kRed, 1.00);
		}
	}

}

//------------------------------------------------------------------
// AddKinematicDataTrack 
//------------------------------------------------------------------
void MyProcessor::AddKinematicDataTrack(const DKinematicData* kd, int color, double size)
{
	// Create a reference trajectory with the given kinematic data and swim
	// it through the detector.
	DReferenceTrajectory rt(Bfield);
	rt.Swim(kd->position(), kd->momentum(), kd->charge());

	// Create a new graphics set and fill it with all of the trajectory points
	DGraphicSet gset(color, kLine, size);
	DReferenceTrajectory::swim_step_t *step = rt.swim_steps;
	for(int i=0; i<rt.Nswim_steps; i++, step++){
		gset.points.push_back(step->origin);
	}
	
	// Push the graphics set onto the stack
	graphics.push_back(gset);
}

