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

#include <particleType.h>
#include "hdview2.h"
#include "hdv_mainframe.h"
#include "MyProcessor.h"
#include "TRACKING/DTrackHit.h"
#include "TRACKING/DQuickFit.h"
#include "TRACKING/DMagneticFieldStepper.h"
#include "TRACKING/DTrackCandidate_factory.h"
#include "TRACKING/DMCTrackHit.h"
#include "TRACKING/DMCThrown.h"
#include "TRACKING/DTrackWireBased.h"
#include "TRACKING/DTrackTimeBased.h"
#include "PID/DParticle.h"
#include "PID/DChargedTrack.h"
#include "TRACKING/DReferenceTrajectory.h"
#include "JANA/JGeometry.h"
#include "TRACKING/DMCTrajectoryPoint.h"
#include "FCAL/DFCALHit.h"
#include "FDC/DFDCGeometry.h"
#include "CDC/DCDCTrackHit.h"
#include "FDC/DFDCPseudo.h"
#include "FDC/DFDCIntersection.h"
#include "HDGEOMETRY/DGeometry.h"
#include "FCAL/DFCALGeometry.h"
#include "FCAL/DFCALHit.h"
#include "PID/DPhoton.h"
#include "PID/DTwoGammaFit.h"
#include "BCAL/DBCALHit.h"
#include "DVector2.h"

extern hdv_mainframe *hdvmf;

// These are declared in hdv_mainframe.cc, but as static so we need to do it here as well (yechh!)
static float FCAL_Zmin = 622.8;
static float FCAL_Rmin = 6.0;
static float FCAL_Rmax = 212.0/2.0;
static float BCAL_Rmin=65.0;
static float BCAL_Zlen = 390.0;
static float BCAL_Zmin = 212.0 - BCAL_Zlen/2.0;


static vector<vector <DFDCWire *> >fdcwires;


bool DMCTrajectoryPoint_track_cmp(const DMCTrajectoryPoint *a,const DMCTrajectoryPoint *b){
	// sort by track number and then by particle type, then by z-coordinate
	// (yes, I saw the same track number have different particle types!)
	if(a->track != b->track)return a->track < b->track;
	if(a->part  != b->part )return a->part < b->part;
	return a->z < b->z;
}


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
		hdvmf->SetCandidateFactories(facnames);
		hdvmf->SetWireBasedTrackFactories(facnames);
		hdvmf->SetTimeBasedTrackFactories(facnames);
		hdvmf->SetReconstructedFactories(facnames);
		hdvmf->SetChargedTrackFactories(facnames);
		fulllistmf = hdvmf->GetFullListFrame();
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
	const DGeometry *dgeom  = dapp->GetDGeometry(runnumber);
	dgeom->GetFDCWires(fdcwires);

	RootGeom = dapp->GetRootGeom();
	geom = dapp->GetDGeometry(runnumber);
	
	MATERIAL_MAP_MODEL="DGeometry";
	gPARMS->SetDefaultParameter("TRKFIT:MATERIAL_MAP_MODEL",			MATERIAL_MAP_MODEL);

	eventloop->GetCalib("PID/photon_track_matching", photon_track_matching);
	DELTA_R_FCAL = photon_track_matching["DELTA_R_FCAL"];

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
	
	string source = "<no source>";
	if(last_jevent.GetJEventSource())source = last_jevent.GetJEventSource()->GetSourceName();
	
	cout<<"----------- New Event "<<eventnumber<<" -------------"<<endl;
	hdvmf->SetEvent(eventnumber);
	hdvmf->SetSource(source.c_str());
	hdvmf->DoMyRedraw();	

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
	graphics_xyA.clear(); // The objects placed in these will be deleted by hdv_mainframe
	graphics_xyB.clear(); // The objects placed in these will be deleted by hdv_mainframe
	graphics_xz.clear();  // The objects placed in these will be deleted by hdv_mainframe
	graphics_yz.clear();  // The objects placed in these will be deleted by hdv_mainframe
	
	if(!loop)return;
	
	// BCAL hits
	if(hdvmf->GetCheckButton("bcal")){
		vector<const DBCALHit*> bcalhits;
		loop->Get(bcalhits);
		
		for(unsigned int i=0; i<bcalhits.size(); i++){
			const DBCALHit *hit = bcalhits[i];
			TPolyLine *poly = hdvmf->GetBCALPolyLine(hit->module, hit->layer, hit->sector);
			if(!poly)continue;

			double a = hit->E/0.02;
			double f = sqrt(a>1.0 ? 1.0:a<0.0 ? 0.0:a);
			double grey = 0.8;
			double s = 1.0 - f;

			float r = s*grey;
			float g = s*grey;
			float b = f*(1.0-grey) + grey;

			poly->SetFillColor(TColor::GetColor(r,g,b));
			//poly->SetLineColor(TColor::GetColor(r,g,b));
			//poly->SetLineWidth(3.0);
			poly->SetFillStyle(3001);
		}
	}	

	// FCAL hits
	if(hdvmf->GetCheckButton("fcal")){
		vector<const DFCALHit*> fcalhits;
		loop->Get(fcalhits);
		
		for(unsigned int i=0; i<fcalhits.size(); i++){
			const DFCALHit *hit = fcalhits[i];
			TPolyLine *poly = hdvmf->GetFCALPolyLine(hit->x, hit->y);
			if(!poly)continue;

#if 0			
			double a = hit->E/0.005;
			double f = sqrt(a>1.0 ? 1.0:a<0.0 ? 0.0:a);
			double grey = 0.8;
			double s = 1.0 - f;

			float r = s*grey;
			float g = s*grey;
			float b = f*(1.0-grey) + grey;
#endif
			double s = log10(hit->E/0.005)/log10(1.0/0.005); // s=1 for 1GeV energy deposit
			if(s<0.0) s=0.0;
			float r = s;
			float g = s;
			float b = s;
			

			poly->SetFillColor(TColor::GetColor(r,g,b));
		}
	}	
	
	// CDC hits
	if(hdvmf->GetCheckButton("cdc")){
		vector<const DCDCTrackHit*> cdctrackhits;
		loop->Get(cdctrackhits);
		
		for(unsigned int i=0; i<cdctrackhits.size(); i++){
			const DCDCWire *wire = cdctrackhits[i]->wire;
			
			int color = (cdctrackhits[i]->tdrift>-20 && cdctrackhits[i]->tdrift<400) ? kCyan:kYellow;
			DGraphicSet gset(color, kLine, 1.0);
			DVector3 dpoint=wire->origin-(wire->L/2.0)*wire->udir;
			TVector3 tpoint(dpoint.X(),dpoint.Y(),dpoint.Z());
			gset.points.push_back(tpoint);
			dpoint=wire->origin+(wire->L/2.0)*wire->udir;
			tpoint.SetXYZ(dpoint.X(),dpoint.Y(),dpoint.Z());
			gset.points.push_back(tpoint);
			graphics.push_back(gset);
			
			// Rings for drift times.
			// NOTE: These are not perfect since they still have TOF in them
			if(hdvmf->GetCheckButton("cdcdrift") && wire->stereo==0.0){
				double x = wire->origin.X();
				double y = wire->origin.Y();
				double dist1 = cdctrackhits[i]->dist;
				TEllipse *e = new TEllipse(x, y, dist1, dist1);
				e->SetLineColor(38);
				e->SetFillStyle(0);
				graphics_xyA.push_back(e);

				double dist2 = dist1 - 4.0*55.0E-4; // what if our TOF was 4ns?
				e = new TEllipse(x, y, dist2, dist2);
				e->SetLineColor(38);
				e->SetFillStyle(0);
				graphics_xyA.push_back(e);
			}
		}
	}	

	// FDC wire
	if(hdvmf->GetCheckButton("fdcwire")){
		vector<const DFDCHit*> fdchits;
		loop->Get(fdchits);

		for(unsigned int i=0; i<fdchits.size(); i++){
			const DFDCHit *fdchit = fdchits[i];
			if(fdchit->type!=0)continue;
			const DFDCWire *wire =fdcwires[fdchit->gLayer-1][fdchit->element-1];
			if(!wire){
				_DBG_<<"Couldn't find wire for gLayer="<<fdchit->gLayer<<" and element="<<fdchit->element<<endl;
				continue;
			}
			
			// Wire
			int color = (fdchit->t>-50 && fdchit->t<400) ? kCyan:kYellow;
			DGraphicSet gset(color, kLine, 1.0);
			DVector3 dpoint=wire->origin-(wire->L/2.0)*wire->udir;
			TVector3 tpoint(dpoint.X(),dpoint.Y(),dpoint.Z());
			gset.points.push_back(tpoint);
			dpoint=wire->origin+(wire->L/2.0)*wire->udir;
			tpoint.SetXYZ(dpoint.X(),dpoint.Y(),dpoint.Z());
			gset.points.push_back(tpoint);
			graphics.push_back(gset);
		}		
	}

	// FDC intersection hits
	if(hdvmf->GetCheckButton("fdcintersection")){
		vector<const DFDCIntersection*> fdcints;
		loop->Get(fdcints);
		DGraphicSet gsetp(46, kMarker, 0.5);
		
		for(unsigned int i=0; i<fdcints.size(); i++){
		  TVector3 tpos(fdcints[i]->pos.X(),fdcints[i]->pos.Y(),
				fdcints[i]->pos.Z());
		  gsetp.points.push_back(tpos);
		}
		graphics.push_back(gsetp);
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
			int color=14;
			double size=1.5;
			if(mcthrown[i]->charge()==0.0) color = kGreen;
			if(mcthrown[i]->charge() >0.0) color = kBlue;
			if(mcthrown[i]->charge() <0.0) color = kRed;
			switch(mcthrown[i]->type){
				case Gamma:
				case Positron:
				case Electron:
					size = 1.0;
					break;
				case Pi0:
				case PiPlus:
				case PiMinus:
					size = 2.0;
					break;
				case Neutron:
				case Proton:
				case AntiProton:
					size = 3.0;
					break;
			}
			AddKinematicDataTrack(mcthrown[i], color, size);
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
		vector<const DFCALGeometry*> fcalgeometries;
		vector<const DFCALHit*> mcfcalhits;
		loop->Get(fcalgeometries);
		loop->Get(mcfcalhits);
		if(fcalgeometries.size()>0){
			const DFCALGeometry *fgeom = fcalgeometries[0];

			DGraphicSet gset(kBlack, kMarker, 0.25);
			for(unsigned int i=0; i<mcfcalhits.size(); i++){
				const DFCALHit *hit = mcfcalhits[i];

				DVector2 pos_face = fgeom->positionOnFace(hit->row, hit->column);
				TVector3 pos(pos_face.X(), pos_face.Y(), FCAL_Zmin);
				gset.points.push_back(pos);
			}
			graphics.push_back(gset);
		}
	}

	// BCAL reconstructed photons
	if(hdvmf->GetCheckButton("recon_photons_bcal")){
		vector<const DPhoton*> photons;
		loop->Get(photons);
		DGraphicSet gset(kYellow+2, kMarker, 1.25);
		gset.marker_style=21;
		for(unsigned int i=0; i<photons.size(); i++){
			const DPhoton *photon = photons[i];
			if(photon->getTag()!=DPhoton::kBcal)continue;
			
			TVector3 pos(photon->getPositionCal().X(), photon->getPositionCal().Y(), photon->getPositionCal().Z());
			gset.points.push_back(pos);
		}
		graphics.push_back(gset);
	}

	// FCAL reconstructed photons
	if(hdvmf->GetCheckButton("recon_photons_fcal")){
		vector<const DPhoton*> photons;
		loop->Get(photons);
		DGraphicSet gset(kOrange, kMarker, 1.25);
		for(unsigned int i=0; i<photons.size(); i++){
			const DPhoton *photon = photons[i];
			if(photon->getTag()!=DPhoton::kFcal)continue;
			
			TVector3 pos(photon->getPositionCal().X(), photon->getPositionCal().Y(), photon->getPositionCal().Z());
			gset.points.push_back(pos);
			
			double dist2 = 2.0 + 2.0*photon->energy();
			TEllipse *e = new TEllipse(photon->getPositionCal().X(), photon->getPositionCal().Y(), dist2, dist2);
			e->SetLineColor(kOrange);
			e->SetFillStyle(0);
			e->SetLineWidth(2);
			graphics_xyB.push_back(e);
		}
		graphics.push_back(gset);
	}

	// Reconstructed photons matched with tracks
	if(hdvmf->GetCheckButton("recon_photons_track_match")){
		vector<const DPhoton*> photons;
		loop->Get(photons);
		for(unsigned int i=0; i<photons.size(); i++){
			const DPhoton *photon = photons[i];
			if(photon->getTag()!=DPhoton::kCharge)continue;
			
			// Decide if this hit BCAL of FCAL based on z of position on calorimeter
			bool is_bcal = photon->getPositionCal().Z()<(FCAL_Zmin-10.0);
			
			// Draw on all frames except FCAL frame
			DGraphicSet gset(kRed, kMarker, 1.25);
			gset.marker_style = is_bcal ? 22:3;
			TVector3 tpos(photon->getPositionCal().X(),
				      photon->getPositionCal().Y(),
				      photon->getPositionCal().Z());
			gset.points.push_back(tpos);
			graphics.push_back(gset);
			
			// For BCAL hits, don't draw them on FCAL pane
			if(is_bcal)continue;

			// Draw on FCAL pane
			double dist2 = 2.0 + 2.0*photon->energy();
			TEllipse *e = new TEllipse(photon->getPositionCal().X(), photon->getPositionCal().Y(), dist2, dist2);
			e->SetLineColor(gset.color);
			e->SetFillStyle(0);
			e->SetLineWidth(1);
			TMarker *m = new TMarker(photon->getPositionCal().X(), photon->getPositionCal().Y(), gset.marker_style);
			m->SetMarkerColor(gset.color);
			m->SetMarkerSize(1.75);
			graphics_xyB.push_back(e);
			graphics_xyB.push_back(m);
		}
	}
	
	// FCAL and BCAL thrown photon projections
	if(hdvmf->GetCheckButton("thrown_photons_fcal") || hdvmf->GetCheckButton("thrown_photons_bcal")){
		vector<const DMCThrown*> throwns;
		loop->Get(throwns);
		DGraphicSet gset(kSpring, kMarker, 1.25);
		for(unsigned int i=0; i<throwns.size(); i++){
			const DMCThrown *thrown = throwns[i];
			if(thrown->charge()!=0.0)continue;
			
			// This may seem a little funny, but it saves swimming the reference trajectory
			// multiple times. The GetIntersectionWithCalorimeter() method will find the
			// intersection point of the photon with whichever calorimeter it actually hits
			// or if it doesn't hit either of them. Then, we decide to draw the marker based
			// on whether the flag is set to draw the calorimeter it hit or not.
			DVector3 pos;
			DPhoton::PhotonTag who;
			GetIntersectionWithCalorimeter(thrown, pos, who);
			
			if(who!=DPhoton::kFcal && who!=DPhoton::kBcal)continue;
			if(who==DPhoton::kFcal && !hdvmf->GetCheckButton("thrown_photons_fcal"))continue;
			if(who==DPhoton::kBcal && !hdvmf->GetCheckButton("thrown_photons_bcal"))continue;
			TVector3 tpos(pos.X(),pos.Y(),pos.Z());
			gset.points.push_back(tpos);
			
			// Only draw on FCAL pane if photon hits FCAL
			if(who==DPhoton::kBcal)continue;
			
			double dist2 = 2.0 + 2.0*thrown->energy();
			TEllipse *e = new TEllipse(pos.X(), pos.Y(), dist2, dist2);
			e->SetLineColor(kSpring);
			e->SetFillStyle(0);
			e->SetLineWidth(4);
			graphics_xyB.push_back(e);
		}
		graphics.push_back(gset);
	}
	
	// FCAL and BCAL thrown charged particle projections
	if(hdvmf->GetCheckButton("thrown_charged_fcal") || hdvmf->GetCheckButton("thrown_charged_bcal")){
		vector<const DMCThrown*> throwns;
		loop->Get(throwns);
		
		for(unsigned int i=0; i<throwns.size(); i++){
			const DMCThrown *thrown = throwns[i];
			if(thrown->charge()==0.0)continue;
			
			// This may seem a little funny, but it saves swimming the reference trajectory
			// multiple times. The GetIntersectionWithCalorimeter() method will find the
			// intersection point of the photon with whichever calorimeter it actually hits
			// or if it doesn't hit either of them. Then, we decide to draw the marker based
			// on whether the flag is set to draw the calorimeter it hit or not.
			DVector3 pos;
			DPhoton::PhotonTag who;
			GetIntersectionWithCalorimeter(thrown, pos, who);
			
			if(who!=DPhoton::kFcal && who!=DPhoton::kBcal)continue;
			if(who==DPhoton::kFcal && !hdvmf->GetCheckButton("thrown_charged_fcal"))continue;
			if(who==DPhoton::kBcal && !hdvmf->GetCheckButton("thrown_charged_bcal"))continue;
			
			DGraphicSet gset(thrown->charge()>0.0 ? kBlue:kRed, kMarker, 1.25);
			TVector3 tpos(pos.X(),pos.Y(),pos.Z());
			gset.points.push_back(tpos);
			graphics.push_back(gset);
			
			//double dist2 = 6.0 + 2.0*thrown->momentum().Mag();
			double dist2 = DELTA_R_FCAL;
			TEllipse *e = new TEllipse(pos.X(), pos.Y(), dist2, dist2);
			e->SetLineColor(thrown->charge()>0.0 ? kBlue:kRed);
			e->SetFillStyle(0);
			e->SetLineWidth(4);
			graphics_xyB.push_back(e);
		}
	}
	
	// FCAL and BCAL reconstructed charged particle projections
	if(hdvmf->GetCheckButton("recon_charged_bcal") || hdvmf->GetCheckButton("recon_charged_fcal")){
	
		// Here we used the full time-based track list, even though it includes multiple
		// hypotheses for each candidate. This is because currently, the photon/track
		// matching code in PID/DPhoton_factory.cc uses the DTrackTimeBased objects and
		// the current purpose of drawing these is to see matching of reconstructed
		// charged tracks with calorimeter clusters.
		vector<const DTrackTimeBased*> tracks;
		loop->Get(tracks, hdvmf->GetFactoryTag("DTrackTimeBased"));
		
		for(unsigned int i=0; i<tracks.size(); i++){
			const DTrackTimeBased *track = tracks[i];
			
			// See notes for above sections
			DVector3 pos;
			DPhoton::PhotonTag who;
			GetIntersectionWithCalorimeter(track, pos, who);
			
			if(who!=DPhoton::kFcal && who!=DPhoton::kBcal)continue;
			if(who==DPhoton::kFcal && !hdvmf->GetCheckButton("recon_charged_fcal"))continue;
			if(who==DPhoton::kBcal && !hdvmf->GetCheckButton("recon_charged_bcal"))continue;
			
			DGraphicSet gset(track->charge()>0.0 ? kBlue:kRed, kMarker, 1.25);
			TVector3 tpos(pos.X(),pos.Y(),pos.Z());
			gset.points.push_back(tpos);
			graphics.push_back(gset);
			
			if(who==DPhoton::kBcal)continue; // Don't draw tracks hitting BCAL on FCAL pane
			
			//double dist2 = 6.0 + 2.0*track->momentum().Mag();
			double dist2 = DELTA_R_FCAL;
			TEllipse *e = new TEllipse(pos.X(), pos.Y(), dist2, dist2);
			e->SetLineColor(track->charge()>0.0 ? kBlue:kRed);
			e->SetFillStyle(0);
			e->SetLineWidth(4);
			graphics_xyB.push_back(e);
		}
	}

	// DMCTrajectoryPoints
	if(hdvmf->GetCheckButton("trajectories")){
		vector<const DMCTrajectoryPoint*> mctrajectorypoints;
		loop->Get(mctrajectorypoints);
		//sort(mctrajectorypoints.begin(), mctrajectorypoints.end(), DMCTrajectoryPoint_track_cmp);
		
		int colors[] = {kBlack, kMagenta, kGreen, 13, 14, 39};
		int Ncolors = 1;
		int Ntraj=0;
		DGraphicSet gset(colors[Ntraj%Ncolors], kLine, 3.0);
		//gset.marker_style = 7;
		TVector3 last_point;
		int last_track=-1;
		int last_part=-1;
		for(unsigned int i=0; i<mctrajectorypoints.size(); i++){
			const DMCTrajectoryPoint *pt = mctrajectorypoints[i];
			
			switch(pt->part){
				case	Gamma:
					if(!hdvmf->GetCheckButton("trajectories_photon"))continue;
					break;
				case	Electron:
					if(!hdvmf->GetCheckButton("trajectories_electron"))continue;
					break;
				case	Positron:
					if(!hdvmf->GetCheckButton("trajectories_positron"))continue;
					break;
				case	Proton:
					if(!hdvmf->GetCheckButton("trajectories_proton"))continue;
					break;
				case	Neutron:
					if(!hdvmf->GetCheckButton("trajectories_neutron"))continue;
					break;
				case	PiPlus:
					if(!hdvmf->GetCheckButton("trajectories_piplus"))continue;
					break;
				case	PiMinus:
					if(!hdvmf->GetCheckButton("trajectories_piminus"))continue;
					break;
				default:
					if(!hdvmf->GetCheckButton("trajectories_other"))continue;
					break;
			}
						
			TVector3 v(pt->x, pt->y, pt->z);

			if(i>0){
				//if((v-last_point).Mag() > 10.0){
				if(pt->track!=last_track || pt->part!=last_part){
					graphics.push_back(gset);
					gset.points.clear();
					gset.color = colors[(++Ntraj)%Ncolors];
				}
			}
			
			gset.points.push_back(v);
			last_point = v;
			last_track = pt->track;
			last_part = pt->part;
		}
		graphics.push_back(gset);
	}

	// DTrackCandidate
	if(hdvmf->GetCheckButton("candidates")){
		vector<const DTrackCandidate*> trackcandidates;
		loop->Get(trackcandidates, hdvmf->GetFactoryTag("DTrackCandidate"));
		for(unsigned int i=0; i<trackcandidates.size(); i++){
			int color=30;
			double size=2.0;
			//if(trackcandidates[i]->charge() >0.0) color += 100; // lighter shade
			//if(trackcandidates[i]->charge() <0.0) color += 150; // darker shade
			AddKinematicDataTrack(trackcandidates[i], color, size);
		}
	}

	// DTrackWireBased
	if(hdvmf->GetCheckButton("wiretracks")){
		vector<const DTrackWireBased*> wiretracks;
		loop->Get(wiretracks, hdvmf->GetFactoryTag("DTrackWireBased"));
		for(unsigned int i=0; i<wiretracks.size(); i++){
			AddKinematicDataTrack(wiretracks[i], (wiretracks[i]->charge()>0.0 ? kBlue:kRed)+2, 1.25);
		}
	}

	// DTrackTimeBased
	if(hdvmf->GetCheckButton("timetracks")){
		vector<const DTrackTimeBased*> timetracks;
		loop->Get(timetracks, hdvmf->GetFactoryTag("DTrackTimeBased"));
		for(unsigned int i=0; i<timetracks.size(); i++){
			AddKinematicDataTrack(timetracks[i], (timetracks[i]->charge()>0.0 ? kBlue:kRed)+0, 1.00);
		}
	}

	// DChargedTrack
	if(hdvmf->GetCheckButton("chargedtracks")){
	  vector<const DChargedTrack*> chargedtracks;
		loop->Get(chargedtracks, hdvmf->GetFactoryTag("DChargedTrack"));
		for(unsigned int i=0; i<chargedtracks.size(); i++){
		  int color=kViolet-3;
		  double size=1.25;
		  if (chargedtracks[i]->hypotheses[0]->charge()>0) color=kMagenta;
		
		  if (chargedtracks[i]->hypotheses[0]->mass()>0.9) size=2.5;
		  AddKinematicDataTrack(chargedtracks[i]->hypotheses[0],color,size);
		}
	}

	// DPhoton
	if(hdvmf->GetCheckButton("photon")){
		vector<const DPhoton*> photons;
		loop->Get(photons, hdvmf->GetFactoryTag("DPhoton"));
		for(unsigned int i=0; i<photons.size(); i++){
			int color = kBlack;
			if(photons[i]->getTag()==DPhoton::kFcal)color = kOrange;
			if(photons[i]->getTag()==DPhoton::kBcal)color = kYellow+2;
			if(photons[i]->getTag()==DPhoton::kCharge)color = kRed;
			AddKinematicDataTrack(photons[i], color, 1.00);
		}
	}
}

//------------------------------------------------------------------
// UpdateTrackLabels 
//------------------------------------------------------------------
void MyProcessor::UpdateTrackLabels(void)
{
	// Get the label pointers
	string name, tag;
	map<string, vector<TGLabel*> > &thrownlabs = hdvmf->GetThrownLabels();
	map<string, vector<TGLabel*> > &reconlabs = hdvmf->GetReconstructedLabels();
	hdvmf->GetReconFactory(name, tag);
	
	// Get Thrown particles
	vector<const DMCThrown*> throwns;
	if(loop)loop->Get(throwns);

	// Get the track info as DKinematicData objects
	vector<const DKinematicData*> trks;
	if(name=="DChargedTrack"){
		vector<const DChargedTrack*> chargedtracks;
		if(loop)loop->Get(chargedtracks, tag.c_str());
		for(unsigned int i=0; i<chargedtracks.size(); i++){
		  trks.push_back(chargedtracks[i]->hypotheses[0]);
		}
	}	
	if(name=="DTrackTimeBased"){
		vector<const DTrackTimeBased*> timetracks;
		if(loop)loop->Get(timetracks, tag.c_str());
		for(unsigned int i=0; i<timetracks.size(); i++)trks.push_back(timetracks[i]);
	}
	if(name=="DTrackWireBased"){
		vector<const DTrackWireBased*> wiretracks;
		if(loop)loop->Get(wiretracks, tag.c_str());
		for(unsigned int i=0; i<wiretracks.size(); i++)trks.push_back(wiretracks[i]);
	}
	if(name=="DTrackCandidate"){
		vector<const DTrackCandidate*> candidates;
		if(loop)loop->Get(candidates, tag.c_str());
		for(unsigned int i=0; i<candidates.size(); i++)trks.push_back(candidates[i]);
	}
	if(name=="DPhoton"){
		vector<const DPhoton*> photons;
		if(loop)loop->Get(photons, tag.c_str());
		for(unsigned int i=0; i<photons.size(); i++)trks.push_back(photons[i]);
	}
	if(name=="DTwoGammaFit"){
		vector<const DTwoGammaFit*> twogammafits;
		if(loop)loop->Get(twogammafits, tag.c_str());
		for(unsigned int i=0; i<twogammafits.size(); i++)trks.push_back(twogammafits[i]);
	}
	
	// Clear all labels (i.e. draw ---- in them)
	map<string, vector<TGLabel*> >::iterator iter;
	for(iter=reconlabs.begin(); iter!=reconlabs.end(); iter++){
		vector<TGLabel*> &labs = iter->second;
		for(unsigned int i=1; i<labs.size(); i++){
			labs[i]->SetText("--------");
		}
	}
	for(iter=thrownlabs.begin(); iter!=thrownlabs.end(); iter++){
		vector<TGLabel*> &labs = iter->second;
		for(unsigned int i=1; i<labs.size(); i++){
			labs[i]->SetText("--------");
		}
	}

	// Loop over thrown particles and fill in labels
	int ii=0;
	for(unsigned int i=0; i<throwns.size(); i++){
		const DMCThrown *trk = throwns[i];
		if(trk->type==0)continue;
		int row = thrownlabs["trk"].size()-(ii++)-1;
		if(row<1)break;
		
		stringstream trkno, type, p, theta, phi, z;
		trkno<<setprecision(4)<<i+1;
		thrownlabs["trk"][row]->SetText(trkno.str().c_str());

		thrownlabs["type"][row]->SetText(ParticleType((Particle_t)trk->type));

		p<<setprecision(4)<<trk->momentum().Mag();
		thrownlabs["p"][row]->SetText(p.str().c_str());

		theta<<setprecision(4)<<trk->momentum().Theta()*TMath::RadToDeg();
		thrownlabs["theta"][row]->SetText(theta.str().c_str());

		double myphi = trk->momentum().Phi();
		if(myphi<0.0)myphi+=2.0*M_PI;
		phi<<setprecision(4)<<myphi;
		thrownlabs["phi"][row]->SetText(phi.str().c_str());

		z<<setprecision(4)<<trk->position().Z();
		thrownlabs["z"][row]->SetText(z.str().c_str());
	}

	// Loop over tracks and fill in labels
	for(unsigned int i=0; i<trks.size(); i++){
		const DKinematicData *trk = trks[i];
		int row = reconlabs["trk"].size()-i-1;
		if(row<1)break;
		
		stringstream trkno, type, p, theta, phi, z, chisq_per_dof, Ndof;
		stringstream fom;
		trkno<<setprecision(4)<<i+1;
		reconlabs["trk"][row]->SetText(trkno.str().c_str());
		
		double mass = trk->mass();
		if(fabs(mass-0.13957)<1.0E-4)type<<"pi";
		else if(fabs(mass-0.93827)<1.0E-4)type<<"proton";
		else if(fabs(mass-0.493677)<1.0E-4)type<<"K";
		else if(fabs(mass-0.000511)<1.0E-4)type<<"e";
		else if (fabs(mass)<1.e-4 && fabs(trk->charge())<1.e-4) type << "gamma";
		else type<<"q=";
		if (fabs(trk->charge())>1.e-4){
		  type<<(trk->charge()>0 ? "+":"-");
		}
		reconlabs["type"][row]->SetText(type.str().c_str());

		p<<setprecision(4)<<trk->momentum().Mag();
		reconlabs["p"][row]->SetText(p.str().c_str());

		theta<<setprecision(4)<<trk->momentum().Theta()*TMath::RadToDeg();
		reconlabs["theta"][row]->SetText(theta.str().c_str());

		double myphi = trk->momentum().Phi();
		if(myphi<0.0)myphi+=2.0*M_PI;
		phi<<setprecision(4)<<myphi;
		reconlabs["phi"][row]->SetText(phi.str().c_str());

		z<<setprecision(4)<<trk->position().Z();
		reconlabs["z"][row]->SetText(z.str().c_str());

		// Get chisq and Ndof for DTrackTimeBased or DTrackWireBased objects
		const DTrackTimeBased *timetrack=dynamic_cast<const DTrackTimeBased*>(trk);
		const DTrackWireBased *track=dynamic_cast<const DTrackWireBased*>(trk);	
		const DTwoGammaFit *twogammafit=dynamic_cast<const DTwoGammaFit*>(trk);	
	
		if(timetrack){
			chisq_per_dof<<setprecision(4)<<timetrack->chisq/timetrack->Ndof;
			Ndof<<timetrack->Ndof;
			fom << timetrack->FOM;
		}else if(track){
			chisq_per_dof<<setprecision(4)<<track->chisq/track->Ndof;
			Ndof<<track->Ndof;
			fom << "N/A";
		}else if(twogammafit){
			chisq_per_dof<<setprecision(4)<<twogammafit->getChi2();
			Ndof<<twogammafit->getNdf();
			fom << twogammafit->getProb();
		}else{
			chisq_per_dof<<"N/A";
			Ndof<<"N/A";
			fom << "N/A";
		}
		reconlabs["chisq/Ndof"][row]->SetText(chisq_per_dof.str().c_str());
		reconlabs["Ndof"][row]->SetText(Ndof.str().c_str());
		reconlabs["FOM"][row]->SetText(fom.str().c_str());
	}
	
	// Have the pop-up window with the full particle list update it's labels
	fulllistmf->UpdateTrackLabels(throwns, trks);
}

//------------------------------------------------------------------
// AddKinematicDataTrack 
//------------------------------------------------------------------
void MyProcessor::AddKinematicDataTrack(const DKinematicData* kd, int color, double size)
{
	// Create a reference trajectory with the given kinematic data and swim
	// it through the detector.
	DReferenceTrajectory rt(Bfield);

	if(MATERIAL_MAP_MODEL=="DRootGeom"){
		rt.SetDRootGeom(RootGeom);
		rt.SetDGeometry(NULL);
	}else if(MATERIAL_MAP_MODEL=="DGeometry"){
		rt.SetDRootGeom(NULL);
		rt.SetDGeometry(geom);
	}else if(MATERIAL_MAP_MODEL!="NONE"){
		_DBG_<<"WARNING: Invalid value for TRKFIT:MATERIAL_MAP_MODEL (=\""<<MATERIAL_MAP_MODEL<<"\")"<<endl;
	}

	rt.SetMass(kd->mass());
	rt.Swim(kd->position(), kd->momentum(), kd->charge());

	// Create a new graphics set and fill it with all of the trajectory points
	DGraphicSet gset(color, kLine, size);
	DReferenceTrajectory::swim_step_t *step = rt.swim_steps;
	for(int i=0; i<rt.Nswim_steps; i++, step++){
	  TVector3 tpoint(step->origin.X(),step->origin.Y(),step->origin.Z());
		gset.points.push_back(tpoint);
	}
	
	// Push the graphics set onto the stack
	graphics.push_back(gset);
}

//------------------------------------------------------------------
// GetIntersectionWithCalorimeter 
//------------------------------------------------------------------
void MyProcessor::GetIntersectionWithCalorimeter(const DKinematicData* kd, DVector3 &pos, DPhoton::PhotonTag &who)
{
	// Create a reference trajectory with the given kinematic data and swim
	// it through the detector.
	DReferenceTrajectory rt(Bfield);

	if(MATERIAL_MAP_MODEL=="DRootGeom"){
		rt.SetDRootGeom(RootGeom);
		rt.SetDGeometry(NULL);
	}else if(MATERIAL_MAP_MODEL=="DGeometry"){
		rt.SetDRootGeom(NULL);
		rt.SetDGeometry(geom);
	}else if(MATERIAL_MAP_MODEL!="NONE"){
		_DBG_<<"WARNING: Invalid value for TRKFIT:MATERIAL_MAP_MODEL (=\""<<MATERIAL_MAP_MODEL<<"\")"<<endl;
	}

	rt.SetMass(kd->mass());
	rt.Swim(kd->position(), kd->momentum(), kd->charge());

	// Find intersection with FCAL
	DVector3 pos_fcal;
	double s_fcal = 1.0E6;
	DVector3 origin(0.0, 0.0, FCAL_Zmin);
	DVector3 norm(0.0, 0.0, -1.0);
	rt.GetIntersectionWithPlane(origin, norm, pos_fcal, &s_fcal);
	if(pos_fcal.Perp()<FCAL_Rmin || pos_fcal.Perp()>FCAL_Rmax)s_fcal = 1.0E6;
	
	// Find intersection with BCAL
	DVector3 pos_bcal;
	double s_bcal = 1.0E6;
	rt.GetIntersectionWithRadius(BCAL_Rmin, pos_bcal, &s_bcal);
	if(pos_bcal.Z()<BCAL_Zmin || pos_bcal.Z()>(BCAL_Zmin+BCAL_Zlen))s_bcal = 1.0E6;

	if(s_fcal>1000.0 && s_bcal>1000.0){
		// neither calorimeter hit
		who = DPhoton::kDefaultTag;
		pos.SetXYZ(0.0,0.0,0.0);
	}else if(s_fcal<s_bcal){
		// FCAL hit
		who = DPhoton::kFcal;
		pos = pos_fcal;
	}else{
		// BCAL hit
		who = DPhoton::kBcal;
		pos = pos_bcal;
	}
}

//------------------------------------------------------------------
// GetFactoryNames 
//------------------------------------------------------------------
void MyProcessor::GetFactoryNames(vector<string> &facnames)
{
	vector<JEventLoop*> loops = app->GetJEventLoops();
	if(loops.size()>0){
		vector<string> facnames;
		loops[0]->GetFactoryNames(facnames);
	}
}

//------------------------------------------------------------------
// GetFactories 
//------------------------------------------------------------------
void MyProcessor::GetFactories(vector<JFactory_base*> &factories)
{
	vector<JEventLoop*> loops = app->GetJEventLoops();
	if(loops.size()>0){
		factories = loops[0]->GetFactories();
	}
}

//------------------------------------------------------------------
// GetNrows 
//------------------------------------------------------------------
unsigned int MyProcessor::GetNrows(const string &factory, string tag)
{
	vector<JEventLoop*> loops = app->GetJEventLoops();
	if(loops.size()>0){
		// Here is something a little tricky. The GetFactory() method of JEventLoop
		// gets the factory of the specified data name and tag, but without trying
		// to substitute a user-specified tag (a'la -PDEFTAG:XXX=YYY) as is done
		// on normal Get() method calls. Therefore, we have to check for the default
		// tags ourselves and substitute it "by hand".
		if(tag==""){
			map<string,string> default_tags = loops[0]->GetDefaultTags();
			map<string, string>::const_iterator iter = default_tags.find(factory);
			if(iter!=default_tags.end())tag = iter->second.c_str();
		}
		JFactory_base *fac = loops[0]->GetFactory(factory, tag.c_str());

		// Since calling GetNrows will cause the program to quit if there is
		// not a valid event, then first check that there is one before calling it
		if(loops[0]->GetJEvent().GetJEventSource() == NULL)return 0;
		
		return fac==NULL ? 0:(unsigned int)fac->GetNrows();
	}
	
	return 0;
}

//------------------------------------------------------------------
// GetDReferenceTrajectory 
//------------------------------------------------------------------
void MyProcessor::GetDReferenceTrajectory(string dataname, string tag, unsigned int index, DReferenceTrajectory* &rt, vector<const DCDCTrackHit*> &cdchits)
{
_DBG__;
	// initialize rt to NULL in case we don't find the one requested
	rt = NULL;
	cdchits.clear();

	// Get pointer to the JEventLoop so we can get at the data
	vector<JEventLoop*> loops = app->GetJEventLoops();
	if(loops.size()==0)return;
	JEventLoop* &loop = loops[0];
	
	// Variables to hold track parameters
	DVector3 pos, mom(0,0,0);
	double q=0.0;
	double mass;

	// Find the specified track	
	if(dataname=="DChargedTrack"){
		vector<const DChargedTrack*> chargedtracks;
		loop->Get(chargedtracks, tag.c_str());
		if(index>=chargedtracks.size())return;
		q = chargedtracks[index]->hypotheses[0]->charge();
		pos = chargedtracks[index]->hypotheses[0]->position();
		mom = chargedtracks[index]->hypotheses[0]->momentum();
		chargedtracks[index]->hypotheses[0]->Get(cdchits);
		mass = chargedtracks[index]->hypotheses[0]->mass();
	}

	if(dataname=="DTrackTimeBased"){
		vector<const DTrackTimeBased*> timetracks;
		loop->Get(timetracks, tag.c_str());
		if(index>=timetracks.size())return;
		q = timetracks[index]->charge();
		pos = timetracks[index]->position();
		mom = timetracks[index]->momentum();
		timetracks[index]->Get(cdchits);
		mass = timetracks[index]->mass();
	}

	if(dataname=="DTrackWireBased"){
		vector<const DTrackWireBased*> wiretracks;
		loop->Get(wiretracks, tag.c_str());
		if(index>=wiretracks.size())return;
		q = wiretracks[index]->charge();
		pos = wiretracks[index]->position();
		mom = wiretracks[index]->momentum();
		wiretracks[index]->Get(cdchits);
		mass = wiretracks[index]->mass();
	}

	if(dataname=="DTrackCandidate"){
		vector<const DTrackCandidate*> tracks;
		loop->Get(tracks, tag.c_str());
		if(index>=tracks.size())return;
		q = tracks[index]->charge();
		pos = tracks[index]->position();
		mom = tracks[index]->momentum();
		tracks[index]->Get(cdchits);
		mass = tracks[index]->mass();
	}

	if(dataname=="DMCThrown"){
		vector<const DMCThrown*> tracks;
		loop->Get(tracks, tag.c_str());
		if(index>=tracks.size())return;
		const DMCThrown *t = tracks[index];
		q = t->charge();
		pos = t->position();
		mom = t->momentum();
		tracks[index]->Get(cdchits);
		mass = tracks[index]->mass();
_DBG_<<"mass="<<mass<<endl;
	}

	// Make sure we found a charged particle we can track
	if(q==0.0 || mom.Mag()<0.01)return;
	
	// Create a new DReference trajectory object. The caller takes
	// ownership of this and so they are responsible for deleting it.
	rt = new DReferenceTrajectory(Bfield);
	if(MATERIAL_MAP_MODEL=="DRootGeom"){
		rt->SetDRootGeom(RootGeom);
		rt->SetDGeometry(NULL);
	}else if(MATERIAL_MAP_MODEL=="DGeometry"){
		rt->SetDRootGeom(NULL);
		rt->SetDGeometry(geom);
	}else if(MATERIAL_MAP_MODEL!="NONE"){
		_DBG_<<"WARNING: Invalid value for TRKFIT:MATERIAL_MAP_MODEL (=\""<<MATERIAL_MAP_MODEL<<"\")"<<endl;
	}
	rt->Swim(pos, mom, q);
}

//------------------------------------------------------------------
// GetAllWireHits
//------------------------------------------------------------------
void MyProcessor::GetAllWireHits(vector<pair<const DCoordinateSystem*,double> > &allhits)
{
	/// Argument is vector of pairs that contain a pointer to the
	/// DCoordinateSystem representing a wire and a double that
	/// represents the drift distance. To get info on the specific
	/// wire, one needs to attempt a dynamic_cast to both a DCDCWire
	/// and a DFDCWire and access the parameters of whichever one succeeds.
	
	// Get pointer to the JEventLoop so we can get at the data
	vector<JEventLoop*> loops = app->GetJEventLoops();
	if(loops.size()==0)return;
	JEventLoop* &loop = loops[0];

	// Get CDC wire hits
	vector<const DCDCTrackHit*> cdchits;
	loop->Get(cdchits);
	for(unsigned int i=0; i<cdchits.size(); i++){
		pair<const DCoordinateSystem*,double> hit;
		hit.first = cdchits[i]->wire;
		hit.second = cdchits[i]->dist;
		allhits.push_back(hit);
	}
	
	// Get FDC wire hits
	vector<const DFDCPseudo*> fdchits;
	loop->Get(fdchits);
	for(unsigned int i=0; i<fdchits.size(); i++){
		pair<const DCoordinateSystem*,double> hit;
		hit.first = fdchits[i]->wire;
		hit.second = 0.0055*fdchits[i]->time;
		allhits.push_back(hit);
	}
}



