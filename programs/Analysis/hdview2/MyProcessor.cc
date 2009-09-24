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
#include "TRACKING/DTrack.h"
#include "PID/DParticle.h"
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
#include "BCAL/DHDDMBCALHit.h"

extern hdv_mainframe *hdvmf;

// This is declared in hdv_mainframe.cc, but as static so we need to do it here as well (yechh!)
static float FCAL_Zmin = 622.8;


static vector<vector <DFDCWire *> >fdcwires;


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
		hdvmf->SetTrackFactories(facnames);
		hdvmf->SetParticleFactories(facnames);
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
	const DGeometry *dgeom  = dapp->GetDGeometry(runnumber);
	dgeom->GetFDCWires(fdcwires);

	RootGeom = dapp->GetRootGeom();
	geom = dapp->GetDGeometry(runnumber);
	
	MATERIAL_MAP_MODEL="DGeometry";
	gPARMS->SetDefaultParameter("TRKFIT:MATERIAL_MAP_MODEL",			MATERIAL_MAP_MODEL);

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
	hdvmf->SetSource(last_jevent.GetJEventSource()->GetSourceName());
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
		vector<const DHDDMBCALHit*> bcalhits;
		loop->Get(bcalhits);
		
		for(unsigned int i=0; i<bcalhits.size(); i++){
			const DHDDMBCALHit *hit = bcalhits[i];
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
			
			double a = hit->E/0.005;
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
			
			int color = (cdctrackhits[i]->tdrift>-20 && cdctrackhits[i]->tdrift<400) ? kCyan:kYellow;
			DGraphicSet gset(color, kLine, 1.0);
			gset.points.push_back(wire->origin-(wire->L/2.0)*wire->udir);
			gset.points.push_back(wire->origin+(wire->L/2.0)*wire->udir);
			graphics.push_back(gset);
			
			// Rings for drift times.
			// NOTE: These are not perfect since they still have TOF in them
			if(hdvmf->GetCheckButton("cdcdrift") && wire->stereo==0.0){
				double x = wire->origin.X();
				double y = wire->origin.Y();
				double dist1 = cdctrackhits[i]->dist;
				TEllipse *e = new TEllipse(x, y, dist1, dist1);
				e->SetLineColor(38);
				graphics_xyA.push_back(e);

				double dist2 = dist1 - 4.0*55.0E-4; // what if our TOF was 4ns?
				e = new TEllipse(x, y, dist2, dist2);
				e->SetLineColor(38);
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
			gset.points.push_back(wire->origin-(wire->L/2.0)*wire->udir);
			gset.points.push_back(wire->origin+(wire->L/2.0)*wire->udir);
			graphics.push_back(gset);
		}		
	}

	// FDC intersection hits
	if(hdvmf->GetCheckButton("fdcintersection")){
		vector<const DFDCIntersection*> fdcints;
		loop->Get(fdcints);
		DGraphicSet gsetp(46, kMarker, 0.5);
		
		for(unsigned int i=0; i<fdcints.size(); i++){
			gsetp.points.push_back(fdcints[i]->pos);
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

				TVector2 pos_face = fgeom->positionOnFace(hit->row, hit->column);
				TVector3 pos(pos_face.X(), pos_face.Y(), FCAL_Zmin);
				gset.points.push_back(pos);
			}
			graphics.push_back(gset);
		}
	}

	// DMCTrajectoryPoints
	if(hdvmf->GetCheckButton("trajectories")){
		vector<const DMCTrajectoryPoint*> mctrajectorypoints;
		loop->Get(mctrajectorypoints);
		
		int colors[] = {kBlack, kMagenta, kGreen, 13, 14, 39};
		int Ncolors = 1;
		int Ntraj=0;
		DGraphicSet gset(colors[Ntraj%Ncolors], kLine, 3.0);
		TVector3 last_point;
		for(unsigned int i=0; i<mctrajectorypoints.size(); i++){
			const DMCTrajectoryPoint *pt = mctrajectorypoints[i];
						
			TVector3 v(pt->x, pt->y, pt->z);

			if(i>0){
				if((v-last_point).Mag() > 10.0){
					graphics.push_back(gset);
					gset.points.clear();
					gset.color = colors[(++Ntraj)%Ncolors];
				}
			}
			
			gset.points.push_back(v);
			last_point = v;
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

	// DTrack
	if(hdvmf->GetCheckButton("tracks")){
		vector<const DTrack*> tracks;
		loop->Get(tracks, hdvmf->GetFactoryTag("DTrack"));
		for(unsigned int i=0; i<tracks.size(); i++){
			AddKinematicDataTrack(tracks[i], 28, 1.25);
		}
	}

	// DParticle
	if(hdvmf->GetCheckButton("particles")){
		vector<const DParticle*> particles;
		loop->Get(particles, hdvmf->GetFactoryTag("DParticle"));
		for(unsigned int i=0; i<particles.size(); i++){
			AddKinematicDataTrack(particles[i], 46, 1.00);
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
	if(name=="DParticle"){
		vector<const DParticle*> particles;
		if(loop)loop->Get(particles, tag.c_str());
		for(unsigned int i=0; i<particles.size(); i++)trks.push_back(particles[i]);
	}
	if(name=="DTrack"){
		vector<const DTrack*> tracks;
		if(loop)loop->Get(tracks, tag.c_str());
		for(unsigned int i=0; i<tracks.size(); i++)trks.push_back(tracks[i]);
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
		trkno<<setprecision(4)<<i+1;
		reconlabs["trk"][row]->SetText(trkno.str().c_str());
		
		double mass = trk->mass();
		if(fabs(mass-0.13957)<1.0E-4)type<<"pi";
		else if(fabs(mass-0.93827)<1.0E-4)type<<"proton";
		else if(fabs(mass-0.493677)<1.0E-4)type<<"K";
		else if(fabs(mass-0.000511)<1.0E-4)type<<"e";
		else type<<"q=";
		type<<(trk->charge()>0 ? "+":"-");
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

		// Get chisq and Ndof for DParticle or DTrack objects
		const DParticle *part=dynamic_cast<const DParticle*>(trk);
		const DTrack *track=dynamic_cast<const DTrack*>(trk);
		if(part){
			chisq_per_dof<<setprecision(4)<<part->chisq/part->Ndof;
			Ndof<<part->Ndof;
		}else if(track){
			chisq_per_dof<<setprecision(4)<<track->chisq/track->Ndof;
			Ndof<<track->Ndof;
		}else{
			chisq_per_dof<<"N/A";
			Ndof<<"N/A";
		}
		reconlabs["chisq/Ndof"][row]->SetText(chisq_per_dof.str().c_str());
		reconlabs["Ndof"][row]->SetText(Ndof.str().c_str());
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
		gset.points.push_back(step->origin);
	}
	
	// Push the graphics set onto the stack
	graphics.push_back(gset);
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
	if(dataname=="DParticle"){
		vector<const DParticle*> particles;
		loop->Get(particles, tag.c_str());
		if(index>=particles.size())return;
		q = particles[index]->charge();
		pos = particles[index]->position();
		mom = particles[index]->momentum();
		particles[index]->Get(cdchits);
		mass = particles[index]->mass();
	}

	if(dataname=="DTrack"){
		vector<const DTrack*> tracks;
		loop->Get(tracks, tag.c_str());
		if(index>=tracks.size())return;
		q = tracks[index]->charge();
		pos = tracks[index]->position();
		mom = tracks[index]->momentum();
		tracks[index]->Get(cdchits);
		mass = tracks[index]->mass();
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
		hit.second = fdchits[i]->dist;
		allhits.push_back(hit);
	}
}



