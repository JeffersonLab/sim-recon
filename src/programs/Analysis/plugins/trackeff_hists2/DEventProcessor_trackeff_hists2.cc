// $Id: DEventProcessor_trackeff_hists2.cc 2774 2007-07-19 15:59:02Z davidl $
//
//    File: DEventProcessor_trackeff_hists2.cc
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include <iostream>
#include <cmath>
using namespace std;

#include <TThread.h>

#include <TROOT.h>

#include "DEventProcessor_trackeff_hists2.h"
#include "DTrackingResolutionGEANT.h"

#include <JANA/JApplication.h>
#include <JANA/JEventLoop.h>

#include <DANA/DApplication.h>
#include <TRACKING/DMCThrown.h>
#include <TRACKING/DTrackCandidate.h>
#include <TRACKING/DTrackWireBased.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <PID/DParticle.h>
#include <FDC/DFDCGeometry.h>
#include <DVector2.h>
#include <particleType.h>
#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCPseudo.h>
#include <FDC/DFDCHit.h>

// Routine used to create our DEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_trackeff_hists2());
}
} // "C"


//------------------
// DEventProcessor_trackeff_hists2
//------------------
DEventProcessor_trackeff_hists2::DEventProcessor_trackeff_hists2()
{
	trk_ptr = &trk;

	//trkres = new DTrackingResolutionGEANT();

	pthread_mutex_init(&mutex, NULL);
	
}

//------------------
// ~DEventProcessor_trackeff_hists2
//------------------
DEventProcessor_trackeff_hists2::~DEventProcessor_trackeff_hists2()
{

}

//------------------
// init
//------------------
jerror_t DEventProcessor_trackeff_hists2::init(void)
{
	// Create TRACKING directory
	TDirectory *dir = (TDirectory*)gROOT->FindObject("TRACKING");
	if(!dir)dir = new TDirectoryFile("TRACKING","TRACKING");
	dir->cd();

	// Create Trees
	trkeff = new TTree("trkeff2","Tracking Efficiency");
	trkeff->Branch("E","track2",&trk_ptr);

	dir->cd("../");
	
	JParameterManager *parms = app->GetJParameterManager();

	DEBUG = 1;
	
	parms->SetDefaultParameter("TRKEFF:DEBUG", DEBUG);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_trackeff_hists2::brun(JEventLoop *loop, int runnumber)
{
	DApplication *dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
	const DGeometry *dgeom = dapp->GetDGeometry(runnumber);
	
	rt_thrown = new DReferenceTrajectory(dgeom->GetBfield());
	
	double dz, rmin, rmax;
	dgeom->GetCDCEndplate(CDCZmax, dz, rmin, rmax);

	double cdc_axial_length;
	dgeom->GetCDCAxialLength(cdc_axial_length);
	CDCZmin = CDCZmax-cdc_axial_length;

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_trackeff_hists2::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_trackeff_hists2::fini(void)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_trackeff_hists2::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DCDCTrackHit*> cdctrackhits;
	vector<const DFDCHit*> fdchits;
	vector<const DParticle*> particles;
	vector<const DMCThrown*> mcthrowns;
	vector<const DMCTrajectoryPoint*> mctrajpoints;
	
	// Bail quick on events with too many or too few CDC hits
	loop->Get(cdctrackhits);
	//if(cdctrackhits.size()>30 || cdctrackhits.size()<6)return NOERROR;

	// Bail quick on events with too many FDC hits
	loop->Get(fdchits);
	//if(fdchits.size()>30)return NOERROR;

	loop->Get(particles);
	loop->Get(mcthrowns);
	loop->Get(mctrajpoints);

	// Lock mutex
	pthread_mutex_lock(&mutex);
	
	// Get hit list for all throwns
	for(unsigned int i=0; i<mcthrowns.size(); i++){
		const DMCThrown *mcthrown = mcthrowns[i];
		
		// If there aren't enough DMCTrajectoryPoint objects then we will need to
		// get the LR information by swimming the thrown value ourself.
		bool use_rt_thrown = true; //mctrajpoints.size()<20;
		if(use_rt_thrown)rt_thrown->Swim(mcthrown->position(), mcthrown->momentum(), mcthrown->charge());

		// if this isn't a charged track, then skip it
		if(fabs(mcthrowns[i]->charge())==0.0)continue;

		// Momentum of thrown particle
		DVector3 pthrown = mcthrown->momentum();
		trk.pthrown = pthrown;

		// Initialize with the "not found" values
		trk.pfit.SetXYZ(0,0,0);
		trk.pfit_wire.SetXYZ(0,0,0);
		trk.pcan.SetXYZ(0,0,0);
		trk.trk_chisq=1.0E20;
		trk.trk_Ndof=1;
		trk.trk_chisq_wb=1.0E20;
		trk.trk_Ndof_wb=1;
		trk.delta_pt_over_pt=1.0E20;
		trk.delta_theta=1.0E20;
		trk.delta_phi=1.0E20;
		trk.isreconstructable = isReconstructable(mcthrown, mctrajpoints);
		trk.Nstereo = 0;
		trk.Ncdc = 0;
		trk.Nfdc = 0;
		trk.NLR_bad_stereo = 0;
		trk.NLR_bad = 0;
		trk.event = eventnumber;
		
		double fom_best = 1.0E8;

		// Loop over found/fit tracks
		for(unsigned int j=0; j<particles.size(); j++){
			const DParticle *particle = particles[j];
			
			// Get DTrackWireBased and DTrackCandidate objects for this DParticle
			vector<const DTrackWireBased*> tracks;
			particle->Get(tracks);
			const DTrackWireBased *track = tracks.size()==1 ? tracks[0]:NULL;
			vector<const DTrackCandidate*> trackcandidates;
			if(track)track->Get(trackcandidates);
			const DTrackCandidate *trackcandidate = trackcandidates.size()==1 ? trackcandidates[0]:NULL;

			// Copy momentum vectors to convenient local variables
			DVector3 pfit  = particle->momentum();
			DVector3 pfit_wire  = track ? track->momentum():DVector3(0,0,0);
			DVector3 pcan  = trackcandidate ? trackcandidate->momentum():DVector3(0,0,0);
			
			// Calculate residuals from momentum parameters from DParticle
			double delta_pt_over_pt = (pfit.Perp() - pthrown.Perp())/pthrown.Perp();
			double delta_theta = (pfit.Theta() - pthrown.Theta())*1000.0;
			double delta_phi = (pfit.Phi() - pthrown.Phi())*1000.0;

			// Formulate a figure of merit to decide if this fit track is closer to
			// the thrown track than the best one we found so far. We hardwire
			// dpt/pt=2%, dtheta=20mrad and dphi=20mrad for now.
			double fom = pow(delta_pt_over_pt/0.02, 2.0) + pow(delta_theta/20.0, 2.0) + pow(delta_phi/20.0, 2.0);
			if(fom<fom_best){
				fom_best = fom;
				
				trk.pfit = pfit;
				trk.pfit_wire = pfit_wire;
				trk.pcan = pcan;
				trk.trk_chisq = particle->chisq;
				trk.trk_Ndof = particle->Ndof;
				trk.trk_chisq_wb = track!=NULL ? track->chisq:1.0E6;
				trk.trk_Ndof_wb = track!=NULL ? track->Ndof:0;
				trk.delta_pt_over_pt = delta_pt_over_pt;
				trk.delta_theta = delta_theta;
				trk.delta_phi = delta_phi;

				// Get Nstereo, Ncdc, and Nfdc
				vector<const DCDCTrackHit*> cdchits;
				particle->Get(cdchits);
				trk.Nstereo = 0;
				for(unsigned int k=0; k<cdchits.size(); k++)if(cdchits[k]->wire->stereo!=0.0)trk.Nstereo++;
				trk.Ncdc = cdchits.size();
				vector<const DFDCPseudo*> fdchits;
				particle->Get(fdchits);
				trk.Nfdc = fdchits.size();
				
				// Get the number LR signs for all hits used on this track. We have to 
				// do this for the thrown track here so we get the list for the same hits
				// used in fitting this track.
				vector<int> LRthrown;
				vector<int> LRfit;
				vector<const DCoordinateSystem*> wires;
				for(unsigned int k=0; k<cdchits.size(); k++)wires.push_back(cdchits[k]->wire);
				for(unsigned int k=0; k<fdchits.size(); k++)wires.push_back(fdchits[k]->wire);
				if(use_rt_thrown){
					FindLR(wires, rt_thrown, LRthrown);
				}else{
					FindLR(wires, mctrajpoints, LRthrown);
				}
				FindLR(wires, particle->rt, LRfit);
				
				// Make sure the number of entries is the same for both LR vectors
				if(LRfit.size()!=LRthrown.size() || LRfit.size()!=wires.size()){
					_DBG_<<"LR vector sizes do not match! LRfit.size()="<<LRfit.size()<<" LRthrown.size()="<<LRthrown.size()<<" wires.size()="<<wires.size()<<endl;
					continue;
				}
				
				// count total number of incorrect LR choices and how many are stereo
				trk.NLR_bad_stereo = trk.NLR_bad = 0;
				for(unsigned int k=0; k<wires.size(); k++){
					if(LRfit[k] == LRthrown[k])continue;
					trk.NLR_bad++;
					bool is_stereo = (wires[k]->udir.Theta()*57.3)>2.0;
					if(is_stereo)trk.NLR_bad_stereo++;
				}
			}
		}		
		
		trkeff->Fill();
	}

	// Unlock mutex
	pthread_mutex_unlock(&mutex);
	
	return NOERROR;
}

//------------------
// isReconstructable
//------------------
bool DEventProcessor_trackeff_hists2::isReconstructable(const DMCThrown *mcthrown, vector<const DMCTrajectoryPoint*> &mctrajpoints)
{
	/// In order to test the efficiency of the finder/fitter, we must first
	/// determine whether a track is "reconstructible" or not. (See
	/// COMPASS-Note 2004-1 section 5.1 and
	/// Mankel Rep. Prog. Phys. 67 (2004) 553-622 section 2.5.1)
	///
	/// We do this here by checking if the "death" point of the 
	/// track is inside or outside of some inner volume of the detector.
	/// Specifically, if the death point is outside of a cylinder
	/// of radius 45cm and z-extent covering from the upstream end
	/// of the CDC to the downstream end of the 2nd FDC package.
	/// This ensures that it either passes through the first 2 layers
	/// of the outermost axial superlayer of the CDC or, the first
	/// 2 packages of the FDC.
	///
	/// Note that the beam hole is not considered here. It is assumed
	/// that efficiency plots will be made either with a cut on theta
	/// or plotted against theta in order to accomodate that area.

	const DMCTrajectoryPoint *mctraj_first = NULL;
	const DMCTrajectoryPoint *mctraj_last = NULL;
	for(unsigned int i=0; i<mctrajpoints.size(); i++){
		const DMCTrajectoryPoint *mctraj = mctrajpoints[i];
		if(mctraj->track != mcthrown->myid)continue;
		if(mctraj->primary_track != mcthrown->myid)continue; // redundant?
		if(mctraj_first==NULL)mctraj_first = mctraj;
		mctraj_last = mctraj;
	}
	
	if(mctraj_last!=NULL){
		double r = sqrt(pow((double)mctraj_last->x,2.0) + pow((double)mctraj_last->y,2.0));
		if(r>45.0)return true;
		if(mctraj_last->z<17.0)return true;
		if(mctraj_last->z<17.0)return true;
	}
	
	return false;
}

//------------------
// FindLR
//------------------
void DEventProcessor_trackeff_hists2::FindLR(vector<const DCoordinateSystem*> &wires, const DReferenceTrajectory *crt, vector<int> &LRhits)
{
	/// Fill the vector LRhits with +1 or -1 values indicating the side of each wire in the "wires"
	/// vector the given reference trajectory passed on.
	
	LRhits.clear();
	
	// This first bit is shameful. In order to use the DistToRT methods of the reference trajectory,
	// we have to cast away the const qualifier since those methods must modify the object
	DReferenceTrajectory *rt = const_cast<DReferenceTrajectory*>(crt);
	for(unsigned int i=0; i<wires.size(); i++){
		const DCoordinateSystem *wire = wires[i];

		DVector3 pos_doca, mom_doca;
		rt->DistToRT(wire);
		rt->GetLastDOCAPoint(pos_doca, mom_doca);
		DVector3 shift = wire->udir.Cross(mom_doca);
		double u = rt->GetLastDistAlongWire();
		DVector3 pos_wire = wire->origin + u*wire->udir;
		DVector3 pos_diff = pos_doca-pos_wire;
		
		LRhits.push_back(shift.Dot(pos_diff)<0.0 ? -1:1);
	}
}

//------------------
// FindLR
//------------------
void DEventProcessor_trackeff_hists2::FindLR(vector<const DCoordinateSystem*> &wires, vector<const DMCTrajectoryPoint*> &trajpoints, vector<int> &LRhits)
{		
	/// Fill the vector LRhits with +1 or -1 values indicating the side of each wire in the "wires"
	/// vector the particle swum by GEANT passed on according to the closest DMCTrajectoryPoint.

}

