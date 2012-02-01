// $Id: DEventProcessor_trackeff_hists2.cc 2774 2007-07-19 15:59:02Z davidl $
//
//    File: DEventProcessor_trackeff_hists2.cc
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include "DEventProcessor_trackeff_hists2.h"

using namespace std;

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

  // Get the particle ID algorithms
	vector<const DParticleID *> locPIDAlgorithms;
	loop->Get(locPIDAlgorithms);
	if(locPIDAlgorithms.size() < 1){
		_DBG_<<"Unable to get a DParticleID object! NO PID will be done!"<<endl;
		return RESOURCE_UNAVAILABLE;
	}
	// Drop the const qualifier from the DParticleID pointer (I'm surely going to hell for this!)
	dPIDAlgorithm = const_cast<DParticleID*>(locPIDAlgorithms[0]);
  
	// Warn user if something happened that caused us NOT to get a dPIDAlgorithm object pointer
	if(!dPIDAlgorithm){
		_DBG_<<"Unable to get a DParticleID object! NO PID will be done!"<<endl;
		return RESOURCE_UNAVAILABLE;
	}

	use_rt_thrown = true; //mctrajpoints.size()<20;

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
	
	// Bail quick on events with too many or too few CDC hits
	vector<const DCDCTrackHit*> cdctrackhits;
	loop->Get(cdctrackhits);
	//if(cdctrackhits.size()>30 || cdctrackhits.size()<6)return NOERROR;

	// Bail quick on events with too many FDC hits
	vector<const DFDCHit*> fdchits;
	loop->Get(fdchits);
	//if(fdchits.size()>30)return NOERROR;

	vector<const DMCThrown*> mcthrowns;
	vector<const DMCTrajectoryPoint*> mctrajpoints;
	loop->Get(mctrajpoints);
	loop->Get(mcthrowns);

	// Lock mutex
	pthread_mutex_lock(&mutex);
	
	// Get hit list for all throwns
	for(unsigned int i=0; i<mcthrowns.size(); i++){
		const DMCThrown *mcthrown = mcthrowns[i];
		
		// If there aren't enough DMCTrajectoryPoint objects then we will need to
		// get the LR information by swimming the thrown value ourself.
		if(use_rt_thrown)rt_thrown->Swim(mcthrown->position(), mcthrown->momentum(), mcthrown->charge());

		// if this isn't a charged track, then skip it
		if(fabs(mcthrowns[i]->charge())==0.0)continue;

		// Momentum of thrown particle
		DVector3 pthrown = mcthrown->momentum();
		trk.pthrown = pthrown;
		trk.q_thrown = mcthrown->charge(); //make sure this value isn't bogus!!
		trk.PID_thrown = mcthrown->type;

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
		trk.delta_pt_over_pt_wire=1.0E20;
		trk.delta_theta_wire=1.0E20;
		trk.delta_phi_wire=1.0E20;
		trk.delta_pt_over_pt_can=1.0E20;
		trk.delta_theta_can=1.0E20;
		trk.delta_phi_can=1.0E20;
		trk.isreconstructable = isReconstructable(mcthrown, mctrajpoints);
		trk.Nstereo = 0;
		trk.Ncdc = 0;
		trk.Nfdc = 0;
		trk.NLR_bad_stereo = 0;
		trk.NLR_bad = 0;
		trk.event = eventnumber;
		trk.dTrackReconstructedFlag_Candidate = false;
		trk.dTrackReconstructedFlag_WireBased = false;
		trk.dTrackReconstructedFlag_TimeBased = false;
		trk.q_candidate = 0.0;
		trk.q_wirebased = 0.0;
		trk.q_timebased = 0.0;
		trk.PID_candidate = 0;
		trk.PID_wirebased = 0;
		trk.PID_timebased = 0;
		trk.PID_hypothesized = 0;
		trk.q_hypothesized = 0.0;
		trk.FOM_hypothesized = 0.0;

		bool locFoundFlag;
/*
		locFoundFlag = Search_ChargedTrackHypotheses(loop, eventnumber, mcthrown);
		if(locFoundFlag == false){
			trk.dTrackReconstructedFlag_TimeBased = false;
			locFoundFlag = Search_WireBasedTracks(loop, eventnumber, mcthrown);
			if(locFoundFlag == false){
				trk.dTrackReconstructedFlag_WireBased = false;
				locFoundFlag = Search_TrackCandidates(loop, eventnumber, mcthrown);
				if(locFoundFlag == false)
					trk.dTrackReconstructedFlag_Candidate = false;
			}
		}
*/

		locFoundFlag = Search_ChargedTrackHypotheses(loop, eventnumber, mcthrown);
		if(locFoundFlag == false)
			trk.dTrackReconstructedFlag_Candidate = false;

		locFoundFlag = Search_WireBasedTracks(loop, eventnumber, mcthrown);
		if(locFoundFlag == false)
			trk.dTrackReconstructedFlag_WireBased = false;

		locFoundFlag = Search_TrackCandidates(loop, eventnumber, mcthrown);
		if(locFoundFlag == false)
			trk.dTrackReconstructedFlag_Candidate = false;

		trkeff->Fill();
	}

	// Unlock mutex
	pthread_mutex_unlock(&mutex);
	
	return NOERROR;
}

bool DEventProcessor_trackeff_hists2::Search_ChargedTrackHypotheses(JEventLoop *loop, int eventnumber, const DMCThrown *mcthrown){
	vector<const DChargedTrackHypothesis*> locChargedTrackHypotheses;
	vector<const DMCTrajectoryPoint*> mctrajpoints;
	vector<const DCDCTrackHit*> cdctrackhits;
	vector<const DFDCHit*> fdchits;

	loop->Get(cdctrackhits);
	loop->Get(fdchits);
	loop->Get(locChargedTrackHypotheses);
	loop->Get(mctrajpoints);

	DVector3 pthrown = trk.pthrown;

	bool locFoundFlag = false;
	double fom_best = 1.0E8;

	trk.num_timebased = locChargedTrackHypotheses.size();

	// Loop over found/fit tracks
	for(unsigned int j=0; j<locChargedTrackHypotheses.size(); j++){
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrackHypotheses[j];
		const DTrackTimeBased *locTimeBasedTrack = locChargedTrackHypothesis->dTrackTimeBased;
	
		// Get DTrackWireBased and DTrackCandidate objects for this DTrackTimeBased
		vector<const DTrackWireBased*> tracks;
		locTimeBasedTrack->Get(tracks);
		const DTrackWireBased *track = tracks.size()==1 ? tracks[0]:NULL;
		vector<const DTrackCandidate*> trackcandidates;
		if(track)track->Get(trackcandidates);

		// Copy momentum vectors to convenient local variables
		DVector3 pfit  = locTimeBasedTrack->momentum();
	
		// Calculate residuals from momentum parameters from DTrackTimeBased
		double delta_pt_over_pt = (pfit.Perp() - pthrown.Perp())/pthrown.Perp();
		double delta_theta = (pfit.Theta() - pthrown.Theta())*1000.0;
		double delta_phi = (pfit.Phi() - pthrown.Phi())*1000.0;

		// Formulate a figure of merit to decide if this fit track is closer to
		// the thrown track than the best one we found so far. We hardwire
		// dpt/pt=2%, dtheta=20mrad and dphi=20mrad for now.
		double fom = pow(delta_pt_over_pt/0.02, 2.0) + pow(delta_theta/20.0, 2.0) + pow(delta_phi/20.0, 2.0);
		if(fom<fom_best){
			fom_best = fom;
			locFoundFlag = true;
			trk.PID_hypothesized = int(locChargedTrackHypothesis->dPID);
			trk.q_hypothesized = locChargedTrackHypothesis->charge();
			trk.FOM_hypothesized = locChargedTrackHypothesis->dFOM;

			trk.dTrackReconstructedFlag_TimeBased = true;

			trk.q_timebased = locTimeBasedTrack->charge();

			trk.PID_timebased = int(dPIDAlgorithm->IDTrack(locTimeBasedTrack->charge(), locTimeBasedTrack->mass()));

			trk.pfit = pfit;
			trk.trk_chisq = locTimeBasedTrack->chisq;
			trk.trk_Ndof = locTimeBasedTrack->Ndof;
			trk.delta_pt_over_pt = delta_pt_over_pt;
			trk.delta_theta = delta_theta;
			trk.delta_phi = delta_phi;

			// Get Nstereo, Ncdc, and Nfdc
			vector<const DCDCTrackHit*> cdchits;
			locTimeBasedTrack->Get(cdchits);
			trk.Nstereo = 0;
			for(unsigned int k=0; k<cdchits.size(); k++)if(cdchits[k]->wire->stereo!=0.0)trk.Nstereo++;
			trk.Ncdc = cdchits.size();
			vector<const DFDCPseudo*> fdchits;
			locTimeBasedTrack->Get(fdchits);
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
			FindLR(wires, locTimeBasedTrack->rt, LRfit);
		
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
	return locFoundFlag;	
}

bool DEventProcessor_trackeff_hists2::Search_WireBasedTracks(JEventLoop *loop, int eventnumber, const DMCThrown *mcthrown){
	vector<const DMCTrajectoryPoint*> mctrajpoints;
	vector<const DCDCTrackHit*> cdctrackhits;
	vector<const DFDCHit*> fdchits;
	vector<const DTrackWireBased*> tracks;

	loop->Get(cdctrackhits);
	loop->Get(fdchits);
	loop->Get(tracks);
	loop->Get(mctrajpoints);

	DVector3 pthrown = trk.pthrown;
	double fom_best = 1.0E8;

	trk.num_wirebased = tracks.size();

	bool locFoundFlag = false;
	// Loop over found/fit tracks
	for(unsigned int j=0; j<tracks.size(); j++){
		const DTrackWireBased *track = tracks[j];

		// Copy momentum vectors to convenient local variables
		DVector3 pfit_wire = track->momentum();
	
		// Calculate residuals from momentum parameters from DTrackTimeBased
		double delta_pt_over_pt_wire = (pfit_wire.Perp() - pthrown.Perp())/pthrown.Perp();
		double delta_theta_wire = (pfit_wire.Theta() - pthrown.Theta())*1000.0;
		double delta_phi_wire = (pfit_wire.Phi() - pthrown.Phi())*1000.0;

		// Formulate a figure of merit to decide if this fit track is closer to
		// the thrown track than the best one we found so far. We hardwire
		// dpt/pt=2%, dtheta=20mrad and dphi=20mrad for now.
		double fom = pow(delta_pt_over_pt_wire/0.02, 2.0) + pow(delta_theta_wire/20.0, 2.0) + pow(delta_phi_wire/20.0, 2.0);
		if(fom<fom_best){
			fom_best = fom;
			locFoundFlag = true;
			trk.dTrackReconstructedFlag_WireBased = true;

			trk.q_wirebased = track->charge();

			trk.PID_wirebased = int(dPIDAlgorithm->IDTrack(track->charge(), track->mass()));

			trk.pfit_wire = pfit_wire;
			trk.trk_chisq_wb = track->chisq;
			trk.trk_Ndof_wb = track->Ndof;
			trk.delta_pt_over_pt_wire = delta_pt_over_pt_wire;
			trk.delta_theta_wire = delta_theta_wire;
			trk.delta_phi_wire = delta_phi_wire;

			// Get Nstereo, Ncdc, and Nfdc
			vector<const DCDCTrackHit*> cdchits;
			track->Get(cdchits);
			trk.Nstereo = 0;
			for(unsigned int k=0; k<cdchits.size(); k++)if(cdchits[k]->wire->stereo!=0.0)trk.Nstereo++;
			trk.Ncdc = cdchits.size();
			vector<const DFDCPseudo*> fdchits;
			track->Get(fdchits);
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
			FindLR(wires, track->rt, LRfit);
		
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
	return locFoundFlag;	
}

bool DEventProcessor_trackeff_hists2::Search_TrackCandidates(JEventLoop *loop, int eventnumber, const DMCThrown *mcthrown){
	vector<const DMCTrajectoryPoint*> mctrajpoints;
	vector<const DCDCTrackHit*> cdctrackhits;
	vector<const DFDCHit*> fdchits;
	vector<const DTrackCandidate*> trackcandidates;

	loop->Get(cdctrackhits);
	loop->Get(fdchits);
	loop->Get(mctrajpoints);
	loop->Get(trackcandidates);

	DVector3 pthrown = trk.pthrown;
	double fom_best = 1.0E8;

	trk.num_candidates = trackcandidates.size();

	bool locFoundFlag = false;
	// Loop over found/fit tracks
	for(unsigned int j=0; j<trackcandidates.size(); j++){
		const DTrackCandidate *trackcandidate = trackcandidates[j];

		// Copy momentum vectors to convenient local variables
		DVector3 pcan  = trackcandidate->momentum();
	
		// Calculate residuals from momentum parameters from DTrackTimeBased
		double delta_pt_over_pt_can = (pcan.Perp() - pthrown.Perp())/pthrown.Perp();
		double delta_theta_can = (pcan.Theta() - pthrown.Theta())*1000.0;
		double delta_phi_can = (pcan.Phi() - pthrown.Phi())*1000.0;

		// Formulate a figure of merit to decide if this fit track is closer to
		// the thrown track than the best one we found so far. We hardwire
		// dpt/pt=2%, dtheta=20mrad and dphi=20mrad for now.
		double fom = pow(delta_pt_over_pt_can/0.02, 2.0) + pow(delta_theta_can/20.0, 2.0) + pow(delta_phi_can/20.0, 2.0);
		if(fom<fom_best){
			fom_best = fom;
			locFoundFlag = true;
			trk.dTrackReconstructedFlag_Candidate = true;
			trk.q_candidate = trackcandidate->charge();
			trk.PID_candidate = int(dPIDAlgorithm->IDTrack(trackcandidate->charge(), trackcandidate->mass()));
			trk.pcan = pcan;
			trk.delta_pt_over_pt_can = delta_pt_over_pt_can;
			trk.delta_theta_can = delta_theta_can;
			trk.delta_phi_can = delta_phi_can;

			// Get Nstereo, Ncdc, and Nfdc
			vector<const DCDCTrackHit*> cdchits;
			trackcandidate->Get(cdchits);
			trk.Nstereo = 0;
			for(unsigned int k=0; k<cdchits.size(); k++)if(cdchits[k]->wire->stereo!=0.0)trk.Nstereo++;
			trk.Ncdc = cdchits.size();
			vector<const DFDCPseudo*> fdchits;
			trackcandidate->Get(fdchits);
			trk.Nfdc = fdchits.size();
		}
	}
	return locFoundFlag;	
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

