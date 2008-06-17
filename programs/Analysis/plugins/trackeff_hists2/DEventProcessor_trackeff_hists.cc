// $Id: DEventProcessor_trackeff_hists.cc 2774 2007-07-19 15:59:02Z davidl $
//
//    File: DEventProcessor_trackeff_hists.cc
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include <iostream>
#include <cmath>
using namespace std;

#include <TThread.h>

#include "DEventProcessor_trackeff_hists.h"
#include "DTrackingResolutionGEANT.h"

#include <JANA/JApplication.h>
#include <JANA/JEventLoop.h>

#include <DANA/DApplication.h>
#include <TRACKING/DMCThrown.h>
#include <TRACKING/DTrackCandidate.h>
#include <TRACKING/DTrack.h>
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
	app->AddProcessor(new DEventProcessor_trackeff_hists());
}
} // "C"


//------------------
// DEventProcessor_trackeff_hists
//------------------
DEventProcessor_trackeff_hists::DEventProcessor_trackeff_hists()
{
	trk_ptr = &trk;

	trkres = new DTrackingResolutionGEANT();

	pthread_mutex_init(&mutex, NULL);
	
}

//------------------
// ~DEventProcessor_trackeff_hists
//------------------
DEventProcessor_trackeff_hists::~DEventProcessor_trackeff_hists()
{

}

//------------------
// init
//------------------
jerror_t DEventProcessor_trackeff_hists::init(void)
{
	// Create TRACKING directory
	TDirectory *dir = new TDirectory("TRACKING","TRACKING");
	dir->cd();

	// Create Trees
	trkeff = new TTree("trkeff","Tracking Efficiency");
	trkeff->Branch("E","track",&trk_ptr);

	dir->cd("../");
	
	JParameterManager *parms = app->GetJParameterManager();

	DEBUG = 1;
	
	parms->SetDefaultParameter("TRKEFF:DEBUG", DEBUG);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_trackeff_hists::brun(JEventLoop *loop, int runnumber)
{
	DApplication *dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
	const DGeometry *dgeom = dapp->GetDGeometry(runnumber);
	
	double dz, rmin, rmax;
	dgeom->GetCDCEndplate(CDCZmax, dz, rmin, rmax);

	double cdc_axial_length;
	dgeom->GetCDCAxialLength(cdc_axial_length);
	CDCZmin = CDCZmax-cdc_axial_length;
	
_DBG_<<"CDCZmin="<<CDCZmin<<"  CDCZmax="<<CDCZmax<<endl;

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_trackeff_hists::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_trackeff_hists::fini(void)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_trackeff_hists::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DCDCTrackHit*> cdctrackhits;
	vector<const DFDCHit*> fdchits;
	vector<const DTrack*> tracks;
	vector<const DMCThrown*> mcthrowns;
	vector<const DMCTrajectoryPoint*> mctrajpoints;
	
	// Bail quick on events with too many CDC hits
	loop->Get(cdctrackhits);
	if(cdctrackhits.size()>30 || cdctrackhits.size()<6)return NOERROR;

	// Bail quick on events with too many FDC hits
	loop->Get(fdchits);
	if(fdchits.size()>30)return NOERROR;

	loop->Get(tracks,"ALT1");
	loop->Get(mcthrowns);
	loop->Get(mctrajpoints);
	
	// Lock mutex
	pthread_mutex_lock(&mutex);

	// Get hit list for all throwns
	for(unsigned int i=0; i<mcthrowns.size(); i++){
		const DMCThrown *mcthrown = mcthrowns[i];
		
		// if this isn't a charged track, then skip it
		if(fabs(mcthrowns[i]->charge())==0.0)continue;

		// Momentum of thrown particle
		DVector3 pthrown = mcthrown->momentum();
		trk.pthrown = pthrown;
		
		// Get resolutions for this thrown track
		double pt_res, theta_res, phi_res;
		trkres->GetResolution(8 , pthrown, pt_res, theta_res, phi_res);

		// Initialize with the "no candidate" values
		trk.likelihood = 0.0;
		trk.chisq=1.0E20;
		trk.pt_pull=1.0E20;
		trk.theta_pull=1.0E20;
		trk.phi_pull=1.0E20;
		trk.isreconstructable = isReconstructable(mcthrown, mctrajpoints);

		// Loop over found/fit tracks and calculate chisq and likelihood
		// for each to be from this thrown track
		for(unsigned int j=0; j<tracks.size(); j++){
			const DTrack *track = tracks[j];
			
			// For places outside of our acceptance, the resolutions will be
			// returned as 0. For these (and any other case when a resolution
			// is 0) we automatically assume the track wasn't found.
			//
			// NOTE: There is a potential gotcha here. The resolution function
			// is set to return zeros for the resolutions for anything outside
			// of the ROOT histogram used to hold the resolutions. At this point,
			// the limits are 0-150 degrees in theta and 0.2-7.1 GeV/c in 
			// total momentum.
			if(DEBUG>10)_DBG_<<"pt_res="<<pt_res<<" theta_res="<<theta_res<<" phi_res="<<phi_res<<endl;
			if(DEBUG>10)_DBG_<<"p="<<pthrown.Mag()<<" theta="<<pthrown.Theta()*57.3<<" phi="<<pthrown.Phi()*57.3<<endl;
			if(pt_res==0.0 || theta_res==0.0 || phi_res==0.0)continue;

			DVector3 pfit = track->momentum();
			double pt_pull = (pfit.Perp() - pthrown.Perp())/pthrown.Perp()/pt_res;
			double theta_pull = (pfit.Theta() - pthrown.Theta())*1000.0/theta_res;
			double phi_pull = (pfit.Phi() - pthrown.Phi())*1000.0/phi_res;

			double chisq = (pow(pt_pull, 2.0) + pow(theta_pull, 2.0) + pow(phi_pull, 2.0))/3.0;
			double likelihood = exp(-chisq*3.0/2.0);
			if(DEBUG>10)_DBG_<<"chisq="<<chisq<<" likelihood="<<likelihood<<endl;
			if(chisq<trk.chisq){
				trk.likelihood = likelihood;
				trk.chisq = chisq;
				trk.pt_pull = pt_pull;
				trk.theta_pull = theta_pull;
				trk.phi_pull = phi_pull;
				trk.pfit = pfit;
				
				// Get # of CDC and FDC hits
				vector<const DCDCTrackHit*> cdchits;
				track->Get(cdchits);
				trk.Nstereo = 0;
				for(unsigned int k=0; k<cdchits.size(); k++)if(cdchits[k]->wire->stereo!=0.0)trk.Nstereo++;
				trk.Ncdc = cdchits.size();
				vector<const DFDCPseudo*> fdchits;
				track->Get(fdchits);
				trk.Nfdc = fdchits.size();
			}
		}
		
		if(DEBUG>5)if(trk.chisq>=3.0 && trk.pthrown.Theta()*57.3<150.0)_DBG_<<"Event:"<<eventnumber<<" thrown track "<<i<<"  trk.chisq="<<trk.chisq<<"  theta="<<trk.pthrown.Theta()*57.3<<" p="<<trk.pthrown.Mag()<<endl;
		if(DEBUG>0){
			//double sigma_pt_pull=0.811035, sigma_theta_pull=0.636278, sigma_phi_pull=1.65046;
			double sigma_pt_pull=1.0, sigma_theta_pull=1.0, sigma_phi_pull=1.0;
			double pt_pull = trk.pt_pull/sigma_pt_pull;
			double theta_pull = trk.theta_pull/sigma_theta_pull;
			double phi_pull = trk.phi_pull/sigma_phi_pull;
			double chisq_corrected = (pow(pt_pull,2.0) + pow(theta_pull,2.0) + pow(phi_pull, 2.0))/3.0;

			if(trk.isreconstructable && chisq_corrected>1000.0){

				_DBG_<<" Reconstructable event not found/fit: "<<eventnumber<<" (chisq_corrected="<<chisq_corrected<<")";
				if(DEBUG>1)cerr<<" pt_pull="<<pt_pull<<" theta_pull="<<theta_pull<<" phi_pull="<<phi_pull;
				cerr<<endl;
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
bool DEventProcessor_trackeff_hists::isReconstructable(const DMCThrown *mcthrown, vector<const DMCTrajectoryPoint*> &mctrajpoints)
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


