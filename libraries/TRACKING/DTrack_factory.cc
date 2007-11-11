// $Id$
//
//    File: DTrack_factory.cc
// Created: Sun Apr  3 12:28:45 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include <math.h>

#include <TROOT.h>

#include <DVector3.h>
#include <DMatrix.h>

#include <JANA/JEventLoop.h>

#include "GlueX.h"
#include "DANA/DApplication.h"
#include "DMagneticFieldStepper.h"
#include "DTrackCandidate.h"
#include "DTrack_factory.h"
#include "CDC/DCDCTrackHit.h"
#include "FDC/DFDCPseudo.h"
#include "DReferenceTrajectory.h"
#include "DMCThrown.h"

#define NaN std::numeric_limits<double>::quiet_NaN()

// The GNU implementation of STL includes definitions of "greater" and "less"
// but the SunOS implementation does not. Since it is a bit of a pain to
// define this only for SunOS, we just define "greaterthan" and use it for
// all platforms/compilers. Note that this is essentially the same as the
// GNU definition from stl_function.h, except it does not derive from the
// templated "binary_function" class.
template<typename T>
class greaterthan{
	public: bool operator()(const T &a, const T &b) const {return a>b;}
};

#if 0
bool CDCTrkHitSort_C(DTrack_factory::cdc_hit_on_track_t const &hit1, DTrack_factory::cdc_hit_on_track_t const &hit2) {
	// These swim steps come from the same array so the addresses of the swim_step
	// structures should be sequential, starting with closest to the target.
	return (unsigned long)hit1.swim_step > (unsigned long)hit2.swim_step;
}
bool FDCTrkHitSort_C(DTrack_factory::fdc_hit_on_track_t const &hit1, DTrack_factory::fdc_hit_on_track_t const &hit2) {
	// These swim steps come from the same array so the addresses of the swim_step
	// structures should be sequential, starting with closest to the target.
	return (unsigned long)hit1.swim_step > (unsigned long)hit2.swim_step;
}
#endif

//------------------
// DTrack_factory   (Constructor)
//------------------
DTrack_factory::DTrack_factory()
{
	// This gets allocated in brun once the bfield is known
	tmprt=NULL;

	// Define target center
	target = new DCoordinateSystem();
	target->origin.SetXYZ(0.0, 0.0, 65.0);
	target->sdir.SetXYZ(1.0, 0.0, 0.0);
	target->tdir.SetXYZ(0.0, 1.0, 0.0);
	target->udir.SetXYZ(0.0, 0.0, 1.0);
	target->L = 30.0;

	MAX_HIT_DIST = 2.0; // cm
	DEBUG_HISTS = false;
	USE_CDC = true;
	USE_FDC_ANODE = true;
	USE_FDC_CATHODE = true;
	MAX_CHISQ_DIFF = 1.0E-6;
	MAX_FIT_ITERATIONS = 10;
	SIGMA_CDC = 0.0200;
	SIGMA_FDC_ANODE = 0.0200;
	SIGMA_FDC_CATHODE = 0.0200;
	CHISQ_MAX_RESI_SIGMAS = 20.0;
	CHISQ_GOOD_LIMIT = 2.0;
	LEAST_SQUARES_DP = 0.0001;
	LEAST_SQUARES_DX = 0.010;
	LEAST_SQUARES_MIN_HITS = 3;
	LEAST_SQUARES_MAX_E2NORM = 1.0E6;
	CANDIDATE_TAG = "CLASSIC";
	DEFAULT_STEP_SIZE = 0.5;
	MIN_CDC_HIT_PROB = 0.2;
	MAX_CDC_DOUBLE_HIT_PROB = 0.1;
	MIN_FDC_HIT_PROB = 0.2;
	MAX_FDC_DOUBLE_HIT_PROB = 0.1;
	TOF_MASS = 0.13957018;
	
	gPARMS->SetDefaultParameter("TRKFIT:MAX_HIT_DIST",				MAX_HIT_DIST);
	gPARMS->SetDefaultParameter("TRKFIT:DEBUG_HISTS",				DEBUG_HISTS);
	gPARMS->SetDefaultParameter("TRKFIT:USE_CDC",					USE_CDC);
	gPARMS->SetDefaultParameter("TRKFIT:USE_FDC_ANODE",			USE_FDC_ANODE);
	gPARMS->SetDefaultParameter("TRKFIT:USE_FDC_CATHODE",			USE_FDC_CATHODE);
	gPARMS->SetDefaultParameter("TRKFIT:MAX_CHISQ_DIFF",			MAX_CHISQ_DIFF);
	gPARMS->SetDefaultParameter("TRKFIT:MAX_FIT_ITERATIONS",		MAX_FIT_ITERATIONS);
	gPARMS->SetDefaultParameter("TRKFIT:SIGMA_CDC",					SIGMA_CDC);
	gPARMS->SetDefaultParameter("TRKFIT:SIGMA_FDC_ANODE",			SIGMA_FDC_ANODE);
	gPARMS->SetDefaultParameter("TRKFIT:SIGMA_FDC_CATHODE",		SIGMA_FDC_CATHODE);
	gPARMS->SetDefaultParameter("TRKFIT:CHISQ_MAX_RESI_SIGMAS",	CHISQ_MAX_RESI_SIGMAS);
	gPARMS->SetDefaultParameter("TRKFIT:CHISQ_GOOD_LIMIT",		CHISQ_GOOD_LIMIT);
	gPARMS->SetDefaultParameter("TRKFIT:LEAST_SQUARES_DP",		LEAST_SQUARES_DP);
	gPARMS->SetDefaultParameter("TRKFIT:LEAST_SQUARES_DX",		LEAST_SQUARES_DX);
	gPARMS->SetDefaultParameter("TRKFIT:LEAST_SQUARES_MIN_HITS",LEAST_SQUARES_MIN_HITS);
	gPARMS->SetDefaultParameter("TRKFIT:LEAST_SQUARES_MAX_E2NORM",LEAST_SQUARES_MAX_E2NORM);		
	gPARMS->SetDefaultParameter("TRKFIT:CANDIDATE_TAG",			CANDIDATE_TAG);
	gPARMS->SetDefaultParameter("TRKFIT:DEFAULT_STEP_SIZE",		DEFAULT_STEP_SIZE);
	gPARMS->SetDefaultParameter("TRKFIT:MIN_CDC_HIT_PROB",			MIN_CDC_HIT_PROB);
	gPARMS->SetDefaultParameter("TRKFIT:MAX_CDC_DOUBLE_HIT_PROB",	MAX_CDC_DOUBLE_HIT_PROB);
	gPARMS->SetDefaultParameter("TRKFIT:MIN_FDC_HIT_PROB",			MIN_FDC_HIT_PROB);
	gPARMS->SetDefaultParameter("TRKFIT:MAX_FDC_DOUBLE_HIT_PROB",	MAX_FDC_DOUBLE_HIT_PROB);
	gPARMS->SetDefaultParameter("TRKFIT:TOF_MASS",					TOF_MASS);
	
}

//------------------
// DTrack_factory   (Destructor)
//------------------
DTrack_factory::~DTrack_factory()
{
	if(tmprt)delete tmprt;
	for(unsigned int i=0; i<rtv.size(); i++)delete rtv[i];
}

//------------------
// init
//------------------
jerror_t DTrack_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrack_factory::brun(JEventLoop *loop, int runnumber)
{
	// Get pointer to DMagneticFieldMap field object
	//gPARMS->SetParameter("GEOM:BZ_CONST",  -2.0);	
	DApplication* dapp = dynamic_cast<DApplication*>(eventLoop->GetJApplication());
	bfield = dapp->GetBfield(); // temporary until new geometry scheme is worked out
	
	// (re)-allocate temporary DReferenceTrajectory used in
	// numerical derivatives
	if(tmprt)delete tmprt;
	tmprt = new DReferenceTrajectory(bfield);	
	
	// Set limits on CDC. (This should eventually come from HDDS)
	CDC_Z_MIN = 17.0;
	CDC_Z_MAX = CDC_Z_MIN + 175.0;
	hit_based = false;
	//cout<<__FILE__<<":"<<__LINE__<<"-------------- Helical Fitter TRACKING --------------"<<endl;
	cout<<__FILE__<<":"<<__LINE__<<"-------------- Least Squares TRACKING --------------"<<endl;
	
	if(DEBUG_HISTS){
		dapp->Lock();
		
		// Histograms may already exist. (Another thread may have created them)
		// Try and get pointers to the existing ones.
		cdcdoca_vs_dist = (TH2F*)gROOT->FindObject("cdcdoca_vs_dist");
		cdcdoca_vs_dist_vs_ring = (TH3F*)gROOT->FindObject("cdcdoca_vs_dist_vs_ring");
		fdcdoca_vs_dist = (TH2F*)gROOT->FindObject("fdcdoca_vs_dist");
		fdcu_vs_s = (TH2F*)gROOT->FindObject("fdcu_vs_s");
		dist_axial = (TH1F*)gROOT->FindObject("dist_axial");
		doca_axial = (TH1F*)gROOT->FindObject("doca_axial");
		dist_stereo = (TH1F*)gROOT->FindObject("dist_stereo");
		doca_stereo = (TH1F*)gROOT->FindObject("doca_stereo");
		chisq_final_vs_initial = (TH2F*)gROOT->FindObject("chisq_final_vs_initial");
		nhits_final_vs_initial = (TH2F*)gROOT->FindObject("nhits_final_vs_initial");
		Npasses = (TH1F*)gROOT->FindObject("Npasses");
		ptotal = (TH1F*)gROOT->FindObject("ptotal");
		residuals_cdc = (TH2F*)gROOT->FindObject("residuals_cdc");
		residuals_fdc_anode = (TH2F*)gROOT->FindObject("residuals_fdc_anode");
		residuals_fdc_cathode = (TH2F*)gROOT->FindObject("residuals_fdc_cathode");
		residuals_cdc_vs_s = (TH3F*)gROOT->FindObject("residuals_cdc_vs_s");
		residuals_fdc_anode_vs_s = (TH3F*)gROOT->FindObject("residuals_fdc_anode_vs_s");
		residuals_fdc_cathode_vs_s = (TH3F*)gROOT->FindObject("residuals_fdc_cathode_vs_s");
		initial_chisq_vs_Npasses = (TH2F*)gROOT->FindObject("initial_chisq_vs_Npasses");
		chisq_vs_pass = (TH2F*)gROOT->FindObject("chisq_vs_pass");
		dchisq_vs_pass = (TH2F*)gROOT->FindObject("dchisq_vs_pass");
		cdc_single_hit_prob = (TH1F*)gROOT->FindObject("cdc_single_hit_prob");
		cdc_double_hit_prob = (TH1F*)gROOT->FindObject("cdc_double_hit_prob");
		fdc_single_hit_prob = (TH1F*)gROOT->FindObject("fdc_single_hit_prob");
		fdc_double_hit_prob = (TH1F*)gROOT->FindObject("fdc_double_hit_prob");

		if(!cdcdoca_vs_dist)cdcdoca_vs_dist = new TH2F("cdcdoca_vs_dist","DOCA vs. DIST",300, 0.0, 1.2, 300, 0.0, 1.2);
		if(!cdcdoca_vs_dist_vs_ring)cdcdoca_vs_dist_vs_ring = new TH3F("cdcdoca_vs_dist_vs_ring","DOCA vs. DIST vs. ring",300, 0.0, 1.2, 300, 0.0, 1.2,23,0.5,23.5);
		if(!fdcdoca_vs_dist)fdcdoca_vs_dist = new TH2F("fdcdoca_vs_dist","DOCA vs. DIST",500, 0.0, 2.0, 500, 0.0, 2.0);
		if(!fdcu_vs_s)fdcu_vs_s = new TH2F("fdcu_vs_s","DOCA vs. DIST along wire",500, -60.0, 60.0, 500, -60.0, 60.0);
		if(!dist_axial)dist_axial = new TH1F("dist_axial","Distance from drift time for axial CDC wires",300,0.0,3.0);
		if(!doca_axial)doca_axial = new TH1F("doca_axial","DOCA of track for axial CDC wires",300,0.0,3.0);
		if(!dist_stereo)dist_stereo = new TH1F("dist_stereo","Distance from drift time for stereo CDC wires",300,0.0,3.0);
		if(!doca_stereo)doca_stereo = new TH1F("doca_stereo","DOCA of track for axial CDC wires",300,0.0,3.0);
		if(!chisq_final_vs_initial)chisq_final_vs_initial = new TH2F("chisq_final_vs_initial","Final vs. initial chi-sq.",200, 0.0, 10.0,50, 0.0, 2.5);
		if(!nhits_final_vs_initial)nhits_final_vs_initial = new TH2F("nhits_final_vs_initial","Final vs. initial nhits used in chi-sq.",30, -0.5, 29.5, 30, -0.5, 29.5);
		if(!Npasses)Npasses = new TH1F("Npasses","Npasses", 21, -0.5, 20.5);
		if(!ptotal)ptotal = new TH1F("ptotal","ptotal",1000, 0.1, 8.0);
		if(!residuals_cdc)residuals_cdc = new TH2F("residuals_cdc","Residuals in CDC",1000,-2.0,2.0,24,0.5,24.5);
		if(!residuals_fdc_anode)residuals_fdc_anode = new TH2F("residuals_fdc_anode","Residuals in FDC anodes",1000,-2.0,2.0,24,0.5,24.5);
		if(!residuals_fdc_cathode)residuals_fdc_cathode = new TH2F("residuals_fdc_cathode","Residuals in FDC cathode",1000,-2.0,2.0,24,0.5,24.5);
		if(!residuals_cdc_vs_s)residuals_cdc_vs_s = new TH3F("residuals_cdc_vs_s","Residuals in CDC vs. pathlength",1000,-2.0,2.0,24,0.5,24.5,100, 0.0, 800);
		if(!residuals_fdc_anode_vs_s)residuals_fdc_anode_vs_s = new TH3F("residuals_fdc_anode_vs_s","Residuals in FDC anode vs. pathlength",1000,-2.0,2.0,24,0.5,24.5,100, 0.0, 800);
		if(!residuals_fdc_cathode_vs_s)residuals_fdc_cathode_vs_s = new TH3F("residuals_fdc_cathode_vs_s","Residuals in FDC cathode vs. pathlength",1000,-2.0,2.0,24,0.5,24.5,100, 0.0, 800);
		if(!initial_chisq_vs_Npasses)initial_chisq_vs_Npasses = new TH2F("initial_chisq_vs_Npasses","Initial chi-sq vs. number of iterations", 25, -0.5, 24.5, 400, 0.0, 40.0);
		if(!chisq_vs_pass)chisq_vs_pass = new TH2F("chisq_vs_pass","Chi-sq vs. fit iteration", 25, -0.5, 24.5, 400, 0.0, 40.0);
		if(!dchisq_vs_pass)dchisq_vs_pass = new TH2F("dchisq_vs_pass","Change in Chi-sq vs. fit iteration", 25, -0.5, 24.5, 400, 0.0, 8.0);
		if(!cdc_single_hit_prob)cdc_single_hit_prob = new TH1F("cdc_single_hit_prob","Probability a CDC hit belongs to a track for all tracks",200,0.0,1.0);
		if(!cdc_double_hit_prob)cdc_double_hit_prob = new TH1F("cdc_double_hit_prob","Probability a CDC hit belongs to two tracks",200,0.0,1.0);
		if(!fdc_single_hit_prob)fdc_single_hit_prob = new TH1F("fdc_single_hit_prob","Probability a FDC hit belongs to a track for all tracks",200,0.0,1.0);
		if(!fdc_double_hit_prob)fdc_double_hit_prob = new TH1F("fdc_double_hit_prob","Probability a FDC hit belongs to two tracks",200,0.0,1.0);

		dapp->Unlock();		
	}

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrack_factory::evnt(JEventLoop *loop, int eventnumber)
{
	// Get the thrown values (this is just temporary)
	vector<const DMCThrown*> mcthrowns;
	loop->Get(mcthrowns);
	trackcandidates.clear();
	cdctrackhits.clear();
	fdctrackhits.clear();
	loop->Get(trackcandidates,CANDIDATE_TAG.c_str());
	loop->Get(cdctrackhits);
	loop->Get(fdctrackhits);

	// Assign hits to track candidates
	AssignHitsToCandidates();

	// Loop over track candidates
	for(unsigned int i=0; i<trackcandidates.size(); i++){
	
		// Copy the hits assigned to this track into the
		// cdchits_on_track and fdchits_on_track vectors
		cdchits_on_track.clear();
		for(unsigned int j=0; j<trackassignmentcdc.size(); j++){
			if(trackassignmentcdc[j]==(int)i)cdchits_on_track.push_back(cdctrackhits[j]);
		}
		fdchits_on_track.clear();
		for(unsigned int j=0; j<trackassignmentfdc.size(); j++){
			if(trackassignmentfdc[j]==(int)i)fdchits_on_track.push_back(fdctrackhits[j]);
		}
//_DBG_<<"Track Number: "<<eventnumber+i<<"  cdchits_on_track="<<cdchits_on_track.size()<<"  fdchits_on_track="<<fdchits_on_track.size()<<endl;
		// Fit the track
		DTrack *track = FitTrack(rtv[i], trackcandidates[i]->id,  NULL);

		// If fit is successful, then store the track
		if(track)_data.push_back(track);
	}

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTrack_factory::fini(void)
{
	return NOERROR;
}

//------------------
// FitTrack
//------------------
DTrack* DTrack_factory::FitTrack(DReferenceTrajectory* rt, int candidateid, const DMCThrown *thrown)
{
	/// Fit a track candidate
//_DBG__;
//_DBG_<<"cdchits_on_track.size="<<cdchits_on_track.size()<<"  fdchits_on_track.size="<<fdchits_on_track.size()<<endl;
//_DBG_<<"cdctrackhits.size="<<cdctrackhits.size()<<"  fdctrackhits.size="<<fdctrackhits.size()<<endl;
	// Get starting position and momentum from reference trajectory
	DVector3 start_pos = rt->swim_steps[0].origin;
	DVector3 start_mom = rt->swim_steps[0].mom;
	DVector3 pos = start_pos;
	DVector3 mom = start_mom;

	// Fit the track and get the results in the form of the
	// position/momentum
	DVector3 vertex_pos=pos; // to hold fitted values on return
	DVector3 vertex_mom=mom; // to hold fitted values on return
	double initial_chisq=1.0E6, chisq=1.0E6, last_chisq;
	int Niterations;
	for(Niterations=0; Niterations<MAX_FIT_ITERATIONS; Niterations++){
		last_chisq = chisq;
		fit_status_t status = LeastSquares(pos, mom, rt, vertex_pos, vertex_mom, chisq);
		if(status != FIT_OK){
			vertex_pos = pos;
			vertex_mom = mom;
			rt->Swim(vertex_pos, vertex_mom);
			break;
		}

		if(Niterations==0)initial_chisq = chisq;
		if(DEBUG_HISTS){
			chisq_vs_pass->Fill(Niterations+1, chisq);
			dchisq_vs_pass->Fill(Niterations+1, last_chisq-chisq);			
		}
		if(fabs(last_chisq-chisq) < MAX_CHISQ_DIFF)break;
		if(chisq>1.0E4)break;
		
		// If the fit succeeded (which it had to if we got here)
		// then rt should already reflect the track at
		// vertex_pos, vertex_mom.
		pos = vertex_pos;
		mom = vertex_mom;
	}
	
	// At this point we must decided whether the fit succeeded or not.
	// We'll consider the fit a success if:
	//
	// 1. We got through at least one iteration in the above loop
	// 2. The chi-sq is less than CHISQ_GOOD_LIMIT
	// 3. MAX_FIT_ITERATIONS is zero (for debugging)
	bool fit_succeeded = false;
	if(Niterations>0)fit_succeeded = true;
	if(chisq<=CHISQ_GOOD_LIMIT)fit_succeeded = true;
	if(MAX_FIT_ITERATIONS==0)fit_succeeded = true;
	if(!fit_succeeded)return NULL;
	
	// Find point of closest approach to target and use parameters
	// there for vertex position and momentum
	//double s;
	//rt->DistToRT(target, &s);
	//rt->GetLastDOCAPoint(vertex_pos, vertex_mom);
//double dp_over_p = fabs(thrown->p-vertex_mom.Mag())/thrown->p;
//start_pos.Print();
//vertex_pos.Print();
//start_mom.Print();
//vertex_mom.Print();
//_DBG_<<"final chi-sq:"<<chisq<<"   DOF="<<chisqv.size()<<"   Niterations="<<Niterations<<" dp_over_p="<<dp_over_p<<endl;

	// Create new DTrack object and initialize parameters with those
	// from track candidate
	DTrack *track = new DTrack;
	track->q			= rt->q;
	track->p			= vertex_mom.Mag();
	track->theta	= vertex_mom.Theta();
	track->phi		= vertex_mom.Phi();
	if(track->phi<0.0)track->phi+=2.0*M_PI;
	track->x			= vertex_pos.X();
	track->y			= vertex_pos.Y();
	track->z			= vertex_pos.Z();
	track->candidateid = candidateid;
	track->chisq	= chisq;
	track->rt		= rt;
	
	// Fill in DKinematicData part
	track->setMass(0.0);
	track->setMomentum(vertex_mom);
	track->setPosition(vertex_pos);
	track->setCharge(rt->q);
	//track->setErrorMatrix(last_covariance);

	// Fill debugging histos if requested
	if(DEBUG_HISTS){
		Npasses->Fill(Niterations);
		initial_chisq_vs_Npasses->Fill(Niterations, initial_chisq);
		FillDebugHists(rt, vertex_pos, vertex_mom);
	}

	return track;
}

//------------------
// AssignHitsToCandidates
//------------------
void DTrack_factory::AssignHitsToCandidates(void)
{
	/// Sort all CDC and FDC hits by which track candidate they most likely belong to.
	///
	/// For each hit this calculates the probabilities that it came
	/// from each of the track candidates. Hits are then assigned
	/// to a candidate according to highest probability. Hits
	/// are also dropped if the probability of belonging to 
	/// more than one track is sufficiently large.
	
	// Allocate more DReferenceTrajectory objects if needed.
	// These each have a large enough memory footprint that
	// it causes noticable performance problems if we allocated
	// and deallocated them every event. Therefore, we allocate
	// when needed, but recycle them on the next event.
	// They are delete in the factory deconstructor.
	while(rtv.size()<trackcandidates.size())rtv.push_back(new DReferenceTrajectory(bfield));
	
	// Loop over track candidates
	vector<vector<double> > cdcprobs(cdctrackhits.size());
	vector<vector<double> > fdcprobs(fdctrackhits.size());
	for(unsigned int i=0; i<trackcandidates.size(); i++){
		const DTrackCandidate *tc = trackcandidates[i];
		DReferenceTrajectory *rt = rtv[i];

		// Get starting position and momentum from track candidate
		// and swim the intial reference trajectory
		const DVector3 &pos = tc->position();
		const DVector3 &mom = tc->momentum();
		rt->Swim(pos, mom, tc->charge());
		
		// Get probabilities of each CDC hit being on this track
		vector<double> cdcprob;
		GetCDCTrackHitProbabilities(rt, cdcprob);
		for(unsigned int j=0; j<cdcprob.size(); j++)cdcprobs[j].push_back(cdcprob[j]);

		// Get probabilities of each FDC hit being on this track
		vector<double> fdcprob;
		GetFDCTrackHitProbabilities(rt, fdcprob);
		for(unsigned int j=0; j<fdcprob.size(); j++)fdcprobs[j].push_back(fdcprob[j]);
//_DBG_<<"Candidate "<<i<<endl;
//for(unsigned int j=0; j<cdcprob.size(); j++)_DBG_<<"  cdcprob["<<j<<"] = "<<cdcprob[j]<<endl;
//for(unsigned int j=0; j<fdcprob.size(); j++)_DBG_<<"  fdcprob["<<j<<"] = "<<fdcprob[j]<<endl;
	}
	
	// Now we need to loop over the hits and decide which track,
	// if any, this hit should belong to.
	// First, the CDC...
	trackassignmentcdc.clear();
	for(unsigned int i=0; i<cdcprobs.size(); i++){ // Loop over hits
		// The order of the elements in the probs
		// vector corresponds to the order of the hits
		// in the cdctrackhits vector so it must be
		// preserved. We make a copy of the vector
		// to sort.
		vector<double> prob = cdcprobs[i];
		sort(prob.begin(), prob.end(), greaterthan<double>()); // sort in *descending* order

		// get most likely probability
		double p_single = prob.size()>0 ? prob[0]:0.0;
		if(DEBUG_HISTS)cdc_single_hit_prob->Fill(p_single);
		if(p_single<MIN_CDC_HIT_PROB){
			trackassignmentcdc.push_back(-1); // not consistent with any track
			continue;
		}
		
		// Check probability of belonging to more than one track
		double p_double = prob.size()>1 ? prob[1]*p_single:0.0;
		if(DEBUG_HISTS)cdc_double_hit_prob->Fill(p_double);
		if(p_double>MAX_CDC_DOUBLE_HIT_PROB){
			trackassignmentcdc.push_back(-1); // consistent with more than one track
			continue;
		}

		// Find index of hit with highest probability
		for(unsigned int j=0; j<cdcprobs[i].size(); j++){
			if(cdcprobs[i][j] == p_single){
				trackassignmentcdc.push_back(j);
				break;
			}
		}
	}

	// ... Now , the FDC
	trackassignmentfdc.clear();
	for(unsigned int i=0; i<fdcprobs.size(); i++){ // Loop over hits
		// The order of the elements in the probs
		// vector corresponds to the order of the hits
		// in the fdctrackhits vector so it must be
		// preserved. We make a copy of the vector
		// to sort.
		vector<double> prob = fdcprobs[i];
		sort(prob.begin(), prob.end(), greaterthan<double>()); // sort in *descending* order

		// get most likely probability
		double p_single = prob.size()>0 ? prob[0]:0.0;
		if(DEBUG_HISTS)fdc_single_hit_prob->Fill(p_single);
		if(p_single<MIN_FDC_HIT_PROB){
			trackassignmentfdc.push_back(-1); // not consistent with any track
			continue;
		}
		
		// Check probability of belonging to more than one track
		double p_double = prob.size()>1 ? prob[1]*p_single:0.0;
		if(DEBUG_HISTS)fdc_double_hit_prob->Fill(p_double);
		if(p_double>MAX_FDC_DOUBLE_HIT_PROB){
			trackassignmentfdc.push_back(-1); // consistent with more than one track
			continue;
		}

		// Find index of hit with highest probability
		for(unsigned int j=0; j<fdcprobs[i].size(); j++){
			if(fdcprobs[i][j] == p_single){
				trackassignmentfdc.push_back(j);
				break;
			}
		}
	}
}

//------------------
// GetCDCTrackHitProbabilities
//------------------
void DTrack_factory::GetCDCTrackHitProbabilities(DReferenceTrajectory *rt, vector<double> &prob)
{
	/// Determine the probability that for each CDC hit that it came from the track with the given trajectory.
	///
	/// This will calculate a probability for each CDC hit that
	/// it came from the track represented by the given
	/// DReference trajectory. The probability is based on
	/// the residual between the distance of closest approach
	/// of the trajectory to the wire and the drift time.

	// Calculate beta of particle assuming its a pion for now.
	double beta = 1.0/sqrt(1.0+TOF_MASS*TOF_MASS/rt->swim_steps[0].mom.Mag2());
	
	// The error on the residual. This should include the
	// error from measurement,track parameters, and multiple 
	// scattering.
	double sigma = sqrt(pow(SIGMA_CDC,2.0) + pow(0.4000,2.0));

	prob.clear();
	for(unsigned int j=0; j<cdctrackhits.size(); j++){
		const DCDCTrackHit *hit = cdctrackhits[j];
		
		// Find the DOCA to this wire
		double s;
		double doca = rt->DistToRT(hit->wire, &s);

		// Distance using drift time
		// NOTE: Right now we assume pions for the TOF
		// and a constant drift velocity of 55um/ns
		double tof = s/(beta*3E10*1E-9);
		double dist = (hit->tdrift - tof)*55E-4;
		
		// Residual
		double resi = dist - doca;

		// Use an un-normalized gaussian so that for a residual
		// of zero, we get a probability of 1.0.
		prob.push_back(finite(resi) ? exp(-pow(resi/sigma,2.0)):0.0);
	}
}

//------------------
// GetFDCTrackHitProbabilities
//------------------
void DTrack_factory::GetFDCTrackHitProbabilities(DReferenceTrajectory *rt, vector<double> &prob)
{
	/// Determine the probability that for each FDC hit that it came from the track with the given trajectory.
	///
	/// This will calculate a probability for each FDC hit that
	/// it came from the track represented by the given
	/// DReference trajectory. The probability is based on
	/// the residual between the distance of closest approach
	/// of the trajectory to the wire and the drift time
	/// and the distance along the wire.

	// Calculate beta of particle assuming its a pion for now.
	double beta = 1.0/sqrt(1.0+TOF_MASS*TOF_MASS/rt->swim_steps[0].mom.Mag2());
	
	// The error on the residual. This should include the
	// error from measurement,track parameters, and multiple 
	// scattering.
	double sigma_anode = sqrt(pow(SIGMA_FDC_ANODE,2.0) + pow(0.4000,2.0));
	double sigma_cathode = sqrt(pow(SIGMA_FDC_CATHODE,2.0) + pow(0.300,2.0));
	
	prob.clear();
	for(unsigned int j=0; j<fdctrackhits.size(); j++){
		const DFDCPseudo *hit = fdctrackhits[j];
		
		// Find the DOCA to this wire
		double s;
		double doca = rt->DistToRT(hit->wire, &s);

		// Distance using drift time
		// NOTE: Right now we assume pions for the TOF
		// and a constant drift velocity of 55um/ns
		double tof = s/(beta*3E10*1E-9);
		double dist = (hit->time - tof)*55E-4;
		
		// Residual
		double resi = dist - doca;
		
		// Use an un-normalized gaussian so that for a residual
		// of zero, we get a probability of 1.0.
		double p = finite(resi) ? exp(-pow(resi/sigma_anode,2.0)):0.0;

		// Cathode
		double u = rt->GetLastDistAlongWire();
		resi = u - hit->s;

		// Same as for the anode. We multiply the
		// probabilities to get a total probability
		// based on both the anode and cathode hits.
		p *= finite(resi) ? exp(-pow(resi/sigma_cathode,2.0)):0.0;
		
		prob.push_back(p);
	}
}

//------------------
// ChiSq
//------------------
double DTrack_factory::ChiSq(double q, DMatrix &state, const swim_step_t *start_step, DReferenceTrajectory *rt)
{
	DVector3 vdir = start_step->sdir.Cross(start_step->mom);
	vdir.SetMag(1.0);

	DVector3 pos =   start_step->origin
						+ state[state_x ][0]*start_step->sdir
						+ state[state_v ][0]*vdir;
	DVector3 mom =   state[state_px][0]*start_step->sdir
						+ state[state_py][0]*start_step->tdir
						+ state[state_pz][0]*start_step->udir;

	if(rt)rt->Swim(pos,mom);

	return ChiSq(q, pos, mom,rt);
}

//------------------
// ChiSq
//------------------
double DTrack_factory::ChiSq(double q, const DVector3 &pos, const DVector3 &mom, DReferenceTrajectory *rt)
{
	// Swim a reference trajectory using the state defined by
	// "state" at "start_step" if one is not provided.
	bool own_rt = false;
	if(!rt){
		rt = new DReferenceTrajectory(bfield, q, NULL, 0, DEFAULT_STEP_SIZE);
		rt->Swim(pos,mom);
		own_rt = true;
	}
	
	double chisq = ChiSq(q, rt);

	// If we created the reference trajectory, then delete
	if(own_rt)delete rt;
	
	return chisq;
}

//------------------
// ChiSq
//------------------
double DTrack_factory::ChiSq(double q, DReferenceTrajectory *rt)
{
	// Clear chisq and sigma vector
	chisqv.clear();
	sigmav.clear();
	
	if(rt->Nswim_steps<3)return 1.0E6;
	
	// Calculate particle beta
	double beta = 1.0/sqrt(1.0+TOF_MASS*TOF_MASS/rt->swim_steps[0].mom.Mag2()); // assume this is a pion for now. This should eventually come from outer detectors
	
	// Add CDC hits (if any)
	for(unsigned int i=0; i<cdchits_on_track.size(); i++){
		const DCDCTrackHit *hit = cdchits_on_track[i];
		const DCoordinateSystem *wire = hit->wire;
		
		// Distance of closest approach for track to wire
		double s;
		double doca = rt->DistToRT(wire, &s);
			
		// Distance from drift time. Hardwired for simulation for now
		double dist, sigma;
		if(hit_based){
			dist=0.400; // Middle of the tube is best guess when no time info is present
			sigma = 2.0*dist/3.464; // sigma of flat distr. is W/sqrt(12) where W is total width
		}else{
			// Calculate time of flight (TOF) so we can subtract it
			double tof = s/(beta*3E10*1E-9);
			dist = (hit->tdrift-tof)*55E-4;
			sigma = SIGMA_CDC; // 200 um
//sigma += 0.010*s/100.0; // add 50 um per meter due to MULS
		}

		// NOTE: Sometimes we push nan or large values on here
		double resi = dist - doca;

		chisqv.push_back(USE_CDC ? resi:NaN);
		sigmav.push_back(sigma);
	}
	
	// Add FDC hits (if any)
	for(unsigned int i=0; i<fdchits_on_track.size(); i++){
		const DFDCPseudo *hit = fdchits_on_track[i];
		const DCoordinateSystem *wire = hit->wire;
		
		// Distance of closest approach for track to wire
		double s;
		double doca = rt->DistToRT(wire,&s);

		// Distance from drift time. Hardwired for simulation for now.
		double dist, sigma;
		if(hit_based){
			dist = 1.116/4.0; // 1/4 of wire spacing is best guess when no time info is present
			sigma = 2.0*dist/3.464; // sigma of flat distr. is W/sqrt(12) where W is total width
		}else{
			// Calculate time of flight (TOF) so we can subtract it
			double tof = s/(beta*3E10*1E-9);
			dist = (hit->time-tof)*55E-4;
			sigma = SIGMA_FDC_ANODE; // 200 um
//sigma += 0.010*s/100.0; // add 50 um per meter due to MULS
		}

		// NOTE: Sometimes we push nan or large values on here
		double resi = dist - doca;
		chisqv.push_back(USE_FDC_ANODE ? resi:NaN);
		sigmav.push_back(sigma);
		
		// For the FDC we also have a measurement along the wire
		// which we include as a separate measurement
		double u = rt->GetLastDistAlongWire();
		resi = u - hit->s;
		sigma = SIGMA_FDC_CATHODE;
//sigma += 0.010*s/100.0; // add 50 um per meter due to MULS
		chisqv.push_back(USE_FDC_CATHODE ? resi:NaN);
		sigmav.push_back(sigma); // 200 um
		
	}
	
	// Include distance from target as element of chi-sq. This will
	// need to be changed for off-beamline vertices.
	double s;
	double d = rt->DistToRT(target, &s);
	chisqv.push_back(d);
	sigmav.push_back(0.1); // 1 mm
	
	// Filter out all but the best hits from the chisq
	// If we have at least 10 hits, then throw away the
	// worst 10% of the hits.
	vector<double> tmp;
	for(unsigned int i=0; i<chisqv.size(); i++) tmp.push_back(fabs(chisqv[i]/sigmav[i]));
	sort(tmp.begin(),tmp.end());
	unsigned int index = tmp.size();
	if(index>10)index = (unsigned int)((float)index*0.9);
	double chisq_max_resi_sigmas = tmp[index];

	// Add "good" hits together to get the chi-squared
	double chisq = 0;
	Ngood_chisq_hits=0.0;
	for(unsigned int i=0; i<chisqv.size(); i++){
		// Check if this hit is close enough to the RT to be on
		// the track. This should be a function of both the distance
		// along the track and the errors of the measurement in 3D.
		// For now, we use a single distance which may be sufficient.
		if(!finite(chisqv[i]))continue;
		if(fabs(chisqv[i]/sigmav[i])>CHISQ_MAX_RESI_SIGMAS)continue;
		if(fabs(chisqv[i]/sigmav[i])>chisq_max_resi_sigmas)continue;

		chisq+=pow(chisqv[i]/sigmav[i], 2.0);
		Ngood_chisq_hits += 1.0;
	}
	chisq/=Ngood_chisq_hits;

	return chisq;
}

//------------------
// LeastSquares
//------------------
DTrack_factory::fit_status_t DTrack_factory::LeastSquares(DVector3 &start_pos, DVector3 &start_mom, DReferenceTrajectory *rt, DVector3 &vertex_pos, DVector3 &vertex_mom, double &chisq)
{
	// Determine the best fit of the track using the least squares method
	// described by R. Mankel Rep. Prog. Phys. 67 (2004) 553-622 pg 565
	const int Nparameters = 5;
	double deltas[Nparameters];

	// For fitting, we want to define a coordinate system very similar to the
	// Reference Trajectory coordinate system. The difference is that we want
	// the position to be defined in a plane perpendicular to the momentum.
	// The RT x-direction is in this plane, but the momentum itself lies
	// somewhere in the Y/Z plane so that neither Y nor Z makes a good choice
	// for the second postion dimension. We will call the second coordinate in 
	// the perpendicular plane "v" which is the coordinate along the axis that
	// is perpendicular to both the x-direction and the momentum direction.
	swim_step_t start_step = rt->swim_steps[0];
	DVector3 vdir = start_step.sdir.Cross(start_step.mom);
	vdir.SetMag(1.0);
	
	// Since we define the state in the X/V plane, we need to propagate the
	// particle from the given starting position up to this plane.
	DVector3 pos = start_pos;
	DVector3 mom = start_mom;
//	DMagneticFieldStepper stepper(bfield, rt->q);
//	if(stepper.SwimToPlane(pos, mom, start_step.origin, start_step.udir)){
//return FIT_FAILED;
//	}

	// Define the particle state in the reference trajectory coordinate
	// system at the start of the RT
	DVector3 pos_diff = pos - start_step.origin;
	DMatrix state(5,1);
	switch(Nparameters){
		case 5: state[state_v	][0] = pos_diff.Dot(vdir);
		case 4: state[state_x	][0] = pos_diff.Dot(start_step.sdir);
		case 3: state[state_pz	][0] = mom.Dot(start_step.udir);
		case 2: state[state_py	][0] = mom.Dot(start_step.tdir);
		case 1: state[state_px	][0] = mom.Dot(start_step.sdir);
	}

	// Create reference trajectory to use in calculating derivatives
	tmprt->Swim(pos,mom, rt->q);

	// Best-guess
	chisq = ChiSq(rt->q, tmprt);
	vector<double> resi = chisqv;
	vector<double> errs = sigmav;

	// Because we have a non-linear system, we must take the derivatives
	// numerically.
	//
	// Note that in the calculations of the deltas below,
	// the change in state should be set first and the value
	// of deltas[...] calculated from that. See Numerical
	// Recipes in C 2nd ed. section 5.7 ppg. 186-189.

	// dpx : tweak by +/- 0.01
	DMatrix state_dpx = state;
	state_dpx[state_px][0] += LEAST_SQUARES_DP;
	deltas[state_px] = state_dpx[state_px][0] - state[state_px][0];
	ChiSq(rt->q, state_dpx, &start_step,tmprt);
	vector<double> resi_dpx_hi = chisqv;
	vector<double> &resi_dpx_lo = resi;

	// dpy : tweak by +/- 0.01
	DMatrix state_dpy = state;
	state_dpy[state_py][0] += LEAST_SQUARES_DP;
	deltas[state_py] = state_dpy[state_py][0] - state[state_py][0];
	ChiSq(rt->q, state_dpy, &start_step,tmprt);
	vector<double> resi_dpy_hi = chisqv;
	vector<double> &resi_dpy_lo = resi;

	// dpz : tweak by +/- 0.01
	DMatrix state_dpz = state;
	state_dpz[state_pz][0] += LEAST_SQUARES_DP;
	deltas[state_pz] = state_dpz[state_pz][0] - state[state_pz][0];
	ChiSq(rt->q, state_dpz, &start_step,tmprt);
	vector<double> resi_dpz_hi = chisqv;
	vector<double> &resi_dpz_lo = resi;

	// dx : tweak by +/- 0.01
	DMatrix state_dx = state;
	state_dx[state_x][0] += LEAST_SQUARES_DX;
	deltas[state_x] = state_dx[state_x][0] - state[state_x][0];
	if(Nparameters>=4)ChiSq(rt->q, state_dx, &start_step,tmprt);
	vector<double> resi_dx_hi = chisqv;
	vector<double> &resi_dx_lo = resi;

	// dv : tweak by +/- 0.01
	DMatrix state_dv = state;
	state_dv[state_v][0] += LEAST_SQUARES_DX;
	deltas[state_v] = state_dv[state_v][0] - state[state_v][0];
	if(Nparameters>=5)ChiSq(rt->q, state_dv, &start_step,tmprt);
	vector<double> resi_dv_hi = chisqv;
	vector<double> &resi_dv_lo = resi;
	
	// Make a list of "clean" hits. Ones with reasonably
	// small, residuals that are not "nan" for the
	// best-guess as well as the tweaked  cases.
	vector<bool> good;
	unsigned int Ngood=0;
	unsigned int Nhits = resi.size();
	for(unsigned int i=0; i<Nhits; i++){
		double res;
		double max=CHISQ_MAX_RESI_SIGMAS;
		res =        resi[i]; if(!finite(res) || fabs(res)>max){good.push_back(false); continue;}
		res = resi_dpx_hi[i]; if(!finite(res) || fabs(res)>max){good.push_back(false); continue;}
		res = resi_dpy_hi[i]; if(!finite(res) || fabs(res)>max){good.push_back(false); continue;}
		res = resi_dpz_hi[i]; if(!finite(res) || fabs(res)>max){good.push_back(false); continue;}

		res = resi_dx_hi[i]; if(!finite(res) || fabs(res)>max){good.push_back(false); continue;}
		res = resi_dv_hi[i]; if(!finite(res) || fabs(res)>max){good.push_back(false); continue;}

		good.push_back(true);
		Ngood++;
	}
	if(Ngood<LEAST_SQUARES_MIN_HITS){
		//cout<<__FILE__<<":"<<__LINE__<<" Bad number of good distance calculations!"<<endl;
		return FIT_FAILED;
	}

	// Build "F" matrix of derivatives
	DMatrix F(Ngood,Nparameters);
	for(unsigned int i=0,j=0; j<Nhits; j++){
		if(!good[j])continue; // skip bad hits
		switch(Nparameters){
			// Note: This is a funny way to use a switch!
			case 5: F[i][state_v ] = (resi_dv_hi[j]-resi_dv_lo[j])/deltas[state_v];
			case 4: F[i][state_x ] = (resi_dx_hi[j]-resi_dx_lo[j])/deltas[state_x];
			case 3: F[i][state_pz] = (resi_dpz_hi[j]-resi_dpz_lo[j])/deltas[state_pz];
			case 2: F[i][state_py] = (resi_dpy_hi[j]-resi_dpy_lo[j])/deltas[state_py];
			case 1: F[i][state_px] = (resi_dpx_hi[j]-resi_dpx_lo[j])/deltas[state_px];
		}
		i++;
	}
	DMatrix Ft(DMatrix::kTransposed, F);

	// V is a diagonal matrix of the measurement errors. In
	// principle, we could fold in the multiple scattering
	// here, but for now, we don't. This is filled in the
	// loop below.
	DMatrix V(Ngood,Ngood);
	V.Zero();
	
	// Measurement vector. This contains the residuals between
	// DOCAs and DISTs.
	DMatrix m(Ngood,1);
	for(unsigned int i=0,j=0; j<Nhits; j++){
		if(!good[j])continue; // skip bad hits
		m[i][0] = -resi[j]; // drift time is already subtracted in ChiSq(...)
		V[i][i] = pow(errs[j], 2.0);
		i++;
	}
	DMatrix Vinv(DMatrix::kInverted, V);
	DMatrix B(DMatrix::kInverted, Ft*Vinv*F);
	
	// If the inversion failed altogether then the invalid flag
	// will be set on the matrix. In these cases, were dead.
	if(!B.IsValid())return FIT_FAILED;

	// The "B" matrix happens to be the covariance matrix of the
	// state parameters. A problem sometimes occurs where one or
	// more elements of B are very large. This tends to happen
	// when a column of F is essentially all zeros making the
	// matrix un-invertable.What we should really do in these
	// cases is check beforehand and drop the bad column(s)
	// before trying to invert. That will add complication that
	// I don't want to introduce quite yet. What we do now
	// is check for it and punt rather than return a nonsensical
	// value.
	if(B.E2Norm() > LEAST_SQUARES_MAX_E2NORM)return FIT_FAILED;
	
	// Copy the B matrix into last_covariance to later copy into DTrack
	last_covariance = B;

	// Calculate step direction and magnitude	
	DMatrix delta_state = B*Ft*Vinv*m;

	// The following is based on Numerical Recipes in C 2nd Ed.
	// ppg. 384-385.
	//
	// We now have the "direction" by which to step in the-d
	// parameter space as well as an amplitude by which to
	// step in "delta_state". A potential problem is that
	// we can over-step such that we end up at a worse
	// chi-squared value. To address this, we try taking 
	// the full step, but if we end up at a worse chi-sq
	// then we backtrack and take a smaller step until
	// we see the chi-sq improve.
	double lambda = 1.0;
	int Ntrys = 0;
	double new_chisq;
	int max_trys = 4;
	for(; Ntrys<max_trys; Ntrys++){

		DMatrix new_state(5,1);
		for(int i=0; i<Nparameters; i++)new_state[i][0] = state[i][0] + delta_state[i][0]*lambda;

		// Calculate initial particle position/momentum.
		vertex_pos =     start_step.origin
							+ new_state[state_x ][0]*start_step.sdir
							+ new_state[state_v ][0]*vdir;
		vertex_mom =     new_state[state_px][0]*start_step.sdir
							+ new_state[state_py][0]*start_step.tdir
							+ new_state[state_pz][0]*start_step.udir;

		// Re-swim reference trajectory using these parameters and find Chi-sq
		rt->Swim(vertex_pos, vertex_mom);
		new_chisq = ChiSq(rt->q, rt);
		
		// If we're at a lower chi-sq then we're done
//_DBG_<<"chisq="<<chisq<<"  new_chisq="<<new_chisq<<" nhits="<<chisqv.size()<<"  lambda="<<lambda<<endl;
		if(new_chisq-chisq < 0.01)break;
		
		// Chi-sq was increased, try a smaller step on the next iteration
		lambda/=2.0;
	}
	
	// If we failed to find a better Chi-Sq above, maybe we were looking 
	// in the wrong direction(??) Try looking in the opposite direction.
	if(Ntrys>=max_trys){
		lambda = -lambda;
		for(Ntrys=0; Ntrys<max_trys; Ntrys++){

			DMatrix new_state(5,1);
			for(int i=0; i<Nparameters; i++)new_state[i][0] = state[i][0] + delta_state[i][0]*lambda;

			// Calculate initial particle position/momentum.
			vertex_pos =     start_step.origin
								+ new_state[state_x ][0]*start_step.sdir
								+ new_state[state_v ][0]*vdir;
			vertex_mom =     new_state[state_px][0]*start_step.sdir
								+ new_state[state_py][0]*start_step.tdir
								+ new_state[state_pz][0]*start_step.udir;

			// Re-swim reference trajectory using these parameters and find Chi-sq
			rt->Swim(vertex_pos, vertex_mom);
			new_chisq = ChiSq(rt->q, rt);

			// If we're at a lower chi-sq then we're done
//_DBG_<<"chisq="<<chisq<<"  new_chisq="<<new_chisq<<"  lambda="<<lambda<<endl;
		if(new_chisq-chisq < 0.01)break;

			// Chi-sq was increased, try a smaller step on the next iteration
			lambda*=2.0;
		}
	}

	// If we failed to make a step to a smaller chi-sq then signal
	// that the fit failed completely.
	if(Ntrys>=max_trys)return FIT_FAILED;
	
	chisq = new_chisq;

	return FIT_OK;
}

//------------------
// FillDebugHists
//------------------
void DTrack_factory::FillDebugHists(DReferenceTrajectory *rt, DVector3 &vertex_pos, DVector3 &vertex_mom)
{
	//vertex_mom.SetMagThetaPhi(6.0, 17.2*M_PI/180.0, 90.0*M_PI/180.0);
	//vertex_pos.SetXYZ(0.0,0.0,65.0);
	//rt->Swim(vertex_pos, vertex_mom);
	ptotal->Fill(vertex_mom.Mag());

	// Calculate particle beta
	double beta = 1.0/sqrt(1.0+TOF_MASS*TOF_MASS/vertex_mom.Mag2()); // assume this is a pion for now. This should eventually come from outer detectors

	for(unsigned int j=0; j<cdchits_on_track.size(); j++){
		const DCDCTrackHit *hit = cdchits_on_track[j];
		const DCDCWire *wire = hit->wire;
		
		// Distance of closest approach for track to wire
		double s;
		double doca = rt->DistToRT(wire, &s);
			
		// Distance from drift time. Hardwired for simulation for now
		double tof = s/(beta*3E10*1E-9);
		double dist = (hit->tdrift-tof)*55E-4;

		// NOTE: Sometimes this could be nan
		double resi = dist - doca;
		if(!finite(resi))continue;
		
		// Fill histos
		residuals_cdc->Fill(resi, wire->ring);
		residuals_cdc_vs_s->Fill(resi, wire->ring, s);

		cdcdoca_vs_dist->Fill(dist, doca);
		cdcdoca_vs_dist_vs_ring->Fill(dist, doca, wire->ring);
		if(wire->stereo==0.0){
			dist_axial->Fill(dist);
			doca_axial->Fill(doca);
		}else{
			dist_stereo->Fill(dist);
			doca_stereo->Fill(doca);
		}
	}
	
	for(unsigned int j=0; j<fdchits_on_track.size(); j++){
		const DFDCPseudo *hit = fdchits_on_track[j];
		const DFDCWire *wire = hit->wire;

		// Distance of closest approach for track to wire
		double s;
		double doca = rt->DistToRT(wire,&s);

		double tof = s/(beta*3E10*1E-9);
		double dist = (hit->time-tof)*55E-4;

		// NOTE: Sometimes this is nan
		double resi = dist - doca;
		if(finite(resi)){
			fdcdoca_vs_dist->Fill(dist, doca);
			residuals_fdc_anode->Fill(resi, wire->layer);
			residuals_fdc_anode_vs_s->Fill(resi, wire->layer,s);
		}
		
		double u = rt->GetLastDistAlongWire();
		resi = u - hit->s;
		if(finite(resi)){
			fdcu_vs_s->Fill(u, hit->s);
			residuals_fdc_cathode->Fill(resi, wire->layer);
			residuals_fdc_cathode_vs_s->Fill(resi, wire->layer,s);
		}
	}

	// Find chi-sq for both initial and final values
#if 0
	double chisq_final = ChiSq(rt->q, rt);
	double nhits_final = Ngood_chisq_hits;
	DVector3 mom;
	mom.SetMagThetaPhi(tc->p, tc->theta, tc->phi);
	DVector3 pos(0.0, 0.0, tc->z_vertex);
	rt->Swim(pos, mom);
	double chisq_initial = ChiSq(rt->q, rt);
	double nhits_initial = Ngood_chisq_hits;
	chisq_final_vs_initial->Fill(chisq_initial,chisq_final);
	nhits_final_vs_initial->Fill(nhits_initial, nhits_final);
#endif
}

//------------------
// toString
//------------------
const string DTrack_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	GetNrows();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("row: q:       p:   theta:   phi:    x:    y:    z:");
	
	for(unsigned int i=0; i<_data.size(); i++){

		DTrack *track = _data[i];

		printnewrow();
		
		printcol("%x", i);
		printcol("%+d", (int)track->q);
		printcol("%3.3f", track->p);
		printcol("%1.3f", track->theta);
		printcol("%1.3f", track->phi);
		printcol("%2.2f", track->x);
		printcol("%2.2f", track->y);
		printcol("%2.2f", track->z);

		printrow();
	}
	
	return _table;
}




//======================================================================
//======================================================================
//           UNUSED CODEBELOW HERE
//======================================================================
//======================================================================

#if 0
//------------------
// KalmanFilter
//------------------
void DTrack_factory::KalmanFilter(DMatrix &state
											, DMatrix &P
											, DReferenceTrajectory *rt
											, DVector3 &vertex_pos
											, DVector3 &vertex_mom)
{
	/// Apply the Kalman filter to the current track starting
	/// with the given initial state and the reference trajectory.

	// State vector for Kalman filter. This is kept as a DMatrix so
	// we can use the ROOT linear algebra package. We index the
	// elements via enum for readability.
	// We start at the last point on the reference trajectory so we can
	// swim in towards the target to get the state at the vertex.
	swim_step_t *step = &rt->swim_steps[rt->Nswim_steps-1];
	DMatrix state(5,1);
	state[state_px	][0] = step->mom.Dot(step->xdir);
	state[state_py	][0] = step->mom.Dot(step->ydir);
	state[state_pz	][0] = step->mom.Dot(step->zdir);
	state[state_x	][0] = 0.0;
	state[state_y	][0] = 0.0;

	// The covariance matrix for the state vector. Assume all of
	// the parameters are independent. The values are initialized
	// using the resolutions of the parameters obtained from 
	// track candidates. NOTE: right now, they are just guesses!
	// This will need to be changed!
	DMatrix P(5,5);
	P.UnitMatrix();
	P[state_px ][state_px ] = pow(0.1*fabs(state[state_px][0]) + 0.010 , 2.0);
	P[state_py ][state_py ] = pow(0.1*fabs(state[state_py][0]) + 0.010 , 2.0);
	P[state_pz ][state_pz ] = pow(0.1*fabs(state[state_pz][0]) + 0.010 , 2.0);
	P[state_x  ][state_x  ] = pow(5.0 , 2.0);  // cm
	P[state_y  ][state_y  ] = pow(5.0 , 2.0);  // cm

	// The A matrix propagates the state of the particle from one
	// point to the next. This is recalculated for each.
	DMatrix A(5, 5);
	
	// The covariance matrix representing the measurement errors.
	// For the CDC we have just one measurement "r". For the FDC
	// we have two: "r" and "w", the distance along the wire.
	DMatrix R_cdc(1,1);
	DMatrix R_fdc(2,2);
	
	// The Q matrix represents the "process noise" in our case,
	// this is where multiple scattering comes in. It will need
	// to be calculated at each step to include M.S. for the
	// materials traversed since the last step. To start with,
	// we will set this to zero indicating no M.S., just to keep
	// it as a place holder.
	DMatrix Q(5,5);
	Q = 0.0*A;
	
	// K is the Kalman "gain matrix"
	DMatrix K(5,5);
	
	// H is the matrix that converts the state vector values
	// into (predicted) measurement values.
	DMatrix H_cdc(1,5);
	DMatrix H_fdc(2,5);

	// DMagneticFieldStepper is used to swim the particle
	// through the magnetic field. Instantiate it here so
	// it doesn't need to be done each time through the loop
	// below.
	double q = rt->q;
	swim_step_t step_prev = rt->swim_steps[rt->Nswim_steps-1];
	DMagneticFieldStepper stepper(bfield);

	// Loop over hits. They should already be in order of
	// distance along RT from largest to smallest.
	for(unsigned int i=0; i<cdchits_on_track.size(); i++){
		cdc_hit_on_track_t &hit = cdchits_on_track[i];
		const DCoordinateSystem *wire = hit.cdchit->wire;

		// The curent hit corresponds to the next measurement point
		// i.e. the place we need to project the state to.
		// This is a bit tricky. The problem is in defining exactly
		// where the measurement point is. The DOCA points for a set
		// of various states will not necessarily correspond to the
		// same "s" position on the RT. What we do is the following:
		// 1.) Find the DOCA point of the RT itself. Call it "s". 
		// 2.) Find the point on the swum track on the plane normal
		//     to the RT at point "s".
		// 3.) Use the point on the swum track found in 2.) to 
		//     define the state of the swum track at "s"
		//
		// The "measurement" value we'll use will be the actual DOCA
		// of the track(i.e. not the distance of the point
		// found in 2.) to the hit).
		
		// 1.)
		//------------------
		// During the hit finding, the value of "s" of the DOCA was
		// stored in the cdc_hit_on_track_t structure. What we really
		// want is a swim step at the DOCA point. We get this by
		// setting up the stepper at the swim step from the 
		// hit, setting the step size to the difference
		// of the "s" where the step is at and the "s" where we 
		// want it to be and then taking one step. This may seem
		// a bit convoluted, but it at least guarantees we swim
		// in the same way as everywhere else. Note that the momentum
		// and charge of the particle at hit.swim_step correspond
		// to the particle being swum forward from the target.
		swim_step_t step_docaRT;
		double ds = hit.s - hit.swim_step->s;
		stepper.SetStepSize(ds);
		stepper.SetStartingParams(q, &hit.swim_step->pos, &hit.swim_step->mom);
		stepper.Step(&step_docaRT.pos);
		stepper.GetDirs(step_docaRT.xdir, step_docaRT.ydir, step_docaRT.zdir);
		stepper.GetMomentum(step_docaRT.mom);
		step_docaRT.Ro = stepper.GetRo();
		step_docaRT.s = hit.s;
		
		// Project to get a'priori prediction of state at "s" of step_docaRT
		DMatrix mystate = state;
		double dist = ProjectStateBackwards(q, mystate, &step_prev, &step_docaRT, wire);
		if(!finite(dist))continue; // skip this hit if there's a problem
#if 1
		// Since the relation between the current state and
		// the next is non-linear, we have to estimate it as
		// a linear transformation. This is not trivally calculated
		// especially considering the non-uniform material
		// the track passes through as it goes through the detector.
		// Thus, we do this numerically by tweaking each of the state
		// parameters one at a time and re-swimming to see how it
		// affects the state parameters at the end of the step.
		double deltas[5], dists[5];
		
		// dpx : tweak by 10% + 1MeV/c
		DMatrix state_dpx = state;
		deltas[state_px] = 0.001 + 0.10*state_dpx[state_px][0];
		state_dpx[state_px][0] += deltas[state_px];
		dists[state_px] = ProjectStateBackwards(q, state_dpx, &step_prev, &step_docaRT, wire);
		if(!finite(dists[state_px]))continue; // skip this hit if there's a problem

		// dpy : tweak by 10% + 1MeV/c
		DMatrix state_dpy = state;
		deltas[state_py] = 0.001 + 0.10*state_dpy[state_py][0];
		state_dpy[state_py][0] += deltas[state_py];
		dists[state_py] = ProjectStateBackwards(q, state_dpy, &step_prev, &step_docaRT, wire);
		if(!finite(dists[state_py]))continue; // skip this hit if there's a problem

		// dpz : tweak by 10% + 1MeV/c
		DMatrix state_dpz = state;
		deltas[state_pz] = 0.001 + 0.10*state_dpz[state_pz][0];
		state_dpz[state_pz][0] += deltas[state_pz];
		dists[state_pz] = ProjectStateBackwards(q, state_dpz, &step_prev, &step_docaRT, wire);
		if(!finite(dists[state_pz]))continue; // skip this hit if there's a problem

		// dx : tweak by 250 microns
		DMatrix state_dx = state;
		deltas[state_x] = 0.0250;
		state_dx[state_x][0] += deltas[state_x];
		dists[state_x] = ProjectStateBackwards(q, state_dx, &step_prev, &step_docaRT, wire);
		if(!finite(dists[state_x]))continue; // skip this hit if there's a problem

		// dy : tweak by 500 microns
		DMatrix state_dy = state;
		deltas[state_y] = 0.0250;
		state_dy[state_y][0] += deltas[state_y];
		dists[state_y] = ProjectStateBackwards(q, state_dy, &step_prev, &step_docaRT, wire);
		if(!finite(dists[state_y]))continue; // skip this hit if there's a problem
		
		// Phew!! Believe it or not, it has taken a *very* long
		// time just to get to this point.
		
		// Calculate "A" and "H" matrices as well as residual
		for(int j=0; j<5; j++){
			A[j][state_px ] = (state_dpx[j][0] - state[j][0])/deltas[state_px];
			A[j][state_py ] = (state_dpy[j][0] - state[j][0])/deltas[state_py];
			A[j][state_pz ] = (state_dpz[j][0] - state[j][0])/deltas[state_pz];
			A[j][state_x  ] = ( state_dx[j][0] - state[j][0])/deltas[state_x ];
			A[j][state_y  ] = ( state_dy[j][0] - state[j][0])/deltas[state_y ];
			
			H_cdc[0][j] = (dists[j]-dist)/deltas[j];
cout<<__FILE__<<":"<<__LINE__<<" dists[j]="<<dists[j]<<" dist="<<dist<<" H_cdc[]="<<H_cdc[0][j]<<endl;
			//H_cdc[0][j]=0.0;
		}

		// For Monte Carlo, dist = 55um*tdrift
		DMatrix z_minus_h(1,1);
		double z = hit.cdchit->tdrift*55.0E-4;
		z=0.0; // Hit-based tracking is equivalent to tdrift=0
		z_minus_h[0][0] = z - dist;
		
		// Measurement error on drift distance
		R_cdc[0][0] = pow(200.0E-4,2.0); // use constant 400um measurement error for now
		
		// Hmmm... it seems like V should be a unit matrix, but then, whay have it?
		DMatrix V(DMatrix::kUnit, R_cdc);
		
		// W corresponds to process noise covariance. We set it to the
		// identity matrix, but it is only used to multiply Q which
		// is a NULL matrix so this doesn't really matter.
		DMatrix W(DMatrix::kUnit, Q);

		// Apply Kalman filter equations to include this measurement point.
		KalmanStep(mystate, P, z_minus_h, A, H_cdc, Q, R_cdc, W, V);
#endif
		state = mystate;

		// Make current step_docaRT next iterations step_prev.
		step_prev = step_docaRT;
	}
	
	// Swim particle backwards to the beamline
	DVector3 pos = step_prev.pos
						+ state[state_x ][0]*step_prev.xdir
						+ state[state_y ][0]*step_prev.ydir;
	DVector3 mom =   state[state_px][0]*step_prev.xdir
						+ state[state_py][0]*step_prev.ydir
						+ state[state_pz][0]*step_prev.zdir;
	mom = -mom;
	stepper.SetStepSize(0.1);
	stepper.SetStartingParams(-q, &pos, &mom);
_DBG_;
pos.Print();
mom.Print();
	double min_diff2 = 1.0E6;
	for(int i=0; i<1000; i++){  // for loop guarantees no infinite loop
		stepper.GetPosition(vertex_pos);
		stepper.GetMomentum(vertex_mom);
		stepper.Step(&pos);
		double diff2 = pos.Perp2();
cout<<__FILE__<<":"<<__LINE__<<" diff2="<<diff2<<endl;
		if(diff2>min_diff2)break;
		min_diff2 = diff2;
	}
	vertex_mom = -vertex_mom; // flip momentum back to forward direction
	
cout<<__FILE__<<":"<<__LINE__<<"Before p="<<rt->swim_steps[0].mom.Mag()<<"  ChiSq="<<ChiSq(rt->q, rt->swim_steps[0].pos, rt->swim_steps[0].mom)<<endl;
_DBG_;
cout<<__FILE__<<":"<<__LINE__<<"After p="<<vertex_mom.Mag()<<"  Chisq="<<ChiSq(q, vertex_pos, vertex_mom)<<endl;

}

//------------------
// KalmanStep
//------------------
void DTrack_factory::KalmanStep(	DMatrix &x,
											DMatrix &P,
											DMatrix &z_minus_h,
											DMatrix &A,
											DMatrix &H,
											DMatrix &Q,
											DMatrix &R,
											DMatrix &W,
											DMatrix &V)
{
	/// Update the state vector x and its covariance matrix P
	/// to include one more measurement point using the Extended
	/// Kalman filter equations. The symbols follow the notation
	/// of Welch and Bishop TR95-041. The notation used in Mankel
	/// Rep. Prog. Phys. 67 (2004) 553-622 uses the following:
	///
	///  Mankel    W.B.  Description
	///  ------   -----  -------------
	///   F         A    Propagation matrix
	///   C         P    Covariance of state
	///   R         -
	///   H         H    Projection matrix (state onto measurement space)
	///   K         K    "Gain" matrix
	///   V         R    Covariance of measurement
	///   Q         Q    Process noise
	///
	///
	/// Upon entry, x should already represent the state projected
	/// up to this measurement point. P, however, should contain
	/// the covariance at the previous measurement point. This is
	/// because we're using the EKF and the calculation of x is
	/// done through a non-linear function. P, however is propagated
	/// using linear transformations.
	///
	/// The value of z_minus_h should be the "residual" between the 
	/// measurement vector and the predicted measurement vector.
	/// This also contains a non-linear transformation (in the
	/// predicted measurement).
	///
	/// The values of A, H, W, and V are all Jacobian matrices.
	/// 
	/// Q and R represent process noise and measurement covariance
	/// respectively. Under ideal conditions, both of these could
	/// be NULL matrices.
	
	DMatrix At(DMatrix::kTransposed, A);
	DMatrix Ht(DMatrix::kTransposed, H);
	DMatrix Vt(DMatrix::kTransposed, V);
	DMatrix Wt(DMatrix::kTransposed, W);

_DBG_;	
DMatrix HPHt = H*P*Ht;
A.Print();
P.Print();
H.Print();
HPHt.Print();

DMatrix Atmp(DMatrix::kUnit, A);
A=Atmp;
At=Atmp;

	P = A*P*At + W*Q*Wt;
	
	DMatrix B(DMatrix::kInverted, H*P*Ht + V*R*Vt);
cout<<__FILE__<<":"<<__LINE__<<" B[0][0] = "<<B[0][0]<<"  sqrt(1/B[0][0])="<<sqrt(1/B[0][0])<<endl;

	DMatrix K = P*Ht*B;
K.Print();
	DMatrix I(DMatrix::kUnit, P);
x.Print();
	x = x + K*(z_minus_h);
x.Print();
	P = (I - K*H)*P;
}

//------------------
// ProjectStateBackwards
//------------------
double DTrack_factory::ProjectStateBackwards(double q, DMatrix &state
											, const swim_step_t *stepRT_prev
											, const swim_step_t *stepRT
											, const DCDCWire *wire)
{
	/// Swim a particle defined by q/state at the given swim step
	/// (stepRT_prev) on the RT backwards to the point corresponding to
	/// "s" at stepRT. The values in the state vector are updated to  
	/// reflect the state of the particle at the "s" value of stepRT.
	/// Note that stepRT_prev should be at a larger "s" value than
	/// stepRT.

	// First, set up a stepper at the starting point. We need to flip
	// both q and mom to swim backwards
	swim_step_t step;
	step.pos =    stepRT_prev->pos
					+ state[state_x ][0]*stepRT_prev->xdir
					+ state[state_y ][0]*stepRT_prev->ydir;
	step.mom =	  state[state_px][0]*stepRT_prev->xdir
					+ state[state_py][0]*stepRT_prev->ydir
					+ state[state_pz][0]*stepRT_prev->zdir;
	step.mom = -step.mom;
	DMagneticFieldStepper stepper(bfield, -q, &step.pos, &step.mom);
	// Swim until we go just past the point at stepRT
	double min_diff2 = 1.0E6;
	double dz_dphi;
	for(int i=0; i<500; i++){ // for-loop avoids infinite loop
		// Get values of current (and presumably closest so far) step
		stepper.GetDirs(step.xdir, step.ydir, step.zdir);
		stepper.GetMomentum(step.mom);
		stepper.GetPosition(step.pos);
		dz_dphi = stepper.Getdz_dphi();
		step.Ro = stepper.GetRo();

		// Step to next point and check if we've started to go away.
		DVector3 mypos;
		stepper.Step(&mypos);
		DVector3 diff = mypos - stepRT->pos;
		double diff2 = diff.Mag2();
		if(diff2>min_diff2)break;
		min_diff2 = diff2;
	}
	
	// At this point, xdir,ydir,zdir,pos,Ro, and dz_dphi should
	// contain the step on the track we're swimming that is
	// closest to stepRT point.
	
	// Find the point on the current track at "s". In other words,
	// we need to find the point on the track we're swimming that
	// is at the same "s" value as stepRT. This is so we can
	// define the state of the track's x/y values in terms of the RT.
	//
	// This "s" value is represented by a plane defined by the origin
	// of the stepRT with the momentum (of stepRT) defining its normal.
	// To define the state at "s", we need to find the point
	// where the track we're swimming intersects this plane.
	// We do this using a similar trick to what is done in
	// GetDistToRT(...) above. Namely, define a point on the helical
	// track segment in terms of the phi angle of the helix
	// such that phi=0 corresponds to track step itself. Call this point
	// h. 
	//
	//  h = Ro*(cos(phi) -1)*xdir  ~= -Ro*phi^2/2 *xdir
	//     + Ro*sin(phi)*ydir      ~=    Ro*phi   *ydir
	//     + dz/dphi*phi*zdir       = dz/dphi*phi *zdir
	//     + hpos
	//
	// where the xdir, ydir, zdir, hpos vectors define the
	// coordinate system at the current, "closest" step.
	//
	// Since the momentum vector at stepRT is normal to
	// this plane, the dot product  of it and the vector
	// pointing from stepRT to h will be zero. i.e.:
	//
	//   p.(h-pos) = 0
	//
	// where the period(.) means dot product, and pos is from
	// stepRT. Using the above definition of h, we get:
	//
	// p.(h-pos) = 0 = -(Ro/2)*(p.xdir)*phi^2
	//                + (Ro*p.ydir + dz/dphi*p.zdir)*phi
	//                + p.(hpos-pos)
	//
	// This formula is quadratic in phi (yea!)
	const DVector3 &p = stepRT->mom;
	double a = -(step.Ro/2.0)*p.Dot(step.xdir);
	double b = step.Ro*p.Dot(step.ydir) + dz_dphi*p.Dot(step.zdir);
	double c = p.Dot(step.pos - stepRT->pos);
	double phi1 = (-b + sqrt(b*b - 4.0*a*c))/(2.0*a);
	double phi2 = (-b - sqrt(b*b - 4.0*a*c))/(2.0*a);
	
	// Since the answer is not positive-definite, then we
	// have to choose which root to use. Choose the one
	// closest to zero (assuming both are finite)
	double phi = !finite(phi2) ? phi1:!finite(phi1) ? phi2:fabs(phi1)<fabs(phi2) ? phi1:phi2;

	// Now, use phi to calculate the point's coordinates
	// in the stepRT coordinate system at "s".
	// To get the momentum vector at "s", we need to
	// swim the particle from the "closest" swim
	// step to the actual point. Note that phi is defined
	// in the "step" coordinate system which is taken from
	// a backwards swimming negative particle.
	double dz = dz_dphi*phi;
	double Rodphi = step.Ro*phi;
	double ds = sqrt(dz*dz + Rodphi*Rodphi);
	if(phi<0.0)ds = -ds;
	stepper.SetStartingParams(-q, &step.pos, &step.mom);
	stepper.SetStepSize(ds);
	stepper.Step(&step.pos);
	stepper.GetMomentum(step.mom);

	// Flip momentum to get it in the "forward" direction for the state
	step.mom = -step.mom;
	DVector3 pos_diff = step.pos - stepRT->pos;
	state[state_px][0] = step.mom.Dot(stepRT->xdir);
	state[state_py][0] = step.mom.Dot(stepRT->ydir);
	state[state_pz][0] = step.mom.Dot(stepRT->zdir);
	state[state_x ][0] = pos_diff.Dot(stepRT->xdir);
	state[state_y ][0] = pos_diff.Dot(stepRT->ydir);

	// Return the distance from the wire
	step.Ro = stepper.GetRo();

	return GetDistToRT(wire, &step, step.s);
}

//------------------
// GetDistToRT
//------------------
DVector3 DTrack_factory::GetDistToRT(DVector3 &hit, swim_step_t *s2)
{
	/// Calculate the distance of the given hit vector(in the lab
	/// reference frame) to the Reference Trajectory which the
	/// given swim step belongs to. This uses the momentum directions
	/// and positions of the swim_step and it's nearest neighbors
	/// to define a curve and calculate the distance of the hit
	/// from it. The swim step should be the closest one to the hit.
	/// It should also be part of an array of swim steps so that there
	/// both one before and one after it in memory.
	
	// Get the distances of the swim step's nearest neighbors
	// on either side.
	swim_step_t *s1=s2, *s3=s2;
	s1--;	// step before
	s3++;	// step after
	DVector3 delta_s1 = s1->pos - hit;
	DVector3 delta_s3 = s3->pos - hit;
	double delta2_s1 = delta_s1.Mag2();
	double delta2_s3 = delta_s3.Mag2();

	// s_nn is pointer to nearest neighbor of s2 that is closest to hit
	swim_step_t *s_nn = delta2_s1<delta2_s3 ? s1:s3;

	// The point of closest approach on the R.T. for the given hit
	// should now be between the swim steps s1 and s_nn. Each swim
	// step has a momentum vector which is pointing in the
	// same direction as the tangent of the R.T. at that step.
	//
	// We split this up into two 2-D problems: Defining the 
	// parameters of a parabola in the x/y plane and separately
	// in the y/z plane. Note that the values of the hit vector
	// and of the steps are all given in the lab frame. However,
	// here we need to work in the R.T. frame of s2.
	//
	// We assume smooth functions of the form
	//
	// x = Q*y^2 + R*y + S
	//
	// z = T*y^2 + U*y +V
	//
	// Note that since the curve must go through the swim step at
	// s2, the values of S  and V are zero by definition.
	// Also, because the slope of the curve in the x/y plane
	// at the origin is zero, the value of R is also zero.
	//
	// Thus, the only parameter we need to find to define the shape
	// of the R.T. in the X/Y plane is "Q". This is given by:
	//
	// Q = (k1 - k2)/(2*(y1 - y2)) = k1/(2 * y1)
	//
	// where k1,k2 are the slopes of the tangents in the X/Y plane
	// at steps s_nn and s2 respectively. y1 and y2 are the "y"
	// positions of the  steps in the s2 R.T. frame. Note that
	// use of the facts that y2=0, k2=0 were used.
	//
	// Therefore we need to find y1 and k1 in the s2 R.T. frame.
	// k1 comes just from the momentum at s_nn. y1 is just the
	// position at s_nn. Both though are calculated in the s2 R.T.
	// frame.
	//
	// The value of U is given by the slope of the trajectory in
	// the y/z plane at s2. Thus comes just from the momentum
	// at s2.
	//
	// Start by finding the projection of the s_nn momentum
	// vector onto the x/y plane.
	double delta_x = s_nn->mom.Dot(s2->xdir);
	double delta_y = s_nn->mom.Dot(s2->ydir);
	double delta_z = s_nn->mom.Dot(s2->zdir);
	double k1_xy = delta_x/delta_y;
	double k1_zy = delta_z/delta_y;

	// Now find the projection of the s2 momentum vector
	// on the y/z plane to get "e"
	delta_y = s2->mom.Dot(s2->ydir);
	delta_z = s2->mom.Dot(s2->zdir);
	double U = delta_z/delta_y;

	// Find y1 in the R.T. frame at s2
	DVector3 tmpv = s_nn->pos - s2->pos; // shift so origin is that of R.T. at s2
	double y1 = tmpv.Dot(s2->ydir);
	
	double Q = k1_xy/2.0/y1;
	double T = (k1_zy - U)/2.0/y1;

	// Now we need to find the coordinates of the point on the
	// RT closest to the hit. We do this by taking the derivative
	// of the distance squared, setting it equal to zero, and
	// solving for y. This leads to a 3rd order polynomial:
	//
	// Ay^3 + By^2 + Cy + D = 0
	//
	// where:
	//     A = 2(Q^2 + T^2)
	//     B = 3TU
	//     C = U^2 - 2Qx - 2Tz +1
	//     D = -(y + Uz)
	//     x,y,z = coordinates of hit in s2 RT frame
	//
	
	// Transform "hit" from lab coordinates to R.T. coordinates at s2
	tmpv = hit - s2->pos;
	DVector3 hit_rt(tmpv.Dot(s2->xdir), tmpv.Dot(s2->ydir), tmpv.Dot(s2->zdir));
	double xh = hit_rt.x();
	double yh = hit_rt.y();
	double zh = hit_rt.z();
	double A = 2.0*(Q*Q + T*T);
	double B = 3.0*T*U;
	double C = U*U - 2.0*Q*xh - 2.0*T*zh +1.0;
	double D = -(yh + U*zh);
	
	// OK, at this point there is a potential problem. If the trajectory's
	// momentum is in essentially the same direction for both s2 and s_nn,
	// then the above problem is really linear and Q and T are both
	// essentially zero. This means both A and B are zero above and
	// the solution for y is just simply y=-D/C. If the momenta are NOT
	// in the same direction, then we must solve the 3rd order polynomial.
	// The problem is that the first step in solving a 3rd order poly
	// is to divide by the cubic term's coefficient so that you get a new
	// equation whos cubic coeeficient is one. In the case that the
	// cubic coefficient is at or near zero, this blows up the the
	// other coefficients leading to an incorrect result. Thus, we
	// need to decide here whether to approximate the trajectory between
	// s2 and s_nn as a line or a parabola. We do this by checking the
	// angle between the two momenta.
	double y;
	if(s2->mom.Angle(s_nn->mom) <= 0.03){ // 0.03rad = 1.7 degrees
		// momentum doesn't really change between these two points.
		// use linear approximation
		y = -D/C;
	}else{
	
		// Now we have to solve the 3rd order polynomial for the y value of
		// the point of closest approach on the RT. This is a well documented
		// procedure. Essentially, when you have an equation of the form:
		//
		//  x^3 + a2*x^2 + a1*x + a0 = 0;
		//
		// a change of variables is made such that w = x + a2/3 which leads
		// to a third order poly with no w^2 term:
		//
		//  w^3 + 3.0*b*w + 2*c = 0 
		//
		// where:
		//    b = a1/3 - (a2^2)/9
		//    c = a0/2 - a1*a2/6  + (a2^3)/27
		//
		// The one real root of this is:
		//
		//  w0 = q - p
		//
		// where:
		//    q^3 = d - c
		//    p^3 = d + c
		//    d^2 = b^3 + c^2
	
		double a0 = D/A;
		double a1 = C/A;
		double a2 = B/A;

		double b = a1/3.0 - a2*a2/9.0;
		double c = a0/2.0 - a1*a2/6.0 + a2*a2*a2/27.0;
	
		double d = sqrt(pow(b, 3.0) + pow(c, 2.0));
		double q = pow(d - c, 1.0/3.0);
		double p = pow(d + c, 1.0/3.0);
	
		double w0 = q - p;
		y = w0 - a2/3.0;
	}
		
	// Calculate the y and z coordinates of the RT DOCA point 
	// from the x coordinate
	double x = Q*y*y;
	double z = T*y*y + U*y;
	// Return vector pointing from DOCA point to hit in RT frame
	DVector3 doca_point(x,y,z);

	return hit_rt - doca_point;
}

#define USE_MINUIT 0

#if USE_MINUIT
// This is for MINUIT which is really not suited for a multi-threaded
// event processing environment. We do a few backflips here ...
static void FCN(int &npar, double *derivatives, double &chisq, double *par, int iflag);

static pthread_mutex_t trk_mutex = PTHREAD_MUTEX_INITIALIZER;
static JEventLoop *g_loop = NULL;
static const DMagneticFieldMap *g_bfield = NULL;
static JFactory<DTrackHit> *g_fac_trackhit;
static const DTrackCandidate *g_trackcandidate;
static double min_chisq, min_par[6];
static TMinuit *minuit=NULL;

typedef struct{
	double x,y,z;
	double dxdz, dydz;
	double d2xdz2, d2ydz2;
}trk_step_t;

#endif // USE_MINUIT

//------------------
// init
//------------------
jerror_t DTrack_factory::init(void)
{
#if USE_MINUIT
	pthread_mutex_lock(&trk_mutex);
	if(minuit==NULL){
		int icondn;
		minuit = new TMinuit();
		minuit->mninit(0,0,0);
		minuit->SetFCN(FCN);
		minuit->mncomd("CLEAR", icondn);
		minuit->mncomd("SET PRI -1", icondn);
		minuit->mncomd("SET NOWARNINGS", icondn);
		minuit->mncomd("SET BATCH", icondn);
		minuit->mncomd("SET STRATEGY 2", icondn);
		
		minuit->mncomd("SET PRI -1", icondn);
		minuit->mncomd("SET NOWARNINGS", icondn);
		minuit->mncomd("SET BATCH", icondn);
		minuit->mncomd("SET STRATEGY 2", icondn);

		minuit->DefineParameter(0,"p",		0.0, 0.02, 0.0, 0.0);
		minuit->DefineParameter(1,"theta",	0.0, 0.17, 0.0, 0.0);
		minuit->DefineParameter(2,"phi",		0.0, 0.17, 0.0, 0.0);
		minuit->DefineParameter(3,"x",		0.0, 0.50, 0.0, 0.0);
		minuit->DefineParameter(4,"y",		0.0, 0.50, 0.0, 0.0);
		minuit->DefineParameter(5,"z",		0.0, 7.0, 0.0, 0.0);

		minuit->mncomd("FIX 4", icondn);
		minuit->mncomd("FIX 5", icondn);
		minuit->mncomd("FIX 6", icondn);
	}
	pthread_mutex_unlock(&trk_mutex);

#endif // USE_MINUIT

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrack_factory::brun(JEventLoop *loop, int runnumber)
{
	//dgeom = loop->GetJApplication()->GetJGeometry(runnumber);
	//bfield = NULL;
	//bfield = dgeom->GetDMagneticFieldMap();
	
	bfield = new DMagneticFieldMap(); // temporary until new geometry scheme is worked out
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrack_factory::evnt(JEventLoop *loop, int eventnumber)
{
	// Get the track candidates and hits
	vector<const DTrackCandidate*> trackcandidates;
	vector<const DTrackHit*> trackhits;
	loop->Get(trackcandidates);
#if USE_MINUIT
	JFactory<DTrackHit> *fac_trackhit =
#endif //USE_MINUIT
	loop->Get(trackhits,"MC");
		
	for(unsigned int i=0; i<trackcandidates.size(); i++){
		const DTrackCandidate *trackcandidate = trackcandidates[i];
		DTrack *track = new DTrack;
		
		// Copy in starting values
		track->q = trackcandidate->q;
		track->p = trackcandidate->p;
		track->theta = trackcandidate->theta;
		track->phi = trackcandidate->phi;
		track->x = 0.0;
		track->y = 0.0;
		track->z = trackcandidate->z_vertex;
		track->candidateid = trackcandidate->id;
		
#if USE_MINUIT
		// Lock mutex so only one thread at a time accesses minuit
		pthread_mutex_lock(&trk_mutex);
		g_loop = loop;
		g_bfield = bfield;
		g_fac_trackhit = fac_trackhit;
		g_trackcandidate = trackcandidate;
		
		// Calculate chi-square value
		double par[6];
		par[0] = track->p;
		par[1] = track->theta;
		par[2] = track->phi;
		par[3] = track->x;
		par[4] = track->y;
		par[5] = track->z;
		int npar=6;
		double dfdpar[6];
		double chisq=0.0;;
		int iflag = 2;
		FCN(npar, dfdpar, chisq, par, iflag);
		
		// Set the initial parameter values
		for(int i=0;i<6;i++){
			int icondn;
			char cmd[256];
			sprintf(cmd, "SET PAR %d %f", i+1, par[i]);
			min_par[i] = par[i];
			minuit->mncomd(cmd, icondn);
		}		
		
		min_chisq = 1000.0;
		//minuit->Migrad();
		//minuit->mncomd("MINIMIZE 100", icondn);
		//int icondn;
		minuit->mncomd("SEEK 100 1", icondn);

//cout<<__FILE__<<":"<<__LINE__<<" initial:"<<chisq<<"  final:"<<min_chisq<<endl;
//for(int i=0; i<6; i++){
//	cout<<__FILE__<<":"<<__LINE__<<" -- p"<<i<<"  initial:"<<par[i]<<"  final:"<<min_par[i]<<endl;
//}

		// Unlock mutex so only one thread at a time accesses minuit
		pthread_mutex_unlock(&trk_mutex);
		
		if(min_chisq < chisq){
			track->p = min_par[0];
			track->theta = min_par[1];
			track->phi = min_par[2];
			track->x = min_par[3];
			track->y = min_par[4];
			track->z = min_par[5];
		}

#endif // USE_MINUIT

		_data.push_back(track);
	}

	return NOERROR;
}


#if USE_MINUIT
//------------------
// FCN   -- for MINUIT
//------------------
void FCN(int &npar, double *derivatives, double &chisq, double *par, int iflag)
{
	// Step a singley charged particle with the momentum and vertex
	// position given by par through the magnetic field, calculating
	// the chisq from the closest hits. We need to keep track of the
	// positions and first and second derivatives of each step along
	// the trajectory. We will first step through field from vertex 
	// until we either hit the approximate BCAL, UPV, or FCAL positions.

	// Minuit will sometimes try zero momentum tracks. Immediately return
	// a large ChiSq when that happens
	if(par[0] < 0.050){chisq=1000.0; return;}
	
	vector<double> chisqv;
	DVector3 p;
	p.SetMagThetaPhi(par[0],par[1],par[2]);
	DVector3 pos(par[3],par[4],par[5]);
	DMagneticFieldStepper *stepper = new DMagneticFieldStepper(g_bfield, g_trackcandidate->q, &pos, &p);
	stepper->SetStepSize(0.1);
	
	// Step until we hit a boundary
	vector<trk_step_t> trk_steps;
	for(int i=0; i<10000; i++){
		stepper->Step(&pos);
		trk_step_t trk_step;
		trk_step.x = pos.X();
		trk_step.y = pos.Y();
		trk_step.z = pos.Z();
		trk_steps.push_back(trk_step);
		
		//if(pos.Perp()>65.0){cout<<__FILE__<<":"<<__LINE__<<" hit BCAL"<<endl;break;} // ran into BCAL
		//if(pos.Z()>650.0){cout<<__FILE__<<":"<<__LINE__<<" hit FCAL"<<endl;break;} // ran into FCAL
		//if(pos.Z()<-50.0){cout<<__FILE__<<":"<<__LINE__<<" hit UPV"<<endl;break;} // ran into UPV
		if(pos.Perp()>65.0){break;} // ran into BCAL
		if(pos.Z()>650.0){break;} // ran into FCAL
		if(pos.Z()<-50.0){break;} // ran into UPV
	}
	delete stepper;

#if 0
	// Calculate first derivatives in trk_steps
	for(unsigned int i=1; i<trk_steps.size()-1; i++){
		// We do this by averaging the slopes of the lines connecting
		// this point to the previous and nexxt ones. We assume
		// that the stepper used small enough steps to vary the
		// positions smoothly.
		trk_step_t *prev = &trk_steps[i-1];
		trk_step_t *curr = &trk_steps[i];
		trk_step_t *next = &trk_steps[i+1];
		curr->dxdz = ((curr->x-prev->x)/(curr->z-prev->z) + (next->x-curr->x)/(next->z-curr->z))/2.0;
		curr->dydz = ((curr->y-prev->y)/(curr->z-prev->z) + (next->y-curr->y)/(next->z-curr->z))/2.0;
	}

	// Calculate second derivatives in trk_steps
	for(unsigned int i=2; i<trk_steps.size()-2; i++){
		// We do this by averaging the slopes of the lines connecting
		// this point to the previous and nexxt ones. We assume
		// that the stepper used small enough steps to vary the
		// positions smoothly.
		trk_step_t *prev = &trk_steps[i-1];
		trk_step_t *curr = &trk_steps[i];
		trk_step_t *next = &trk_steps[i+1];
		curr->d2xdz2 = ((curr->dxdz-prev->dxdz)/(curr->z-prev->z) + (next->dxdz-curr->dxdz)/(next->z-curr->z))/2.0;
		curr->d2ydz2 = ((curr->dydz-prev->dydz)/(curr->z-prev->z) + (next->dydz-curr->dydz)/(next->z-curr->z))/2.0;
	}
#endif //0 
	
	// Loop over the track hits for this candidate using
	// the closest trk_step to determine the distance
	// from the hit to the track.
	vector<const DTrackHit*> trackhits;
	vector<void*>& vptrs = g_fac_trackhit->Get();
	for(unsigned int i=0; i<vptrs.size(); i++){
		const DTrackHit *trackhit = (const DTrackHit*)vptrs[i];
		if(!trackhit)continue;
		
		// Ignore hits outside of the CDC and FDC for now
		if(trackhit->z <= 0.0)continue;
		if(trackhit->z > 405.0)continue;
		double r = sqrt(pow((double)trackhit->x,2.0) + pow((double)trackhit->y,2.0));
		if(r >= 65.0)continue;
		
		// We don't really want the hit closest physically in space.
		// Rather, we want the distance measured in standard deviations
		// of detector resolutions. This depends upon the detector type
		// from which the hit came. Note that this also depends on the
		// step size used in the stepper. The FDC in particular needs
		// much larger sigma values due to detector resolution alone.
		double sigmaxy;
		double sigmaz;
		switch(trackhit->system){
			case SYS_CDC:	sigmaxy = 0.80;	sigmaz = 7.50;	break;
			case SYS_FDC:	sigmaxy = 0.50;	sigmaz = 0.50;	break;
			default:			sigmaxy = 10.00;	sigmaz = 10.00;
		}

		double xh = trackhit->x;
		double yh = trackhit->y;
		double zh = trackhit->z;

		// Loop over all steps to find the closest one to this hit
		double r2_min = 1.0E6;
		trk_step_t *trk_step_min=NULL;
		for(unsigned int j=0; j<trk_steps.size(); j++){
			trk_step_t *trk_step = &trk_steps[j];
			double xs = trk_step->x;
			double ys = trk_step->y;
			double zs = trk_step->z;
			
			double r2 = pow((xh-xs)/sigmaxy,2.0) + pow((yh-ys)/sigmaxy,2.0) + pow((zh-zs)/sigmaz,2.0);
//cout<<__FILE__<<":"<<__LINE__<<" dx="<<xh-xs<<" dy="<<yh-ys<<" dz="<<zh-zs<<"   z="<<zs<<endl;
			if(r2<r2_min){
				r2_min = r2;
				trk_step_min = trk_step;
			}
		}
//cout<<__FILE__<<":"<<__LINE__<<" ##### r_min= "<<sqrt(r2_min)<<"  z="<<trk_step_min->z<<endl;
		
		// Now calculate the distance from the track to the hit
		double d = sqrt(r2_min); // use the distance to the step for now
		
		// If the hit and the step are more than 3 sigma apart, then
		// don't include this hit in the ChiSq.
		if(d>3.0)continue;

		// Add this hit to the chisq vector
		chisqv.push_back(d*d);
	}

	// Sum up total chisq/dof for this event
	chisq = 0.0;
	for(unsigned int i=0; i<chisqv.size(); i++)chisq += chisqv[i];
	chisq /= (double)chisqv.size() - 4.0;

//chisq = pow((par[0]-1.0)/0.100,2.0) + pow((par[5]-65.0)/7.0,2.0);

	// If NO hits were close to this track (how does that even happen?)
	// then the size of chisqv will be zero and chisq will now be
	// "nan". In this case, set the chisq to a very large value to
	// indicate a bad fit.
	if(chisqv.size() < 5)chisq=1000.0;

	if(chisq<min_chisq){
		min_chisq = chisq;
		for(int i=0;i<6;i++)min_par[i] = par[i];
	}

//cout<<__FILE__<<":"<<__LINE__<<" chisq= "<<chisq<<"  chisqv.size()="<<chisqv.size()<<endl;
}
#endif // USE_MINUIT

#endif // 0
