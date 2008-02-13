// $Id$
//
//    File: DTrack_factory_ALT1.cc
// Created: Fri. Feb. 8, 2008
// Creator: davidl
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
#include "DTrack_factory_ALT1.h"
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


//------------------
// DTrack_factory_ALT1   (Constructor)
//------------------
DTrack_factory_ALT1::DTrack_factory_ALT1()
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
	MAX_CHISQ_DIFF = 1.0E-2;
	MAX_FIT_ITERATIONS = 50;
	SIGMA_CDC = 0.0200;
	SIGMA_FDC_ANODE = 0.0200;
	SIGMA_FDC_CATHODE = 0.0200;
	CHISQ_MAX_RESI_SIGMAS = 100.0;
	CHISQ_GOOD_LIMIT = 2.0;
	LEAST_SQUARES_DP = 0.0001;
	LEAST_SQUARES_DX = 0.010;
	LEAST_SQUARES_MIN_HITS = 3;
	LEAST_SQUARES_MAX_E2NORM = 1.0E6;
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
	gPARMS->SetDefaultParameter("TRKFIT:DEFAULT_STEP_SIZE",		DEFAULT_STEP_SIZE);
	gPARMS->SetDefaultParameter("TRKFIT:MIN_CDC_HIT_PROB",			MIN_CDC_HIT_PROB);
	gPARMS->SetDefaultParameter("TRKFIT:MAX_CDC_DOUBLE_HIT_PROB",	MAX_CDC_DOUBLE_HIT_PROB);
	gPARMS->SetDefaultParameter("TRKFIT:MIN_FDC_HIT_PROB",			MIN_FDC_HIT_PROB);
	gPARMS->SetDefaultParameter("TRKFIT:MAX_FDC_DOUBLE_HIT_PROB",	MAX_FDC_DOUBLE_HIT_PROB);
	gPARMS->SetDefaultParameter("TRKFIT:TOF_MASS",					TOF_MASS);

}

//------------------
// DTrack_factory_ALT1   (Destructor)
//------------------
DTrack_factory_ALT1::~DTrack_factory_ALT1()
{
	if(tmprt)delete tmprt;
	for(unsigned int i=0; i<rtv.size(); i++)delete rtv[i];
}

//------------------
// init
//------------------
jerror_t DTrack_factory_ALT1::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrack_factory_ALT1::brun(JEventLoop *loop, int runnumber)
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
		cdc_can_resi = (TH1F*)gROOT->FindObject("cdc_can_resi");
		fdc_can_resi = (TH1F*)gROOT->FindObject("fdc_can_resi");
		fdc_can_resi_cath = (TH1F*)gROOT->FindObject("fdc_can_resi_cath");

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
		if(!Npasses)Npasses = new TH1F("Npasses","Npasses", 101, -0.5, 100.5);
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
		if(!cdc_can_resi)cdc_can_resi = new TH1F("cdc_can_resi","Residual of CDC hits with candidate tracks", 200, -1.0, 1.0);
		if(!fdc_can_resi)fdc_can_resi = new TH1F("fdc_can_resi","Residual of FDC hits with candidate tracks", 200, -1.0, 1.0);
		if(!fdc_can_resi_cath)fdc_can_resi_cath = new TH1F("fdc_can_resi_cath","Residual of FDC cathode hits with candidate tracks", 200, -1.0, 1.0);

		dapp->Unlock();		
	}

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrack_factory_ALT1::evnt(JEventLoop *loop, int eventnumber)
{
	// Store current event number
	this->eventnumber = eventnumber;

	// Get input data
	trackcandidates.clear();
	cdctrackhits.clear();
	fdctrackhits.clear();
	loop->Get(trackcandidates);
	loop->Get(cdctrackhits);
	loop->Get(fdctrackhits);

	// Calculate the probablity of each hit belonging to each candidate
	FindHitCandidateProbabilities();
	
	// Consistency check
	if(cdctrackhits.size()!=cdcprobs.size() || fdctrackhits.size()!=fdcprobs.size()){
		_DBG_<<"Hit vector and assignment vector not the same size!"<<endl;
		_DBG_<<"  cdctrackhits.size()="<<cdctrackhits.size()<<" ; cdcprobs.size()"<<cdcprobs.size()<<endl;
		_DBG_<<"  fdctrackhits.size()="<<fdctrackhits.size()<<" ; fdcprobs.size()"<<fdcprobs.size()<<endl;
		return VALUE_OUT_OF_RANGE;
	}

	// Loop over track candidates
	for(unsigned int i=0; i<trackcandidates.size(); i++){
	
		if(debug_level>1){
			_DBG__;
			_DBG_<<"============ Fitting Candidate "<<i<<" ================"<<endl;
		}

		// Copy the hits consistent with this track into the
		// cdchits_on_track and fdchits_on_track vectors
		cdchits_on_track.clear();
		for(unsigned int j=0; j<cdcprobs.size(); j++){
			if(cdcprobs[j][i]>=MIN_CDC_HIT_PROB)cdchits_on_track.push_back(cdctrackhits[j]);
		}
		fdchits_on_track.clear();
		for(unsigned int j=0; j<fdcprobs.size(); j++){
			if(fdcprobs[j][i]>=MIN_FDC_HIT_PROB)fdchits_on_track.push_back(fdctrackhits[j]);
		}

		// Fit the track
		DTrack *track = FitTrack(rtv[i], trackcandidates[i]->id);

		// If fit is successful, then store the track
		if(track)_data.push_back(track);
	}

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTrack_factory_ALT1::fini(void)
{
	return NOERROR;
}

//------------------
// FitTrack
//------------------
DTrack* DTrack_factory_ALT1::FitTrack(DReferenceTrajectory* rt, int candidateid)
{
	/// Fit a track candidate
	
	// Debug message
	if(debug_level>2){
		_DBG__;
		_DBG_<<"cdchits_on_track.size="<<cdchits_on_track.size()<<"  fdchits_on_track.size="<<fdchits_on_track.size()<<endl;
		_DBG_<<"cdctrackhits.size="<<cdctrackhits.size()<<"  fdctrackhits.size="<<fdctrackhits.size()<<endl;
	}
	
	// Do initial hit-based fit
	if(debug_level>3)_DBG_<<"============= Hit-based ==============="<<endl;
	double chisq_hit_based;
	for(int i=0; i<1; i++){
		if(LeastSquaresB(kHitBased, rt, chisq_hit_based) != FIT_OK){
			if(debug_level>3)_DBG_<<"Hit based fit failed for candidate "<<candidateid<<endl;
			return NULL;
		}
	}
	
	// Iterate over time-based fits until it converges or we reach the 
	// maximum number of iterations
	if(debug_level>3)_DBG_<<"============= Time-based ==============="<<endl;
	double chisq_time_based;
	double last_chisq = 1.0E6;
	int Niterations;
	for(Niterations=0; Niterations<MAX_FIT_ITERATIONS; Niterations++){
		if(LeastSquaresB(kTimeBased, rt, chisq_time_based) != FIT_OK){
			if(debug_level>3)_DBG_<<"Time based fit failed for candidate "<<candidateid<<" on iteration "<<Niterations<<endl;
			return NULL;
		}

		// Fill debug histos
		if(DEBUG_HISTS){
			chisq_vs_pass->Fill(Niterations+1, chisq_time_based);
			dchisq_vs_pass->Fill(Niterations+1, last_chisq-chisq_time_based);			
		}

		// If the chisq is too large, then consider it a hopeless cause
		if(chisq_time_based>1.0E4){
			if(debug_level>3)_DBG_<<"Time based fit chisq too large for candidate "<<candidateid<<" on iteration "<<Niterations<<endl;
			return NULL;
		}
		
		// Check if we converged
		if(fabs(last_chisq-chisq_time_based) < MAX_CHISQ_DIFF)break;
		last_chisq = chisq_time_based;
	}

#if 0
	// Get starting position and momentum from candidate reference trajectory
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
			if(debug_level>1)_DBG_<<"Fit failed for iteration "<<Niterations<<" (not necessarily fatal)"<<endl;
			break;
		}

		if(Niterations==0)initial_chisq = chisq;
		if(debug_level>3)_DBG_<<" ---- iteration "<<Niterations<<"  chisq="<<chisq<<endl;
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
#endif

	if(debug_level>1)_DBG_<<" Niterations="<<Niterations<<"  chisq_time_based="<<chisq_time_based<<" p="<<rt->swim_steps[0].mom.Mag()<<endl;
	
	// At this point we must decide whether the fit succeeded or not.
	// We'll consider the fit a success if any of the following is true:
	//
	// 1. We got through at least one iteration in the above loop
	// 2. The chi-sq is less than CHISQ_GOOD_LIMIT
	// 3. MAX_FIT_ITERATIONS is zero (for debugging)
	bool fit_succeeded = false;
	if(Niterations>0)fit_succeeded = true;
	if(chisq_time_based<=CHISQ_GOOD_LIMIT)fit_succeeded = true;
	if(MAX_FIT_ITERATIONS==0)fit_succeeded = true;
	if(!fit_succeeded)return NULL;
	
	// Find point of closest approach to target and use parameters
	// there for vertex position and momentum
	//double s;
	//rt->DistToRT(target, &s);
	//rt->GetLastDOCAPoint(vertex_pos, vertex_mom);
	//_DBG_<<"final chi-sq:"<<chisq<<"   DOF="<<chisqv.size()<<"   Niterations="<<Niterations<<" dp_over_p="<<dp_over_p<<endl;

	DVector3 &vertex_mom = rt->swim_steps[0].mom;
	DVector3 &vertex_pos = rt->swim_steps[0].origin;

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
	track->chisq	= chisq_time_based;
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
		initial_chisq_vs_Npasses->Fill(Niterations, chisq_hit_based);
		FillDebugHists(rt, vertex_pos, vertex_mom);
	}

	// Debugging messages
	if(debug_level>2)_DBG_<<" -- Fit succeeded: q="<<track->q<<" p="<<track->p<<" theta="<<track->theta<<" phi="<<track->phi<<endl;

	return track;
}

//------------------
// FindHitCandidateProbabilities
//------------------
void DTrack_factory_ALT1::FindHitCandidateProbabilities(void)
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
	
	// Increase capcity (if needed) and clear out hit-track probability matrices.
	if(cdcprobs.size()!=cdctrackhits.size())cdcprobs.resize(cdctrackhits.size());
	if(fdcprobs.size()!=fdctrackhits.size())fdcprobs.resize(fdctrackhits.size());
	for(unsigned int i=0; i<cdcprobs.size(); i++)cdcprobs[i].clear();
	for(unsigned int i=0; i<fdcprobs.size(); i++)fdcprobs[i].clear();
	
	// Loop over track candidates
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
		
		// Debugging messages
		if(debug_level>10){
			_DBG_<<"Candidate "<<i<<endl;
			for(unsigned int j=0; j<cdcprob.size(); j++)_DBG_<<"  cdcprob["<<j<<"] = "<<cdcprob[j]<<endl;
			for(unsigned int j=0; j<fdcprob.size(); j++)_DBG_<<"  fdcprob["<<j<<"] = "<<fdcprob[j]<<endl;
		}
	}
}

//------------------
// GetCDCTrackHitProbabilities
//------------------
void DTrack_factory_ALT1::GetCDCTrackHitProbabilities(DReferenceTrajectory *rt, vector<double> &prob)
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
	//double sigma = sqrt(pow(SIGMA_CDC,2.0) + pow(0.4000,2.0));
	double sigma = 0.261694; // empirically from single pi+ events

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
		double probability = finite(resi) ? exp(-pow(resi/sigma,2.0)):0.0;
		prob.push_back(probability);

		if(DEBUG_HISTS){
			cdc_can_resi->Fill(resi);
			cdc_single_hit_prob->Fill(probability);
		}
	}
}

//------------------
// GetFDCTrackHitProbabilities
//------------------
void DTrack_factory_ALT1::GetFDCTrackHitProbabilities(DReferenceTrajectory *rt, vector<double> &prob)
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
	double sigma_anode = sqrt(pow(SIGMA_FDC_ANODE,2.0) + pow(1.000,2.0));
	double sigma_cathode = sqrt(pow(SIGMA_FDC_CATHODE,2.0) + pow(1.000,2.0));
	
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
		if(DEBUG_HISTS)cdc_can_resi->Fill(resi);
		
		// Use an un-normalized gaussian so that for a residual
		// of zero, we get a probability of 1.0.
		double p = finite(resi) ? exp(-pow(resi/sigma_anode,2.0)):0.0;

		// Cathode
		double u=0;
		double resic=0;
		if(USE_FDC_CATHODE){
			u = rt->GetLastDistAlongWire();
			resic = u - hit->s;
			if(DEBUG_HISTS)fdc_can_resi_cath->Fill(resic);

			// Same as for the anode. We multiply the
			// probabilities to get a total probability
			// based on both the anode and cathode hits.
			p *= finite(resic) ? exp(-pow(resic/sigma_cathode,2.0)):0.0;
		}
		
		if(debug_level>10)_DBG_<<" fdc hit "<<j<<": dist,doca="<<dist<<", "<<doca<<"  u,s="<<u<<", "<<hit->s<<" resi="<<resi<<" resic="<<resic<<endl;
		
		prob.push_back(p);
	}
}

//------------------
// ChiSq
//------------------
double DTrack_factory_ALT1::ChiSq(DMatrix &state, const swim_step_t *start_step, DReferenceTrajectory *rt, vector<const DCoordinateSystem*> &wires, vector<DVector3> &shifts, vector<double> &errs, vector<double> &chisqv)
{
	/// Calculate the chi-squared for a track specified by state relative to the
	/// given reference trajectory. This is just a wrapper for 
	/// ChiSq(DReferenceTrajectory *rt, vector<const DCoordinateSystem*> &wires, vector<DVector3> &shifts, vector<double> &errs, vector<double> &chisqv)
	/// that accepts the state vector and re-swims the trajectory.
	
	// "v" direction is perpendicular to both the rt direction and the
	// x-direction. See LeastSquares() for more.
	DVector3 vdir = start_step->sdir.Cross(start_step->mom);
	vdir.SetMag(1.0);

	DVector3 pos =   start_step->origin
						+ state[state_x ][0]*start_step->sdir
						+ state[state_v ][0]*vdir;
	DVector3 mom =   state[state_px][0]*start_step->sdir
						+ state[state_py][0]*start_step->tdir
						+ state[state_pz][0]*start_step->udir;

	if(!rt){
		_DBG_<<"NULL pointer passed for DReferenceTrajectory object!"<<endl;
		return 1.0E6;
	}
	
	rt->Swim(pos,mom);

	return ChiSq(rt, wires, shifts, errs, chisqv);
}

//------------------
// ChiSq
//------------------
double DTrack_factory_ALT1::ChiSq(DReferenceTrajectory *rt, vector<const DCoordinateSystem*> &wires, vector<DVector3> &shifts, vector<double> &errs, vector<double> &chisqv)
{
	/// Calculate the chisq for a track represented by the given reference
	/// trajectory with the given list of wires. The values in the "shifts"
	/// vector are used to shift the wire positions and can therefore be used
	/// to implement the drift time information. If no shifts are required, then
	/// the shifts vector should be empty. The values in the errs vector
	/// should be the errs for each of the wire measurements.
	///
	/// The return value is the chi-squared per dof. Upon return, the chisqv
	/// vector will be filled with the individual chisq contributions from
	/// each wire hit.
	///
	/// Upon entry, the lengths of the vectors "wires", "errs", and "shifts"
	/// (if non-zero) are checked to make sure they all match. If not, an error
	/// message is printed and a value of 1.0E6 is returned and chisqv
	/// will be empty
	
	// Make sure input vectors match
	chisqv.clear();
	bool use_shifts = shifts.size()==wires.size();
	if(wires.size()!=errs.size() || (shifts.size()!=0 && !use_shifts)){
		_DBG_<<"Error! the wires, errs, and shifts vectors are out of alignment."<<endl;
		_DBG_<<"wires.size()="<<wires.size()<<" errs.size()="<<errs.size()<<" shifts.size()="<<shifts.size()<<endl;
		return 1.0E6;
	}
	
	// Loop over wires
	double chisq = 0.0;
	int dof=0;
	for(unsigned int i=0; i<wires.size(); i++){
		DCoordinateSystem wire = *wires[i];
		if(use_shifts)wire.origin += shifts[i];
		
		// Get distance of the (shifted) wire from the reference trajectory
		double s;
		double d = rt->DistToRT(&wire, &s);
		
		// Add this to the chisq. It may turn out that d is NAN in which case we
		// don't want to include it in the total tally. We still must add it
		// however to the chisqv since that must have an entry for every entry
		// in the wires vector.
		//
		// The value going into chisqv needs to be a signed quantity with positive
		// values meaning the doca is larger than the dist and negative values
		// meaning the doca is smaller than the dist (doca from track, dist
		// from drift time). To decide the sign, we have to 
		double c = d/errs[i];
		chisqv.push_back(c);
		if(finite(c)){
			chisq += c*c;
			dof++;
			if(debug_level>10)_DBG_<<"  chisqv["<<dof<<"] = "<<c<<"  d="<<d<<"  Rwire="<<wires[i]->origin.Perp()<<"  Rshifted="<<wire.origin.Perp()<<endl;
		}
	}

	// If it turns out the dof is zero, return 1.0E6. Otherwise, return
	// the chisq/dof
	return dof<2 ? 1.0E6:chisq/(double)dof;
}

//------------------
// LeastSquaresB
//------------------
DTrack_factory_ALT1::fit_status_t DTrack_factory_ALT1::LeastSquaresB(fit_type_t fit_type, DReferenceTrajectory *rt, double &chisq)
{
	/// Fit the track with starting parameters given in the first step
	/// of the reference trajectory rt. On return, the reference
	/// trajectory rt will represent the final fit parameters and
	/// chisq will contain the chisq/dof of the track.
	///
	/// This determines the best fit of the track using the least squares method
	/// described by R. Mankel Rep. Prog. Phys. 67 (2004) 553-622 pg 565.
	/// Since it uses a linear approximation for the chisq dependance on
	/// the fit parameters, several calls may be required for convergence.
	
	// For fitting, we want to define a coordinate system very similar to the
	// Reference Trajectory coordinate system. The difference is that we want
	// the position to be defined in a plane perpendicular to the momentum.
	// The RT x-direction is in this plane, but the momentum itself lies
	// somewhere in the Y/Z plane so that neither Y nor Z makes a good choice
	// for the second postion dimension. We will call the second coordinate in 
	// the perpendicular plane "v" which is the coordinate along the axis that
	// is perpendicular to both the x-direction and the momentum direction.
	swim_step_t &start_step = rt->swim_steps[0];
	DVector3 pos = start_step.origin;
	DVector3 mom = start_step.mom;
	
	// Make list of wires
	vector<const DCoordinateSystem*> wires;
	vector<DVector3> shifts;
	vector<double> errs;

#if 1
	// --- Target ---
	// Define the target as though it were a wire so it is included
	// in the chi-sq
	DCoordinateSystem target;
	target.origin.SetXYZ(0.0, 0.0, pos.Z());
	target.sdir.SetXYZ(1.0, 0.0, 0.0);
	target.tdir.SetXYZ(0.0, 1.0, 0.0);
	target.udir.SetXYZ(0.0, 0.0, 1.0);
	target.L=3.0;
	wires.push_back(&target);
	errs.push_back(0.1);
	if(fit_type!=kHitBased)shifts.push_back(DVector3(0.0, 0.0, 0.0));
#endif

	// --- CDC ---
	for(unsigned int i=0; i<cdchits_on_track.size(); i++){
		const DCDCTrackHit *hit = cdchits_on_track[i];
		const DCoordinateSystem *wire = hit->wire;
		wires.push_back(wire);

		// Fill in shifts and errs vectors based on whether we're doing
		// hit-based or time-based tracking
		if(fit_type==kHitBased){
			// If we're doing hit-based tracking then only the wire positions
			// are used and the drift time info is ignored. Thus, we don't
			// have to calculate the shift vectors
			errs.push_back(0.261694); // empirically from single pi+ events
		}else{
			// Find the DOCA point for the RT and use the momentum direction
			// there and the wire direction to determine the "shift".
			// Note that whether the shift is to the left or to the right
			// is not decided here. That ambiguity is left to be resolved later
			// by applying a minus sign (or not) to the shift.
			DVector3 pos_doca, mom_doca;
			double s;
			rt->DistToRT(wire, &s);
			rt->GetLastDOCAPoint(pos_doca, mom_doca);
			DVector3 shift = wire->udir.Cross(mom_doca);
			
			// The magnitude of the shift is based on the drift time. The
			// value of the dist member of the DCDCTrackHit object does not
			// subtract out the TOF. This can add 50-100 microns to the
			// resolution in the CDC. Here, we actually can calculate the TOF
			// (for a given mass hypothesis).
			double mass = 0.13957;
			double beta = 1.0/sqrt(1.0 + pow(mass/mom_doca.Mag(), 2.0))*2.998E10;
			double tof = s/beta/1.0E-9; // in ns
			double dist = hit->dist*((hit->tdrift-tof)/hit->tdrift);
			shift.SetMag(dist);

			// If we're doing time-based tracking then there is a sign ambiguity
			// to each of the "shifts". It is not realistic to check every possible
			// solution so we have to make a guess as to what the sign is for each
			// based on the current reference trajectory. Presumably, if, we're
			// doing time-based track fitting then we have already gone through a pass
			// of hit-based track fitting so the initial reference trajectory should be
			// a pretty good indicator as to what side of the wire the track went
			// on. We use this in the assignments for now. In the future, we should
			// try multiple cases for hits close to the wire where the ambiguity 
			// is not clearly resolved.
			
			// shift needs to be in the direction pointing from the avalanche
			// position on the wire to the DOCA point on the rt. We find this
			// by first finding this vector for the rt and then seeing if
			// the "shift" vector is generally pointing in the same direction
			// as it or not by checking if thier dot product is positive or
			// negative.
			double u = rt->GetLastDistAlongWire();
			DVector3 pos_wire = wire->origin + u*wire->udir;
			DVector3 pos_diff = pos_doca-pos_wire;
			if(shift.Dot(pos_diff)<0.0)shift = -shift;
			
			shifts.push_back(shift);
			errs.push_back(SIGMA_CDC);
		}
	}
	
	// Get the chi-squared vector for the initial reference trajectory
	vector<double> chisqv_initial;
	double initial_chisq = ChiSq(rt, wires, shifts, errs, chisqv_initial);

	// Because we have a non-linear system, we must take the derivatives
	// numerically.
	//
	// Note that in the calculations of the deltas below,
	// the change in state should be set first and the value
	// of deltas[...] calculated from that. See Numerical
	// Recipes in C 2nd ed. section 5.7 ppg. 186-189.
	const int Nparameters = 5;
	double deltas[Nparameters];
	DMatrix state(5,1);
	switch(Nparameters){
		case 5: state[state_v	][0] = 0.0;
		case 4: state[state_x	][0] = 0.0;
		case 3: state[state_pz	][0] = mom.Dot(start_step.udir);
		case 2: state[state_py	][0] = mom.Dot(start_step.tdir);
		case 1: state[state_px	][0] = mom.Dot(start_step.sdir);
	}
	
	// For the swimming below, we use a second reference trajectory so as
	// to preserve the original. Set the charge here. The reset of the
	// parameters (starting position and momentum) will be set using
	// values from the state vector.
	tmprt->q = rt->q;
	
	// dpx : tweak by +/- 0.01
	DMatrix state_dpx = state;
	state_dpx[state_px][0] += LEAST_SQUARES_DP;
	deltas[state_px] = state_dpx[state_px][0] - state[state_px][0];
	vector<double> &chisqv_dpx_lo = chisqv_initial;
	vector<double> chisqv_dpx_hi;
	double chisq_dpx = ChiSq(state_dpx, &start_step, tmprt, wires, shifts, errs, chisqv_dpx_hi);

	// dpy : tweak by +/- 0.01
	DMatrix state_dpy = state;
	state_dpy[state_py][0] += LEAST_SQUARES_DP;
	deltas[state_py] = state_dpy[state_py][0] - state[state_py][0];
	vector<double> &chisqv_dpy_lo = chisqv_initial;
	vector<double> chisqv_dpy_hi;
	double chisq_dpy = ChiSq(state_dpy, &start_step, tmprt, wires, shifts, errs, chisqv_dpy_hi);

	// dpz : tweak by +/- 0.01
	DMatrix state_dpz = state;
	state_dpz[state_pz][0] += LEAST_SQUARES_DP;
	deltas[state_pz] = state_dpz[state_pz][0] - state[state_pz][0];
	vector<double> &chisqv_dpz_lo = chisqv_initial;
	vector<double> chisqv_dpz_hi;
	double chisq_dpz = ChiSq(state_dpz, &start_step, tmprt, wires, shifts, errs, chisqv_dpz_hi);

	// dv : tweak by +/- 0.01
	DMatrix state_dv = state;
	state_dv[state_v][0] += LEAST_SQUARES_DX;
	deltas[state_v] = state_dv[state_v][0] - state[state_v][0];
	vector<double> &chisqv_dv_lo = chisqv_initial;
	vector<double> chisqv_dv_hi;
	double chisq_dv = ChiSq(state_dv, &start_step, tmprt, wires, shifts, errs, chisqv_dv_hi);

	// dx : tweak by +/- 0.01
	DMatrix state_dx = state;
	state_dx[state_x][0] += LEAST_SQUARES_DX;
	deltas[state_x] = state_dx[state_x][0] - state[state_x][0];
	vector<double> &chisqv_dx_lo = chisqv_initial;
	vector<double> chisqv_dx_hi;
	double chisq_dx = ChiSq(state_dx, &start_step, tmprt, wires, shifts, errs, chisqv_dx_hi);

	// Make a list of "clean" hits. Ones that are less than 4 sigma and
	// are not "nan" for the initial as well as the tweaked  cases.
	vector<bool> good;
	unsigned int Ngood=0;
	unsigned int Nhits = wires.size();
	for(unsigned int i=0; i<Nhits; i++){
		// Here, we cut hits based on the difference between the track and the 
		// drift time in units of sigmas. However, on the first pass, the fit won't
		// be good and the sigmas which represent the best-fit track to drift
		// time conditions underestimate the quantity. Subsequent passes should
		// have better and better parameters and so the sigmas should be more
		// and more appropriate.
		//
		// To accomodate this, we scale the CHISQ_MAX_RESI_SIGMAS value based
		// on the initial_chisq value since that gives us a handle on how
		// go these hits match to the track overall.
		double max=CHISQ_MAX_RESI_SIGMAS*(initial_chisq>1.0 ? initial_chisq:1.0);
		double sigma;
//_DBG_<<"sigma="<<chisqv_initial[i]<<" : "<<chisqv_dpx_hi[i]<<" : "<<chisqv_dpy_hi[i]<<" : "<<chisqv_dpz_hi[i]<<" : "<<chisqv_dx_hi[i]<<" : "<<chisqv_dv_hi[i]<<endl;
		sigma = chisqv_initial[i]; if(!finite(sigma) || fabs(sigma)>max){good.push_back(false); continue;}
		sigma = chisqv_dpx_hi[i];  if(!finite(sigma) || fabs(sigma)>max){good.push_back(false); continue;}
		sigma = chisqv_dpy_hi[i];  if(!finite(sigma) || fabs(sigma)>max){good.push_back(false); continue;}
		sigma = chisqv_dpz_hi[i];  if(!finite(sigma) || fabs(sigma)>max){good.push_back(false); continue;}
		sigma = chisqv_dx_hi[i];   if(!finite(sigma) || fabs(sigma)>max){good.push_back(false); continue;}
		sigma = chisqv_dv_hi[i];   if(!finite(sigma) || fabs(sigma)>max){good.push_back(false); continue;}

		good.push_back(true);
		Ngood++;
	}
	if(Ngood<LEAST_SQUARES_MIN_HITS){
		if(debug_level>2)_DBG_<<" Bad number of good distance calculations!(Ngood="<<Ngood<<" LEAST_SQUARES_MIN_HITS="<<LEAST_SQUARES_MIN_HITS<<")"<<endl;
		return FIT_FAILED;
	}

	// Build "F" matrix of derivatives
	DMatrix F(Ngood,Nparameters);
	for(unsigned int i=0,j=0; j<Nhits; j++){
		if(!good[j])continue; // skip bad hits
		switch(Nparameters){
			// Note: This is a funny way to use a switch!
			case 5: F[i][state_v ] = pow(errs[j],1.0)*(chisqv_dv_hi[j]-chisqv_dv_lo[j])/deltas[state_v];
			case 4: F[i][state_x ] = pow(errs[j],1.0)*(chisqv_dx_hi[j]-chisqv_dx_lo[j])/deltas[state_x];
			case 3: F[i][state_pz] = pow(errs[j],1.0)*(chisqv_dpz_hi[j]-chisqv_dpz_lo[j])/deltas[state_pz];
			case 2: F[i][state_py] = pow(errs[j],1.0)*(chisqv_dpy_hi[j]-chisqv_dpy_lo[j])/deltas[state_py];
			case 1: F[i][state_px] = pow(errs[j],1.0)*(chisqv_dpx_hi[j]-chisqv_dpx_lo[j])/deltas[state_px];
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
	
	// Measurement vector. This contains the residuals.
	DMatrix m(Ngood,1);
	for(unsigned int i=0,j=0; j<Nhits; j++){
		if(!good[j])continue; // skip bad hits
		// The "measurement" here is the distance of the track from the 
		// shifted wire position where the shift is based on the drift time
		// and a choice of the left-right ambiguity. These are in the 
		// chisqv_XXX vectors in units of sigmas and and are all
		// positive definite. We multiply back in the errors to convert
		// the units back to cm since the errors are handled in the 
		// V matrix. I *think* since we're using a diagonal error matrix
		// we could actually just leave the measurments in units of
		// sigmas and make the V matrix the identity matrix, but we leave
		// it this way for now to make it easier to include covariance
		// later.
		m[i][0] = -chisqv_initial[j]*errs[j];
		V[i][i] = pow(errs[j], 2.0);
		i++;
	}
	DMatrix Vinv(DMatrix::kInverted, V);
	DMatrix B(DMatrix::kInverted, Ft*Vinv*F);
	
	// If the inversion failed altogether then the invalid flag
	// will be set on the matrix. In these cases, were dead.
	if(!B.IsValid()){
		if(debug_level>1)_DBG_<<" -- B matrix invalid"<<endl;
		return FIT_FAILED;
	}

	// The "B" matrix happens to be the covariance matrix of the
	// state parameters. A problem sometimes occurs where one or
	// more elements of B are very large. This tends to happen
	// when a column of F is essentially all zeros making the
	// matrix un-invertable. What we should really do in these
	// cases is check beforehand and drop the bad column(s)
	// before trying to invert. That will add complication that
	// I don't want to introduce quite yet. What we do now
	// is check for it and punt rather than return a nonsensical
	// value.
	if(B.E2Norm() > LEAST_SQUARES_MAX_E2NORM){
		if(debug_level>1)_DBG_<<" -- B matrix E2Norm out of range(E2Norm="<<B.E2Norm()<<" max="<<LEAST_SQUARES_MAX_E2NORM<<")"<<endl;
		return FIT_FAILED;
	}
	
	// Copy the B matrix into last_covariance to later copy into DTrack
	last_covariance.ResizeTo(B);
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
	int max_trys = 4;
	DMatrix new_state(5,1);
	for(; Ntrys<max_trys; Ntrys++){

		for(int i=0; i<Nparameters; i++)new_state[i][0] = state[i][0] + delta_state[i][0]*lambda;

		vector<double> new_chisqv;
		chisq = ChiSq(new_state, &start_step, tmprt, wires, shifts, errs, new_chisqv);
		
		// If we're at a lower chi-sq then we're done
		if(debug_level>4)_DBG_<<" -- initial_chisq="<<initial_chisq<<"  new chisq="<<chisq<<" nhits="<<new_chisqv.size()<<" p="<<tmprt->swim_steps[0].mom.Mag()<<"  lambda="<<lambda<<endl;
		if(chisq-initial_chisq < 0.1)break;
		
		// Chi-sq was increased, try a smaller step on the next iteration
		lambda/=2.0;
	}
	
	// If we failed to find a better Chi-Sq above, maybe we were looking 
	// in the wrong direction(??) Try looking in the opposite direction.
	if(Ntrys>=max_trys){
		lambda = -lambda;
		for(; Ntrys<2*max_trys; Ntrys++){

			for(int i=0; i<Nparameters; i++)new_state[i][0] = state[i][0] + delta_state[i][0]*lambda;

			vector<double> new_chisqv;
			chisq = ChiSq(new_state, &start_step, tmprt, wires, shifts, errs, new_chisqv);
			
			// If we're at a lower chi-sq then we're done
			if(debug_level>4)_DBG_<<" -- initial_chisq="<<initial_chisq<<"  new chisq="<<chisq<<" nhits="<<chisqv.size()<<"  lambda="<<lambda<<endl;
			if(chisq-initial_chisq < 0.1)break;
			
			// Chi-sq was increased, try a smaller step on the next iteration
			lambda/=2.0;
		}
	}

	// If we failed to make a step to a smaller chi-sq then signal
	// that the fit failed completely.
	if(Ntrys>=2*max_trys){
		if(debug_level>1)_DBG_<<"Chisq only increased (both directions searched!)"<<endl;
		return FIT_FAILED;
	}

	// Re-swim reference trajectory using these parameters and re-calc chisq
	vector<double> new_chisqv;
	chisq = ChiSq(new_state, &start_step, rt, wires, shifts, errs, new_chisqv);

//_DBG_<<"chisq="<<chisq<<endl;

	return FIT_OK;
}


#if 0
//------------------
// ChiSq
//------------------
double DTrack_factory_ALT1::ChiSq(double q, DMatrix &state, const swim_step_t *start_step, DReferenceTrajectory *rt)
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
double DTrack_factory_ALT1::ChiSq(double q, const DVector3 &pos, const DVector3 &mom, DReferenceTrajectory *rt)
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
double DTrack_factory_ALT1::ChiSq(double q, DReferenceTrajectory *rt)
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
		//if(fabs(chisqv[i]/sigmav[i])>CHISQ_MAX_RESI_SIGMAS)continue;
		//if(fabs(chisqv[i]/sigmav[i])>chisq_max_resi_sigmas)continue;

		chisq+=pow(chisqv[i]/sigmav[i], 2.0);
		Ngood_chisq_hits += 1.0;
	}
	chisq/=Ngood_chisq_hits;
	if(debug_level>10)_DBG_<<"Chisq:  Ngood_chisq_hits="<<Ngood_chisq_hits<<"  chisq="<<chisq<<endl;

	return chisq;
}

//------------------
// LeastSquares
//------------------
DTrack_factory_ALT1::fit_status_t DTrack_factory_ALT1::LeastSquares(DVector3 &start_pos, DVector3 &start_mom, DReferenceTrajectory *rt, DVector3 &vertex_pos, DVector3 &vertex_mom, double &chisq)
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
		if(debug_level>2)_DBG_<<" Bad number of good distance calculations!"<<endl;
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
	if(!B.IsValid()){
		if(debug_level>1)_DBG_<<" -- B matrix invalid"<<endl;
		return FIT_FAILED;
	}

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
	if(B.E2Norm() > LEAST_SQUARES_MAX_E2NORM){
		if(debug_level>1)_DBG_<<" -- B matrix E2Norm out of range(E2Norm="<<B.E2Norm()<<" max="<<LEAST_SQUARES_MAX_E2NORM<<")"<<endl;
		return FIT_FAILED;
	}
	
	// Copy the B matrix into last_covariance to later copy into DTrack
	last_covariance.ResizeTo(B);
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
		if(debug_level>4)_DBG_<<" -- chisq="<<chisq<<"  new_chisq="<<new_chisq<<" nhits="<<chisqv.size()<<"  lambda="<<lambda<<endl;
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
		if(debug_level>4)_DBG_<<" -- chisq="<<chisq<<"  new_chisq="<<new_chisq<<" nhits="<<chisqv.size()<<"  lambda="<<lambda<<endl;
		if(new_chisq-chisq < 0.01)break;

			// Chi-sq was increased, try a smaller step on the next iteration
			lambda*=2.0;
		}
	}

	// If we failed to make a step to a smaller chi-sq then signal
	// that the fit failed completely.
	if(Ntrys>=max_trys){
		if(debug_level>1)_DBG_<<"Chisq only increased (both directions searched!)"<<endl;
		return FIT_FAILED;
	}
	
	chisq = new_chisq;

	return FIT_OK;
}
#endif

//------------------
// FillDebugHists
//------------------
void DTrack_factory_ALT1::FillDebugHists(DReferenceTrajectory *rt, DVector3 &vertex_pos, DVector3 &vertex_mom)
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
const string DTrack_factory_ALT1::toString(void)
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


