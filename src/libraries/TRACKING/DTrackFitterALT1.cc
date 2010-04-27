// $Id$
//
//    File: DTrackFitterALT1.cc
// Created: Tue Sep  2 11:18:22 EDT 2008
// Creator: davidl
//
namespace{const char* GetMyID(){return "$Id$";}}

#include <math.h>

#include <TROOT.h>

#include <DVector3.h>
#include <DMatrix.h>

#include <JANA/JEventLoop.h>
using namespace jana;

#include "GlueX.h"
#include "DANA/DApplication.h"
#include "DMagneticFieldStepper.h"
#include "DTrackCandidate.h"
#include "DTrackFitterALT1.h"
#include "CDC/DCDCTrackHit.h"
#include "FDC/DFDCPseudo.h"
#include "DReferenceTrajectory.h"
#include "DMCThrown.h"
#include "DMCTrackHit.h"
#include "DTrackHitSelectorTHROWN.h"

extern double GetCDCCovariance(int layer1, int layer2);
extern double GetFDCCovariance(int layer1, int layer2);
extern double GetFDCCathodeCovariance(int layer1, int layer2);


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
// DTrackFitterALT1   (Constructor)
//------------------
DTrackFitterALT1::DTrackFitterALT1(JEventLoop *loop):DTrackFitter(loop)
{
	this->loop = loop;

	// Define target center
	target = new DCoordinateSystem();
	target->origin.SetXYZ(0.0, 0.0, 65.0);
	target->sdir.SetXYZ(1.0, 0.0, 0.0);
	target->tdir.SetXYZ(0.0, 1.0, 0.0);
	target->udir.SetXYZ(0.0, 0.0, 1.0);
	target->L = 30.0;

	// DReferenceTrajectory objects
	rt = new DReferenceTrajectory(bfield);
	tmprt = new DReferenceTrajectory(bfield);
	
	DEBUG_HISTS = false;
	DEBUG_LEVEL = 0;
	MAX_CHISQ_DIFF = 1.0E-2;
	MAX_FIT_ITERATIONS = 50;
	SIGMA_CDC = 0.0150;
	CDC_USE_PARAMETERIZED_SIGMA = true;
	SIGMA_FDC_ANODE = 0.0200;
	SIGMA_FDC_CATHODE = 0.0200;
	CHISQ_GOOD_LIMIT = 2.0;
	LEAST_SQUARES_DP = 0.0001;
	LEAST_SQUARES_DX = 0.010;
	LEAST_SQUARES_MIN_HITS = 3;
	LEAST_SQUARES_MAX_E2NORM = 1.0E20;
	DEFAULT_MASS = rt->GetMass(); // Get default mass from DReferenceTrajectory class itself
	TARGET_CONSTRAINT = false;
	LR_FORCE_TRUTH = false;
	USE_MULS_COVARIANCE = true;
	USE_FDC = true;
	USE_CDC = true;
	USE_FDC_CATHODE = true;
	MATERIAL_MAP_MODEL = "DGeometry";
	
	gPARMS->SetDefaultParameter("TRKFIT:DEBUG_HISTS",					DEBUG_HISTS);
	gPARMS->SetDefaultParameter("TRKFIT:DEBUG_LEVEL",					DEBUG_LEVEL);
	gPARMS->SetDefaultParameter("TRKFIT:MAX_CHISQ_DIFF",				MAX_CHISQ_DIFF);
	gPARMS->SetDefaultParameter("TRKFIT:MAX_FIT_ITERATIONS",			MAX_FIT_ITERATIONS);
	gPARMS->SetDefaultParameter("TRKFIT:SIGMA_CDC",						SIGMA_CDC);
	gPARMS->SetDefaultParameter("TRKFIT:CDC_USE_PARAMETERIZED_SIGMA",	CDC_USE_PARAMETERIZED_SIGMA);
	gPARMS->SetDefaultParameter("TRKFIT:SIGMA_FDC_ANODE",				SIGMA_FDC_ANODE);
	gPARMS->SetDefaultParameter("TRKFIT:SIGMA_FDC_CATHODE",			SIGMA_FDC_CATHODE);
	gPARMS->SetDefaultParameter("TRKFIT:CHISQ_GOOD_LIMIT",			CHISQ_GOOD_LIMIT);
	gPARMS->SetDefaultParameter("TRKFIT:LEAST_SQUARES_DP",			LEAST_SQUARES_DP);
	gPARMS->SetDefaultParameter("TRKFIT:LEAST_SQUARES_DX",			LEAST_SQUARES_DX);
	gPARMS->SetDefaultParameter("TRKFIT:LEAST_SQUARES_MIN_HITS",	LEAST_SQUARES_MIN_HITS);
	gPARMS->SetDefaultParameter("TRKFIT:LEAST_SQUARES_MAX_E2NORM",	LEAST_SQUARES_MAX_E2NORM);		
	gPARMS->SetDefaultParameter("TRKFIT:DEFAULT_MASS",					DEFAULT_MASS);
	gPARMS->SetDefaultParameter("TRKFIT:TARGET_CONSTRAINT",			TARGET_CONSTRAINT);
	gPARMS->SetDefaultParameter("TRKFIT:LR_FORCE_TRUTH",				LR_FORCE_TRUTH);
	gPARMS->SetDefaultParameter("TRKFIT:USE_MULS_COVARIANCE",		USE_MULS_COVARIANCE);
	gPARMS->SetDefaultParameter("TRKFIT:USE_FDC",						USE_FDC);
	gPARMS->SetDefaultParameter("TRKFIT:USE_CDC",						USE_CDC);
	gPARMS->SetDefaultParameter("TRKFIT:USE_FDC_CATHODE",				USE_FDC_CATHODE);
	gPARMS->SetDefaultParameter("TRKFIT:MATERIAL_MAP_MODEL",			MATERIAL_MAP_MODEL);
	
	// Set default mass based on configuration parameter
	rt->SetMass(DEFAULT_MASS);
	tmprt->SetMass(DEFAULT_MASS);
	
	// Set the geometry object pointers based on the material map
	// specified via configuration parameter.
	if(MATERIAL_MAP_MODEL=="DRootGeom"){
		rt->SetDRootGeom(RootGeom);
		tmprt->SetDRootGeom(RootGeom);
	}else if(MATERIAL_MAP_MODEL=="DGeometry"){
		rt->SetDGeometry(geom);
		tmprt->SetDGeometry(geom);
	}else if(MATERIAL_MAP_MODEL=="NONE"){
		rt->SetDRootGeom(NULL);
		tmprt->SetDRootGeom(NULL);
		rt->SetDGeometry(NULL);
		tmprt->SetDGeometry(NULL);
	}else{
		_DBG_<<"ERROR: Invalid value for TRKFIT:MATERIAL_MAP_MODEL (=\""<<MATERIAL_MAP_MODEL<<"\")"<<endl;
		exit(-1);
	}

	cout<<__FILE__<<":"<<__LINE__<<"-------------- Least Squares TRACKING --------------"<<endl;

	DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
	if(!dapp){
		_DBG_<<"Cannot get DApplication from JEventLoop! (are you using a JApplication based program?)"<<endl;
		return;
	}
	
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
		chisq_vs_p_vs_theta = (TH2F*)gROOT->FindObject("chisq_vs_p_vs_theta");
		lambda = (TH1F*)gROOT->FindObject("lambda");

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
		if(!lambda)lambda = new TH1F("lambda","Scaling factor #lambda for Newton-Raphson calculated step", 2048, -2.0, 2.0);

		dapp->Unlock();
	}
}

//------------------
// DTrackFitterALT1   (Destructor)
//------------------
DTrackFitterALT1::~DTrackFitterALT1()
{
	if(rt)delete rt;
	if(tmprt)delete tmprt;
	if(target)delete target;
}

//------------------
// FitTrack
//------------------
DTrackFitter::fit_status_t DTrackFitterALT1::FitTrack(void)
{
	/// Fit a track candidate
	
	// Debug message
	if(DEBUG_LEVEL>2){
		_DBG_<<"cdchits.size="<<cdchits.size()<<"  fdchits.size="<<fdchits.size()<<endl;
		if(fit_type==kTimeBased)_DBG_<<"============= Time-based ==============="<<endl;
		if(fit_type==kWireBased)_DBG_<<"============= Wire-based ==============="<<endl;
	}
	
	// Set mass for our DReferenceTrajectory objects
	rt->SetMass(input_params.mass());
	tmprt->SetMass(input_params.mass());
	
	// Do material boundary checking only for time-based tracking
	// and only if beta*gamma is less than 1.0
	double betagamma = input_params.momentum().Mag()/input_params.mass();
	bool check_material_boundaries = fit_type==kTimeBased && betagamma<=1.0;
	rt->SetCheckMaterialBoundaries(check_material_boundaries);
	tmprt->SetCheckMaterialBoundaries(check_material_boundaries);
	
	// Swim reference trajectory
	fit_status = kFitNotDone; // initialize to a safe default
	DVector3 mom = input_params.momentum();
	if(mom.Mag()>12.0)mom.SetMag(8.0);	// limit crazy big momenta from candidates
	rt->Swim(input_params.position(), mom, input_params.charge());

	if(DEBUG_LEVEL>1){
		double chisq;
		int Ndof;
		ChiSq(fit_type, rt, &chisq, &Ndof);
		_DBG_<<"starting parameters for fit: chisq="<<chisq<<" Ndof="<<Ndof<<" (chisq/Ndof="<<chisq/(double)Ndof<<") p="<<rt->swim_steps[0].mom.Mag()<<endl;
	}

	// Iterate over fitter until it converges or we reach the 
	// maximum number of iterations
	double last_chisq = 1.0E6;
	int last_Ndof=1;
	int Niterations;
	for(Niterations=0; Niterations<MAX_FIT_ITERATIONS; Niterations++){
	
		hitsInfo hinfo;
		GetHits(fit_type, rt, hinfo);
		
		// Optionally use the truth information from the Monte Carlo
		// to force the correct LR choice for each hit. This is
		// for debugging (obviously).
		if(LR_FORCE_TRUTH && fit_type==kTimeBased)ForceLRTruth(loop, rt, hinfo);
	
		switch(fit_status = LeastSquaresB(hinfo, rt)){
			case kFitSuccess:
				if(DEBUG_LEVEL>2)_DBG_<<"Fit succeeded ----- "<<endl;
				break;
			case kFitNoImprovement:
				if(DEBUG_LEVEL>2)_DBG_<<"Unable to improve on input parameters (chisq="<<chisq<<")"<<endl;
				break;
			case kFitFailed:
				if(DEBUG_LEVEL>2)_DBG_<<"Fit failed on iteration "<<Niterations<<endl;
				if(Niterations>0){
					if(DEBUG_LEVEL>2)_DBG_<<"Number of iterations >4. Trying to keep fit from last iteration... "<<endl;
					chisq=last_chisq;
					Ndof = last_Ndof;
					fit_status = kFitNoImprovement;
					break;
				}
				return fit_status;
				break;
			case kFitNotDone: // avoid compiler warnings
				break;
		}

		// Fill debug histos
		if(DEBUG_HISTS){
			chisq_vs_pass->Fill(Niterations+1, chisq/Ndof);
			dchisq_vs_pass->Fill(Niterations+1, last_chisq/last_Ndof - chisq/Ndof);			
		}

		// If the chisq is too large, then consider it a hopeless cause
		if(chisq/Ndof>1.0E4 || !finite(chisq/Ndof)){
			if(DEBUG_LEVEL>3)_DBG_<<"Fit chisq too large on iteration "<<Niterations<<endl;
			return fit_status=kFitFailed;
		}
		
		// Check if we converged
		if(fabs(last_chisq/last_Ndof-chisq/Ndof) < MAX_CHISQ_DIFF){
			if(DEBUG_LEVEL>3)_DBG_<<"Fit chisq change below threshold fabs("<<last_chisq/last_Ndof<<" - "<<chisq/Ndof<<")<"<<MAX_CHISQ_DIFF<<endl;
			break;
		}
		last_chisq = chisq;
		last_Ndof = Ndof;
	}

	if(DEBUG_LEVEL>1)_DBG_<<" Niterations="<<Niterations<<"  chisq="<<chisq<<" Ndof="<<Ndof<<" (chisq/Ndof="<<chisq/(double)Ndof<<") p="<<rt->swim_steps[0].mom.Mag()<<endl;
	
	// At this point we must decide whether the fit succeeded or not.
	// We'll consider the fit a success if any of the following is true:
	//
	// 1. We got through at least one iteration in the above loop
	// 2. The chi-sq is less than CHISQ_GOOD_LIMIT
	// 3. MAX_FIT_ITERATIONS is zero (for debugging)
	bool fit_succeeded = false;
	if(Niterations>0)fit_succeeded = true;
	if(chisq<=CHISQ_GOOD_LIMIT)fit_succeeded = true;
	if(MAX_FIT_ITERATIONS==0){fit_succeeded = true; fit_status=kFitSuccess;}
	if(!fit_succeeded)return fit_status=kFitFailed;
	
	if(DEBUG_LEVEL>9){
		_DBG_<<"-------- Final Chisq for track = "<<this->chisq<<" Ndof="<<this->Ndof<<endl;
		double chisq;
		int Ndof;
		ChiSq(fit_type, rt, &chisq, &Ndof);
		_DBG_<<"-------- Check chisq = "<<chisq<<" Ndof="<<Ndof<<endl;
	}
	
	if(DEBUG_LEVEL>2){
		_DBG_<<"start POCA pos: "; rt->swim_steps->origin.Print();
		_DBG_<<"start POCA mom: "; rt->swim_steps->mom.Print();
		ChiSq(fit_type, rt, &chisq, &Ndof);
		_DBG_<<"start POCA chisq="<<chisq<<" Ndof="<<Ndof<<" (chisq/Ndof="<<chisq/(double)Ndof<<") p="<<rt->swim_steps[0].mom.Mag()<<endl;
	}

	// At this point, the first point on the rt is likely not the POCA to the
	// beamline. To find that, we swim the particle in, towards the target from
	// the swim step closest to the first wire hit. We make sure to set the tmprt to
	// energy *gain* for the backwards swim.
	DVector3 vertex_pos = rt->swim_steps->origin;
	DVector3 vertex_mom = rt->swim_steps->mom;
	const DCoordinateSystem *first_wire = NULL;
	if(cdchits.size()>0) first_wire=cdchits[0]->wire;
	else if(fdchits.size()>0) first_wire=fdchits[0]->wire;
	if(first_wire){
		DReferenceTrajectory::swim_step_t *step = rt->FindClosestSwimStep(first_wire);
		if(step){
#if 0
			tmprt->SetPLossDirection(DReferenceTrajectory::kBackward);
			tmprt->Swim(step->origin, -step->mom, -rt->q, 100.0, target);
			tmprt->DistToRT(target);
			tmprt->GetLastDOCAPoint(vertex_pos, vertex_mom);
			vertex_mom = -vertex_mom;
			tmprt->SetPLossDirection(DReferenceTrajectory::kForward);
#endif
			// Swim back towards target (note: we may already be past it!)
			tmprt->SetPLossDirection(DReferenceTrajectory::kBackward);
			tmprt->Swim(rt->swim_steps->origin, -rt->swim_steps->mom, -rt->q, 100.0, target);

			// Swim swim in farward direction to target (another short swim!)
			tmprt->SetPLossDirection(DReferenceTrajectory::kForward);
			vertex_pos = tmprt->swim_steps[tmprt->Nswim_steps-1].origin;
			vertex_mom = tmprt->swim_steps[tmprt->Nswim_steps-1].mom;
			tmprt->Swim(vertex_pos, -vertex_mom, rt->q, 100.0, target);
			tmprt->DistToRT(target);
			tmprt->GetLastDOCAPoint(vertex_pos, vertex_mom);

			// Now, swim out the rt one last time such that it starts at the POCA to
			// the beamline. We need to do this so that we can calculate the
			// chisq/Ndof based on the vertex parameters
			rt->Swim(vertex_pos, vertex_mom);
			ChiSq(fit_type, rt, &this->chisq, &this->Ndof);
		}
	}else{
		_DBG_<<"NO WIRES IN CANDIDATE!! (event "<<eventnumber<<")"<<endl;
	}

	if(DEBUG_LEVEL>2){
		_DBG_<<"end POCA pos: "; rt->swim_steps->origin.Print();
		_DBG_<<"end POCA mom: "; rt->swim_steps->mom.Print();
		ChiSq(fit_type, rt, &chisq, &Ndof);
		_DBG_<<"end POCA chisq="<<chisq<<" Ndof="<<Ndof<<" (chisq/Ndof="<<chisq/(double)Ndof<<") p="<<rt->swim_steps[0].mom.Mag()<<endl;
	}


#if 0	
	// At this point, the first point on the rt is likely not the POCA to the
	// beamline. To find that, we swim the particle in, towards the target from
	// about 10 cm out along the trajectory. We make sure to set the tmprt to
	// energy gain for the backwards swim.
	DReferenceTrajectory::swim_step_t *step = rt->swim_steps;
	for(int i=0; i<rt->Nswim_steps-1; i++, step++){if(step->s>10.0)break;}
	tmprt->SetPLossDirection(DReferenceTrajectory::kBackward);
	tmprt->Swim(step->origin, -step->mom, -rt->q, 100.0, target);
	tmprt->SetPLossDirection(DReferenceTrajectory::kForward);
	double s;
	tmprt->DistToRT(target, &s);
	DVector3 vertex_mom, vertex_pos;
	tmprt->GetLastDOCAPoint(vertex_pos, vertex_mom);
	vertex_mom = -vertex_mom;
#endif

	// Copy final fit parameters into TrackFitter classes data members. Note that the chisq and Ndof
	// members are copied in during the ChiSq() method call in LeastSquaresB().
	fit_params.setPosition(vertex_pos);
	fit_params.setMomentum(vertex_mom);
	fit_params.setCharge(rt->q);
	fit_params.setMass(input_params.mass());
	cdchits_used_in_fit = cdchits;
	fdchits_used_in_fit = fdchits;

	// Fill debugging histos if requested
	if(DEBUG_HISTS){
		Npasses->Fill(Niterations);
		initial_chisq_vs_Npasses->Fill(Niterations, chisq);
		FillDebugHists(rt, vertex_pos, vertex_mom);
	}
	
	// Debugging messages
	if(DEBUG_LEVEL>1)_DBG_<<" -- Fit succeeded: Final params after vertex adjustment: q="<<rt->q<<" p="<<vertex_mom.Mag()<<" theta="<<vertex_mom.Theta()*57.3<<" phi="<<vertex_mom.Phi()*57.3<<"  chisq="<<chisq<<" Ndof="<<Ndof<<" (chisq/Ndof="<<chisq/(double)Ndof<<")"<<endl;

	return fit_status;
}

//------------------
// ChiSq
//------------------
double DTrackFitterALT1::ChiSq(fit_type_t fit_type, DReferenceTrajectory *rt, double *chisq_ptr, int *dof_ptr, vector<pull_t> *pulls_ptr)
{
	/// Calculate the chisq for the given reference trajectory based on the hits
	/// currently registered through the DTrackFitter base class into the cdchits
	/// and fdchits vectors. This does not get called by the core part of the
	/// fitter, but is used, rather, to give an "independent" value of the
	/// chisq based only on a reference trajectory and the input hits. It is
	/// called after the fit to calculate the final chisq.

	hitsInfo hinfo;
	GetHits(fit_type, rt, hinfo);
	if(LR_FORCE_TRUTH && fit_type==kTimeBased)ForceLRTruth(loop, rt, hinfo);

	vector<resiInfo> residuals;
	GetResiInfo(rt, hinfo, residuals);
	
	double my_chisq = ChiSq(residuals, chisq_ptr, dof_ptr);
	if(pulls_ptr)*pulls_ptr = pulls;
	
	return my_chisq;
}

//------------------
// ChiSq
//------------------
double DTrackFitterALT1::ChiSq(vector<resiInfo> &residuals, double *chisq_ptr, int *dof_ptr)
{
	// The residuals vector contains a list of only good hits, both
	// anode and cathode, along with their measurment errors. Use these to fill
	// the resiv and cov_meas matrices now. Fill in the cov_muls below.
	int Nmeasurements = (int)residuals.size();
	resiv.ResizeTo(Nmeasurements, 1);
	cov_meas.ResizeTo(Nmeasurements, Nmeasurements);
	cov_muls.ResizeTo(Nmeasurements, Nmeasurements);
	pulls.clear();
	
	// Return "infinite" chisq if an empty residuals vector was passed.
	if(Nmeasurements<1){
		if(chisq_ptr)*chisq_ptr=1.0E6;
		if(dof_ptr)*dof_ptr=1;
		if(DEBUG_LEVEL>3)_DBG_<<"Chisq called with empty residuals vector!"<<endl;
		return 1.0E6;
	}
	
	for(unsigned int i=0; i<residuals.size(); i++){
		resiInfo &ri = residuals[i];
		
		resiv[i][0] = ri.resi;
		cov_meas[i][i] = pow(ri.err, 2.0);
	}
	
	// Fill in the cov_muls matrix.
	cov_muls.Zero();
	if(USE_MULS_COVARIANCE){
		// Find the first sensitive detector hit. Ignore
		// the target point if present by checking layer!=0
		const swim_step_t *step0 = NULL;
		for(unsigned int i=0; i<residuals.size(); i++){
			resiInfo &ri = residuals[i];
			if(ri.layer<1)continue;
			if(!ri.step)continue;
			if( (step0==NULL) || (ri.step->s < step0->s) ){
				step0 = ri.step;
			}
		}

		if(step0){
			// Loop over all residuals
			for(unsigned int i=0; i<residuals.size(); i++){
				const swim_step_t *stepA = residuals[i].step;
				if(!stepA)continue;
				resiInfo &riA = residuals[i];
				
				// Loop over all residuals
				for(unsigned int j=i; j<residuals.size(); j++){
					const swim_step_t *stepB = residuals[j].step;
					if(!stepB)continue;
					resiInfo &riB = residuals[j];
					
					// Correlations are very hard to get correct given the left-right
					// ambiguity which determines the sign of the correlation. 
					// Because of this, only include the diagonal elements of the
					// covariance matrix. Attempts were made to do otherwise, but they
					// failed. The code is left here in case we need to try it again
					// at a later time.
					if(i!=j)continue;
					
					// Correlations between A and B are due only to material between
					// the first detector and the most upstream of A or B.
					double sA = stepA->s;
					double sB = stepB->s;
					const swim_step_t *step_end = sA<sB ? stepA:stepB;

					if(step0->s>step_end->s)continue; // Bullet proof
					
					// Need a comment here to explain this ... 
					double itheta02 = step_end->itheta02 - step0->itheta02;
					double itheta02s = step_end->itheta02s - step0->itheta02s;
					double itheta02s2 = step_end->itheta02s2 - step0->itheta02s2;
					double sigmaAB = sA*sB*itheta02 - (sA+sB)*itheta02s + itheta02s2;
					
					// sigmaAB represents the magnitude of the covariance, but not the
					// sign.
					//
					// For wires generally in the same direction, the sign
					// should be positive if they are on the same side of the wire
					// but negative if they are on opposite sides.
					//
					// For wires that are orhogonal (like a CDC is to an FDC wire), the
					// the procedure is not so clear.
					
					// Find LR of this hit.
					// Vector pointing from origin of wire to step position crossed into
					// wire direction gives a vector that will either be pointing
					// generally in the direction of the momentum or opposite to it.
					// Take dot product of above described vectors  for each hit 
					// use it to determine L or R.
					DVector3 crossA = (stepA->origin - riA.hit->wire->origin).Cross(riA.hit->wire->udir);
					DVector3 crossB = (stepB->origin - riB.hit->wire->origin).Cross(riB.hit->wire->udir);
					double sides_angle = crossA.Angle(crossB)*57.3;
					double sign = fabs(sides_angle) < 90.0 ? +1.0:-1.0;
					
					cov_muls[i][j] = cov_muls[j][i] = sign*sigmaAB;
				}
			}
		}
	}

	// Multiply resiv with inverse of total error matrix to get chisq = r^T (cov_meas+cov_muls)^-1 r)
	DMatrix resiv_t(DMatrix::kTransposed, resiv);
	DMatrix W(DMatrix::kInverted, cov_meas+cov_muls);
	DMatrix chisqM(resiv_t*W*resiv);
	double chisq = chisqM[0][0];
	int Ndof = (int)residuals.size() - 5; // assume 5 fit parameters
	weights.ResizeTo(W);
	weights = W;
	
	// Copy pulls into pulls vector
	for(unsigned int i=0; i<residuals.size(); i++){
		double err = sqrt(cov_meas[i][i]+cov_muls[i][i]);
		pulls.push_back(pull_t(resiv[i][0], err));
	}
	
	if(DEBUG_LEVEL>100 || (DEBUG_LEVEL>10 && !finite(chisq)))PrintChisqElements(resiv, cov_meas, cov_muls, weights);

	// If the caller supplied pointers to chisq and dof variables, copy the values into them
	if(chisq_ptr)*chisq_ptr = chisq;
	if(dof_ptr)*dof_ptr = Ndof; // assume 5 fit parameters

	// If it turns out the dof is zero, return 1.0E6. Otherwise, return
	// the chisq/dof
	return Ndof<2 ? 1.0E6:chisq/(double)Ndof;
}

//------------------
// GetHits
//------------------
void DTrackFitterALT1::GetHits(fit_type_t fit_type, DReferenceTrajectory *rt, hitsInfo &hinfo)
{

	// -- Target --
	if(TARGET_CONSTRAINT){
		hitInfo hi;
		hi.wire = target;
		hi.dist = 0.0;
		hi.err = 0.1; // 1mm beam width
		hi.u_dist = 0.0;
		hi.u_err = 0.0;
		hi.good = hi.good_u = false;
		hinfo.push_back(hi);
	}

	// --- CDC ---
	if(USE_CDC){
		for(unsigned int i=0; i<cdchits.size(); i++){
			const DCDCTrackHit *hit = cdchits[i];
			const DCoordinateSystem *wire = hit->wire;
			hitInfo hi;

			hi.wire = wire;
			hi.u_dist = 0.0;
			hi.u_err = 0.0;
			hi.good = hi.good_u = false;

			// Fill in shifts and errs vectors based on whether we're doing
			// hit-based or time-based tracking
			if(fit_type==kWireBased){
				// If we're doing hit-based tracking then only the wire positions
				// are used and the drift time info is ignored.
				// NOTE: The track quality itself goes into the residual resoultion
				// and so we use something larger than the variance over an evenly 
				// illuminated cell size (commented out). The value of 0.31 is
				// empirical from forward (2-40 degree) pi+ tracks when MULS was
				// off. This will likely need to be higher when MULS is on...
				hi.dist = 0.0;
				hi.err = 0.45; // empirical. (see note above)
				//hi.err = 0.8/sqrt(12.0); // variance for evenly illuminated straw
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
				
				// The magnitude of the shift is based on the drift time. The
				// value of the dist member of the DCDCTrackHit object does not
				// subtract out the TOF. This can add 50-100 microns to the
				// resolution in the CDC. Here, we actually can calculate the TOF
				// (for a given mass hypothesis).
				double mass = rt->GetMass();
				double beta = 1.0/sqrt(1.0 + pow(mass/mom_doca.Mag(), 2.0))*2.998E10;
				double tof = s/beta/1.0E-9; // in ns
				hi.dist = hit->dist*((hit->tdrift-tof)/hit->tdrift);
				hi.err = SIGMA_CDC;
				
				// Default is to use parameterized sigma, by allow user to specify
				// using constant sigma.
				if(CDC_USE_PARAMETERIZED_SIGMA){
					// The following is from a fit to Yves' numbers circa Aug 2009. The values fit were
					// resolution (microns) vs. drift distance (mm).
					// par[8] = {699.875, -559.056, 149.391, 25.6929, -22.0238, 4.75091, -0.452373, 0.0163858};
					double x = hi.dist;
					if(x<0.0)x=0.0; // this can be negative due to tof subtraction above.
					//double sigma_d = (699.875) + x*((-559.056) + x*((149.391) + x*((25.6929) + x*((-22.0238) + x*((4.75091) + x*((-0.452373) + x*((0.0163858))))))));
					double sigma_d = 108.55 + 7.62391*x + 556.176*exp(-(1.12566)*pow(x,1.29645));
					hi.err = sigma_d/10000.0;
				}
			}
			hinfo.push_back(hi);
		}
	} // USE_CDC

	// --- FDC ---
	if(USE_FDC){
		for(unsigned int i=0; i<fdchits.size(); i++){
			const DFDCPseudo *hit = fdchits[i];
			const DCoordinateSystem *wire = hit->wire;
			hitInfo hi;

			hi.wire = wire;
			hi.good = hi.good_u = false;

			// Fill in shifts and errs vectors based on whether we're doing
			// hit-based or time-based tracking
			if(fit_type==kWireBased){
				// If we're doing hit-based tracking then only the wire positions
				// are used and the drift time info is ignored.
				// NOTE: The track quality itself goes into the residual resoultion
				// and so we use something larger than the variance over an evenly 
				// illuminated cell size (commented out). The value of 0.30 is
				// empirical from forward (2-40 degree) pi+ tracks when MULS was
				hi.dist = 0.0;
				hi.err = 0.30; // emprical. (see note above)
				//hi.err = 0.5/sqrt(12.0); // variance for evenly illuminated cell - anode
				hi.u_dist = 0.0;
				hi.u_err = 0.0; // variance for evenly illuminated cell - cathode
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
				
				// The magnitude of the shift is based on the drift time. The
				// value of the dist member of the DCDCTrackHit object does not
				// subtract out the TOF. This can add 50-100 microns to the
				// resolution in the CDC. Here, we actually can calculate the TOF
				// (for a given mass hypothesis).
				double mass = rt->GetMass();
				double beta = 1.0/sqrt(1.0 + pow(mass/mom_doca.Mag(), 2.0))*2.998E10;
				double tof = s/beta/1.0E-9; // in ns
				hi.dist = hit->dist*((hit->time-tof)/hit->time);
				hi.err = SIGMA_FDC_ANODE;
				
				if(USE_FDC_CATHODE){
					// Find whether the track is on the "left" or "right" of the wire
					DVector3 shift = wire->udir.Cross(mom_doca);
					if(shift.Mag()!=0.0)shift.SetMag(1.0);
					double u = rt->GetLastDistAlongWire();
					DVector3 pos_wire = wire->origin + u*wire->udir;
					double LRsign = shift.Dot(pos_doca-pos_wire)<0.0 ? +1.0:-1.0;

					// Lorentz corrected poisition along the wire is contained in x,y
					// values, BUT, it is based on a left-right decision of the track
					// segment. This may or may not be the same as the global track. 
					// As such, we have to determine the correction for our track.
					//DVector3 wpos(hit->x, hit->y, wire->origin.Z());
					//DVector3 wdiff = wpos - wire->origin;
					//double u_corr = wire->udir.Dot(wdiff);
					double alpha = mom_doca.Angle(DVector3(0,0,1));
					hi.u_lorentz = LRsign*lorentz_def->GetLorentzCorrection(pos_doca.X(), pos_doca.Y(), pos_doca.Z(), alpha, hi.dist);
					hi.u_dist = hit->s;
					hi.u_err = SIGMA_FDC_CATHODE;
				}else{
					// User specified not to use FDC cathode information in the fit.
					hi.u_dist = 0.0;
					hi.u_err = 0.0; // setting u_err to zero means it's excluded from chi-sq
				}
			}
			hinfo.push_back(hi);
		}
	} // USE_FDC
}


//------------------
// GetResiInfo
//------------------
vector<bool> DTrackFitterALT1::GetResiInfo(DMatrix &state, const swim_step_t *start_step, DReferenceTrajectory *rt, hitsInfo &hinfo, vector<resiInfo> &residuals)
{
	/// Calculate the chi-squared for a track specified by state relative to the
	/// given reference trajectory. This is just a wrapper for 
	/// ChiSq(DReferenceTrajectory *rt, vector<const DCoordinateSystem*> &wires, vector<DVector3> &shifts, vector<double> &errs, vector<double> &chisqv)
	/// that accepts the state vector and re-swims the trajectory.
	
	// In case we need to return early
	int Nhits=0;
	for(unsigned int i=0; i<hinfo.size(); i++){Nhits++; if(hinfo[i].u_err!=0.0)Nhits++;}
	vector<bool> good_none(Nhits, false);
	
	// "v" direction is perpendicular to both the rt direction and the
	// x-direction. See LeastSquares() for more.
	DVector3 vdir = start_step->sdir.Cross(start_step->mom);
	if(vdir.Mag()!=0.0)vdir.SetMag(1.0);

	DVector3 pos =   start_step->origin
						+ state[state_x ][0]*start_step->sdir
						+ state[state_v ][0]*vdir;
	DVector3 mom =   state[state_px][0]*start_step->sdir
						+ state[state_py][0]*start_step->tdir
						+ state[state_pz][0]*start_step->udir;

	if(!rt){
		_DBG_<<"NULL pointer passed for DReferenceTrajectory object!"<<endl;
		return good_none;
	}
	
	if(pos.Mag()>200.0 || fabs(state[state_x ][0])>100.0 || fabs(state[state_v ][0])>100.0){
		if(DEBUG_LEVEL>3)_DBG_<<"state values out of range"<<endl;
		if(DEBUG_LEVEL>6){
			pos.Print();
			mom.Print();
			_DBG_<<"state[state_x ][0]="<<state[state_x ][0]<<"  state[state_v ][0]="<<state[state_x ][0]<<endl;
		}
		return good_none;
	}
	
	// Swim the trajectory with the specified state
	rt->Swim(pos,mom);
	
	// Sometimes, a bad state vector is passed that leads to referance trajectory with no steps.
	if(rt->Nswim_steps<1){
		residuals.clear();
		return good_none;
	}

	return GetResiInfo(rt, hinfo, residuals);
}

//------------------
// GetResiInfo
//------------------
vector<bool> DTrackFitterALT1::GetResiInfo(DReferenceTrajectory *rt, hitsInfo &hinfo, vector<resiInfo> &residuals)
{
	// Make a serialized list of "good" hits as we sparsely fill the residuals
	// container with only the good ones.
	vector<bool> good;

	// Loop over wires hit. Make lists of finite residuals with layer numbers
	// from which to build the actual matrices used to calculate chisq below.
	residuals.clear();
	for(unsigned int i=0; i<hinfo.size(); i++){
		hitInfo &hi = hinfo[i];
		hi.good = hi.good_u = false;

		const DCoordinateSystem *wire = hi.wire;
		
		// Figure out whether this is a CDC or FDC wire. Note that
		// it could be neither if the target constraint is used.
		const DFDCWire *fdcwire = dynamic_cast<const DFDCWire*>(wire);
		const DCDCWire *cdcwire = dynamic_cast<const DCDCWire*>(wire);
	
		// Get distance of the wire from the reference trajectory and the
		// distance s along the track to the point of closest approach.
		double s;
		double d = rt->DistToRT(wire, &s);
		
		// Residual. If we're on the correct side of the wire, then this is
		// dist-doca. If we're on the incorrect side of the wire, then this
		// is -dist-doca. Prior to calling us, the value of hi.dist will have
		// a sign that has already been assigned to indicate the side of the wire
		// the track is believed to be on.
		double resi = hi.dist - d;
		if(finite(resi)){
			resiInfo ri;
			ri.hit = &hi;
			ri.layer = cdcwire ? cdcwire->ring:(fdcwire ? fdcwire->layer:0);
			ri.resi_type = cdcwire ? resi_type_cdc_anode:(fdcwire ? resi_type_fdc_anode:resi_type_other);
			ri.resi = resi;
			ri.err = hi.err;
			ri.step = rt->GetLastSwimStep();
			hi.good = true;
			residuals.push_back(ri);
			good.push_back(true);
		}else{
			good.push_back(false);
		}
		
		// Also add residual along the wire. If the value of u_err is zero
		// that indicates no measurement was made along the wire for this hit.
		if(hi.u_err!=0.0){
			// The sign of hi.dist indicates whether we want to treat this hit as
			// being on the same side of the wire as the track(+) or the opposite
			// side (-). In the latter case, we need to apply the Lorentz correction
			// to the position along the wire in the opposite direction than we
			// would otherwise. Set the sign of the Lorentz deflection based on the
			// sign of hi.dist.
			double LRsign = hi.dist<0.0 ? -1.0:1.0;
		
			double u = rt->GetLastDistAlongWire();
			double u_corrected = hi.u_dist + LRsign*hi.u_lorentz;
			double resic = u - u_corrected;
			if(finite(resic)){
				resiInfo ri;
				ri.hit = &hi;
				ri.layer = fdcwire ? fdcwire->layer:0;
				ri.resi_type =fdcwire ? resi_type_fdc_cathode:resi_type_other;
				ri.resi = resic;
				ri.err = hi.u_err;
				ri.step = rt->GetLastSwimStep();
				hi.good_u = true;
				residuals.push_back(ri);
				good.push_back(true);
			}else{
				good.push_back(false);
			}
		}
	}
	
	return good;
}

//------------------
// LeastSquaresB
//------------------
DTrackFitterALT1::fit_status_t DTrackFitterALT1::LeastSquaresB(hitsInfo &hinfo, DReferenceTrajectory *rt)
{
	/// Fit the track with starting parameters given in the first step
	/// of the reference trajectory rt. On return, the reference
	/// trajectory rt will represent the final fit parameters and
	/// chisq, Ndof, resiv, cov_meas, cov_muls, and cov_parm will be
	/// filled based on the fit results.
	///
	/// This determines the best fit of the track using the least squares method
	/// described by R. Mankel Rep. Prog. Phys. 67 (2004) 553-622 pg 565.
	/// Since it uses a linear approximation for the chisq dependance on
	/// the fit parameters, several calls may be required for convergence.

	// Initialize the chisq and Ndof data members in case we need to bail early
	this->chisq = 1.0E6;
	this->Ndof = 0;

	// Make sure RT is not empty
	if(rt->Nswim_steps<1)return kFitFailed;

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

	// Get the chi-squared vector for the initial reference trajectory
	double initial_chisq;
	int initial_Ndof;
	vector<resiInfo> tmpresiduals;
	vector<bool> good_initial = GetResiInfo(rt, hinfo, tmpresiduals);
	double initial_chisq_per_dof = ChiSq(tmpresiduals, &initial_chisq, &initial_Ndof);
	DMatrix resiv_initial(resiv);
	DMatrix cov_meas_initial(cov_meas);
	DMatrix cov_muls_initial(cov_muls);
	DMatrix weights_initial(weights);
	
	// Check that the initial chisq is reasonable before continuing
	if(initial_Ndof<1){
		if(DEBUG_LEVEL>0)_DBG_<<"Initial Ndof<1. Aborting fit."<<endl;
		return kFitFailed;
	}

	// Because we have a complicated non-linear system, we take the derivatives
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
	// to preserve the original. Set the charge here. The rest of the
	// parameters (starting position and momentum) will be set using
	// values from the state vector.
	tmprt->q = rt->q;
	
	// dpx : tweak by 0.0001
	DMatrix state_dpx = state;
	state_dpx[state_px][0] += LEAST_SQUARES_DP;
	deltas[state_px] = state_dpx[state_px][0] - state[state_px][0];
	vector<bool> good_px = GetResiInfo(state_dpx, &start_step, tmprt, hinfo, tmpresiduals);
	double chisq_dpx = ChiSq(tmpresiduals);
	DMatrix resiv_dpx_hi(resiv);
	DMatrix &resiv_dpx_lo = resiv_initial;

	// dpy : tweak by 0.0001
	DMatrix state_dpy = state;
	state_dpy[state_py][0] += LEAST_SQUARES_DP;
	deltas[state_py] = state_dpy[state_py][0] - state[state_py][0];
	vector<bool> good_py = GetResiInfo(state_dpy, &start_step, tmprt, hinfo, tmpresiduals);
	double chisq_dpy = ChiSq(tmpresiduals);
	DMatrix resiv_dpy_hi(resiv);
	DMatrix &resiv_dpy_lo = resiv_initial;

	// dpz : tweak by 0.0001
	DMatrix state_dpz = state;
	state_dpz[state_pz][0] += LEAST_SQUARES_DP;
	deltas[state_pz] = state_dpz[state_pz][0] - state[state_pz][0];
	vector<bool> good_pz = GetResiInfo(state_dpz, &start_step, tmprt, hinfo, tmpresiduals);
	double chisq_dpz = ChiSq(tmpresiduals);
	DMatrix resiv_dpz_hi(resiv);
	DMatrix &resiv_dpz_lo = resiv_initial;

	// dv : tweak by 0.01
	DMatrix state_dv = state;
	state_dv[state_v][0] += LEAST_SQUARES_DX;
	deltas[state_v] = state_dv[state_v][0] - state[state_v][0];
	vector<bool> good_v = GetResiInfo(state_dv, &start_step, tmprt, hinfo, tmpresiduals);
	double chisq_dv = ChiSq(tmpresiduals);
	DMatrix resiv_dv_hi(resiv);
	DMatrix &resiv_dv_lo = resiv_initial;

	// dx : tweak by 0.01
	DMatrix state_dx = state;
	state_dx[state_x][0] += LEAST_SQUARES_DX;
	deltas[state_x] = state_dx[state_x][0] - state[state_x][0];
	vector<bool> good_x = GetResiInfo(state_dx, &start_step, tmprt, hinfo, tmpresiduals);
	double chisq_dx = ChiSq(tmpresiduals);
	DMatrix resiv_dx_hi(resiv);
	DMatrix &resiv_dx_lo = resiv_initial;

	// Verify all of the "good" vectors are of the same length
	unsigned int size_good = good_initial.size();
	if(   (good_px.size() != size_good)
		|| (good_py.size() != size_good)
		|| (good_pz.size() != size_good)
		|| (good_v.size()  != size_good)
		|| (good_x.size()  != size_good)){
			_DBG_<<"Size of \"good\" vectors don't match!  size_good="<<size_good<<endl;
			_DBG_<<"good_px.size()="<<good_px.size()<<endl;
			_DBG_<<"good_py.size()="<<good_py.size()<<endl;
			_DBG_<<"good_pz.size()="<<good_pz.size()<<endl;
			_DBG_<<" good_v.size()="<<good_v.size()<<endl;
			_DBG_<<" good_x.size()="<<good_x.size()<<endl;
			return kFitFailed;
	}
	
	// We need to get a list of hits that are good for all of the tweaked
	// tracks as well as the initial track.
	unsigned int Ngood = 0;
	vector<bool> good_all;
	for(unsigned int i=0; i<size_good; i++){
		bool isgood = true;
		isgood &= good_initial[i];
		isgood &= good_px[i];
		isgood &= good_py[i];
		isgood &= good_pz[i];
		isgood &= good_v[i];
		isgood &= good_x[i];
		good_all.push_back(isgood);
		if(isgood)Ngood++;
	}
	
	// Make sure there is a minimum number of hits before attempting a fit
	if(Ngood<LEAST_SQUARES_MIN_HITS){
		if(DEBUG_LEVEL>0)_DBG_<<"Not enough good hits to do fit. Aborting..."<<endl;
		return kFitFailed;
	}
	
	// Use the good_all map to filter out hits for each residual vector
	// for which somebody else did not have a good hit.
	FilterGood(resiv_initial, good_initial, good_all);
	FilterGood(resiv_dpx_hi , good_px, good_all);
	FilterGood(resiv_dpy_hi , good_py, good_all);
	FilterGood(resiv_dpz_hi , good_pz, good_all);
	FilterGood(resiv_dv_hi  , good_v , good_all);
	FilterGood(resiv_dx_hi  , good_x , good_all);
	FilterGood(weights_initial  , good_x , good_all);
	
	// Here, we check that the the residual vectors are of the same
	// dimension for all tweaked tracks and the initial track so that we
	// can proceed with building the fit matrices below.
	int Nhits = resiv_initial.GetNrows();
	if(   (resiv_dpx_hi.GetNrows() != Nhits)
		|| (resiv_dpy_hi.GetNrows() != Nhits)
		|| (resiv_dpz_hi.GetNrows() != Nhits)
		|| (resiv_dx_hi.GetNrows()  != Nhits)
		|| (resiv_dv_hi.GetNrows()  != Nhits)){
			_DBG_<<"Size of residual vectors don't match!  Nhits="<<Nhits<<endl;
			_DBG_<<"resiv_dpx_hi.GetNrows()="<<resiv_dpx_hi.GetNrows()<<endl;
			_DBG_<<"resiv_dpy_hi.GetNrows()="<<resiv_dpy_hi.GetNrows()<<endl;
			_DBG_<<"resiv_dpz_hi.GetNrows()="<<resiv_dpz_hi.GetNrows()<<endl;
			_DBG_<<" resiv_dx_hi.GetNrows()="<<resiv_dx_hi.GetNrows()<<endl;
			_DBG_<<" resiv_dv_hi.GetNrows()="<<resiv_dv_hi.GetNrows()<<endl;
			return kFitFailed;
	}

	// Print some debug messages
	if(DEBUG_LEVEL>4){
		_DBG_<<"initial_chisq_per_dof="<<initial_chisq_per_dof<<endl;
		_DBG_<<"chisq_dpx="<<chisq_dpx<<endl;
		_DBG_<<"chisq_dpy="<<chisq_dpy<<endl;
		_DBG_<<"chisq_dpz="<<chisq_dpz<<endl;
		_DBG_<<"chisq_dv="<<chisq_dv<<endl;
		_DBG_<<"chisq_dx="<<chisq_dx<<endl;
		if(DEBUG_LEVEL>10){
			_DBG_<<"hit\tinitial\tpx   \tpy   \tpz   \tx   \tv"<<endl;
			for(int j=0; j<resiv_initial.GetNrows(); j++){
				cout<<j<<"\t";
				cout<<resiv_initial[j][0]/sqrt(cov_meas[j][j])<<"\t";
				cout<<resiv_dpx_hi[j][0]/sqrt(cov_meas[j][j])<<"\t";
				cout<<resiv_dpy_hi[j][0]/sqrt(cov_meas[j][j])<<"\t";
				cout<<resiv_dpz_hi[j][0]/sqrt(cov_meas[j][j])<<"\t";
				cout<<resiv_dv_hi[j][0]/sqrt(cov_meas[j][j])<<"\t";
				cout<<resiv_dx_hi[j][0]/sqrt(cov_meas[j][j])<<endl;
			}
		}
	}

	// Build "F" matrix of derivatives
	DMatrix F(Nhits,Nparameters);
	for(int i=0; i<Nhits; i++){
		switch(Nparameters){
			// Note: This is a funny way to use a switch!
			case 5: F[i][state_v ] = (resiv_dv_hi[i][0]-resiv_dv_lo[i][0])/deltas[state_v];
			case 4: F[i][state_x ] = (resiv_dx_hi[i][0]-resiv_dx_lo[i][0])/deltas[state_x];
			case 3: F[i][state_pz] = (resiv_dpz_hi[i][0]-resiv_dpz_lo[i][0])/deltas[state_pz];
			case 2: F[i][state_py] = (resiv_dpy_hi[i][0]-resiv_dpy_lo[i][0])/deltas[state_py];
			case 1: F[i][state_px] = (resiv_dpx_hi[i][0]-resiv_dpx_lo[i][0])/deltas[state_px];
		}
	}
	
	// Sometimes, "F" has lots of values like 1.44E+09 indicating a problem (I think comming
	// from some nan values floating around.) Anyway, in these cases, the E2Norm value is
	// quite large (>1E+18) which we  can use to punt now. In reality, we do this to avoid
	// ROOT error messages about a matrix being singular when Ft*Vinv*F is inverted below.
	if(F.E2Norm()>1.0E18){
		if(DEBUG_LEVEL>1){
			_DBG_<<" -- F matrix E2Norm out of range(E2Norm="<<F.E2Norm()<<" max="<<1.0E18<<")"<<endl;
		}
		return kFitFailed;
	} 
	
	// Transpose of "F" matrix
	DMatrix Ft(DMatrix::kTransposed, F);

	// Calculate the "B" matrix using the weights from the initial chisq
	DMatrix &Vinv = weights_initial; // Stick with Mankel naming convention
	DMatrix B(DMatrix::kInverted, Ft*Vinv*F);
	
	// If the inversion failed altogether then the invalid flag
	// will be set on the matrix. In these cases, we're dead.
	if(!B.IsValid()){
		if(DEBUG_LEVEL>1)_DBG_<<" -- B matrix invalid"<<endl;
		return kFitFailed;
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
		if(DEBUG_LEVEL>1){
			_DBG_<<" -- B matrix E2Norm out of range(E2Norm="<<B.E2Norm()<<" max="<<LEAST_SQUARES_MAX_E2NORM<<")"<<endl;
			cout<<"--- B ---"<<endl;
			B.Print();
		}
		return kFitFailed;
	}

	// Copy the B matrix into cov_parm to later copy into DTrack
	cov_parm.ResizeTo(B);
	cov_parm = B;

	// Calculate step direction and magnitude	
	DMatrix delta_state = B*Ft*Vinv*resiv_initial;

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
	double min_chisq_per_dof = initial_chisq_per_dof;
	double min_lambda = 0.0;
	double lambda = -1.0;
	int Ntrys = 0;
	int max_trys = 6;
	DMatrix new_state(5,1);
	for(; Ntrys<max_trys; Ntrys++){

		// Scale the delta by lambda to take a partial step (except the 1st iteration where lambda is 1)
		for(int i=0; i<Nparameters; i++)new_state[i][0] = state[i][0] + delta_state[i][0]*lambda;

		GetResiInfo(new_state, &start_step, tmprt, hinfo, tmpresiduals);
		double chisq_per_dof = ChiSq(tmpresiduals);
		
		if(chisq_per_dof<min_chisq_per_dof){
			min_chisq_per_dof = chisq_per_dof;
			min_lambda = lambda;
		}
		
		// If we're at a lower chi-sq then we're done
		if(DEBUG_LEVEL>4)_DBG_<<" -- initial_chisq_per_dof="<<initial_chisq_per_dof<<"  new chisq_per_dof="<<chisq_per_dof<<" nhits="<<resiv.GetNrows()<<" p="<<tmprt->swim_steps[0].mom.Mag()<<"  lambda="<<lambda<<endl;
		if(chisq_per_dof-initial_chisq_per_dof < 0.1 && chisq_per_dof<2.0)break;

		if(chisq_per_dof<initial_chisq_per_dof)break;

		// Chi-sq was increased, try a smaller step on the next iteration
		lambda/=2.0;
	}
	
	// If we failed to find a better Chi-Sq above, maybe we were looking 
	// in the wrong direction(??) Try looking in the opposite direction.
	//if(Ntrys>=max_trys && (min_chisq_per_dof>=initial_chisq_per_dof || min_chisq_per_dof>1.0)){
	if(Ntrys>=max_trys){
		lambda = 1.0/4.0;
		for(int j=0; j<3; j++, Ntrys++){

			// Scale the delta by lambda to take a partial step (except the 1st iteration where lambda is 1)
			for(int i=0; i<Nparameters; i++)new_state[i][0] = state[i][0] + delta_state[i][0]*lambda;

			GetResiInfo(new_state, &start_step, tmprt, hinfo, tmpresiduals);
			double chisq_per_dof = ChiSq(tmpresiduals);
		
			if(chisq_per_dof<min_chisq_per_dof){
				min_chisq_per_dof = chisq_per_dof;
				min_lambda = lambda;
			}
			
			// If we're at a lower chi-sq then we're done
			if(DEBUG_LEVEL>4)_DBG_<<" -- initial_chisq_per_dof="<<initial_chisq_per_dof<<"  new chisq_per_dof="<<chisq_per_dof<<" nhits="<<resiv.GetNrows()<<" p="<<tmprt->swim_steps[0].mom.Mag()<<"  lambda="<<lambda<<endl;
			if(chisq_per_dof-initial_chisq_per_dof < 0.1 && chisq_per_dof<2.0)break;

			if(chisq_per_dof<initial_chisq_per_dof)break;
			
			// Chi-sq was increased, try a smaller step on the next iteration
			lambda/=2.0;
		}
	}

	// If we failed to make a step to a smaller chi-sq then signal
	// that we were unable to make any improvement.
	if(min_lambda==0.0){
		if(DEBUG_LEVEL>1)_DBG_<<"Chisq only increased (both directions searched!)"<<endl;
		GetResiInfo(rt, hinfo, tmpresiduals);
		ChiSq(tmpresiduals, &this->chisq, &this->Ndof); // refill resiv, cov_meas, ...
		return kFitNoImprovement;
	}
	
	// Re-create new_state using min_lambda
	for(int i=0; i<Nparameters; i++)new_state[i][0] = state[i][0] + delta_state[i][0]*min_lambda;
	if(DEBUG_HISTS)this->lambda->Fill(min_lambda);

	// Re-swim reference trajectory using these parameters and re-calc chisq.
	// Note that here we have the chisq and Ndof members set.
	GetResiInfo(new_state, &start_step, rt, hinfo, tmpresiduals);
	ChiSq(tmpresiduals, &this->chisq, &this->Ndof);
	
	if(DEBUG_LEVEL>3){
		DVector3 pos = start_step.origin;
		DVector3 mom = start_step.mom;
		double phi = mom.Phi();
		if(phi<0.0)phi+=2.0*M_PI;
		_DBG_<<"LeastSquaresB succeeded: p="<<mom.Mag()<<" theta="<<mom.Theta()<<" phi="<<phi<<" z="<<pos.Z()<<endl;
	}

	return kFitSuccess;
}

//------------------
// FilterGood
//------------------
void DTrackFitterALT1::FilterGood(DMatrix &my_resiv, vector<bool> &my_good, vector<bool> &good_all)
{
	/// Remove elements from my_resiv that are not "good" according to the good_all list.
	/// The my_good and good_all vectors should be the same size. The number of "true"
	/// entries in my_good should be the size of the my_resiv vector. For entries in the
	/// my_good vector that are true, but have an entry in good_all that is false, the
	/// corresponding entry in my_resiv will be removed. The my_resiv matrix (which should be
	/// N x 1) will be resized upon exit. The my_good vector will be set equal to good_all
	/// upon exit also so that my_resiv and my_good stay in sync.
	
	// Make list of rows (and columns) we should keep
	vector<int> rows_to_keep;
	for(unsigned int i=0, n=0; i<my_good.size(); i++){
		if(my_good[i] && good_all[i])rows_to_keep.push_back(n);
		if(my_good[i])n++;
	}
	
	// Copy my_resiv to a temporary matrix and resize my_resiv to the new size
	int Nrows = (int)rows_to_keep.size();
	int Ncols = my_resiv.GetNcols()>1 ? Nrows:1;
	DMatrix tmp(my_resiv);
	my_resiv.ResizeTo(Nrows, Ncols);
	
	// Loop over rows and columns copying in the elements we're keeping
	for(int i=0; i<Nrows; i++){
		int irow = rows_to_keep[i];
		for(int j=0; j<Ncols; j++){
			int icol = Ncols>1 ? rows_to_keep[j]:0;
			my_resiv[i][j] = tmp[irow][icol];
		}
	}

	my_good = good_all;
}

//------------------
// PrintChisqElements
//------------------
void DTrackFitterALT1::PrintChisqElements(DMatrix &resiv, DMatrix &cov_meas, DMatrix &cov_muls, DMatrix &weights)
{
	/// This is for debugging only.
	int Nhits = resiv.GetNrows();
	double chisq_diagonal = 0.0;
	for(int i=0; i<Nhits; i++){
		_DBG_<<" r/sigma "<<i<<": "<<resiv[i][0]*sqrt(weights[i][i])
				<<"  resi="<<resiv[i][0]
				<<" sigma="<<1.0/sqrt(weights[i][i])
				<<" cov_meas="<<cov_meas[i][i]
				<<" cov_muls="<<cov_muls[i][i]
				<<endl;
		chisq_diagonal += pow(resiv[i][0], 2.0)*weights[i][i];
	}
	DMatrix resiv_t(DMatrix::kTransposed, resiv);
	DMatrix chisqM(resiv_t*weights*resiv);
	int Ndof = Nhits - 5; // assume 5 fit parameters
	
	_DBG_<<" chisq/Ndof: "<<chisqM[0][0]/(double)Ndof<<"  chisq/Ndof diagonal elements only:"<<chisq_diagonal/(double)Ndof<<endl;
}

//------------------
// ForceLRTruth
//------------------
void DTrackFitterALT1::ForceLRTruth(JEventLoop *loop, DReferenceTrajectory *rt, hitsInfo &hinfo)
{
	/// This routine is called when the TRKFIT:LR_FORCE_TRUTH parameters is
	/// set to a non-zero value (e.g. -PTRKFIT:LR_FORCE_TRUTH=1 is passed
	/// on the command line). This is used only for debugging and only with
	/// Monte Carlo data.
	///
	/// This routine will adjust the left-right choice of each hit based
	/// on the current fit track (represented by rt), and the truth information
	/// contained in the the DMCTruthPoint objects. It assumes the hits in
	/// hinfo correspond to the track represented by rt.
	
	// Get Truth hits
	vector<const DMCTrackHit*> mctrackhits;
	loop->Get(mctrackhits);
	
	// Loop over hits
	for(unsigned int i=0; i<hinfo.size(); i++){
		hitInfo &hi = hinfo[i];
		const DCoordinateSystem *wire = hi.wire;
		if(wire==target)continue; // ignore target
		
		// Sometimes dist is NaN
		if(!finite(hi.dist)){
			hi.err = 100.0;
			hi.u_err = 0.0;
			continue;
		}
		
		// Find the truth hit corresponding to this real hit
		const DMCTrackHit *mctrackhit = DTrackHitSelectorTHROWN::GetMCTrackHit(hi.wire, fabs(hi.dist), mctrackhits, 0);
		
		// If no truth hit was found, then this may be a noise hit. Set
		// the error to something large so it is ignored and warn the
		// user.
		if(!mctrackhit){
			if(DEBUG_LEVEL>1)_DBG_<<"No DMCTrackHit found corresponding to hit "<<i<<" in hinfo! (noise hit?)"<<endl;
			hi.err = 100.0;
			hi.u_err = 0.0;
			continue;
		}
		
		// We do this by looking at the direction of the vector pointing from the
		// DOCA point on the wire to that of the fit track and comparing it to
		// a similar vector pointing to the truth point. If they both are generally
		// in the same direction (small angle) then the fit track is considered
		// as being on the correct side of the wire. If they have an angle somewhere
		// in the 180 degree range, the fit is considered to be on the wrong side
		// of the wire and it is "flipped".
		//
		// It can also be that the truth vector and the fit vector are at nearly
		// a right angle with respect to one another. This really only happens in
		// the CDC for tracks going in roughly the same direction as the wire.
		// In these cases, the truth point can't help us solve the LR so we set
		// the drift distance to zero and give it a larger error so this hit will
		// still be included, but appropriately weighted.
						
		// Vector pointing from wire to the fit track
		DVector3 pos_doca, mom_doca;
		rt->DistToRT(wire);
		rt->GetLastDOCAPoint(pos_doca, mom_doca);
		DVector3 shift = wire->udir.Cross(mom_doca);
		double u = rt->GetLastDistAlongWire();
		DVector3 pos_wire = wire->origin + u*wire->udir;

		// Vector pointing from wire to the truth point
		DVector3 pos_truth(mctrackhit->r*cos(mctrackhit->phi), mctrackhit->r*sin(mctrackhit->phi), mctrackhit->z);
		DVector3 pos_wire_truth = wire->origin + (pos_truth - wire->origin).Dot(wire->udir)*wire->udir;

		if(DEBUG_LEVEL>8)_DBG_<<" "<<i+1<<": (pos_doca-pos_wire).Angle(pos_truth-pos_wire_truth)="<<(pos_doca-pos_wire).Angle(pos_truth-pos_wire_truth)*57.3<<endl;

		// Decide what to do based on angle between fit and truth vectors
		double angle_fit_truth = (pos_doca-pos_wire).Angle(pos_truth-pos_wire_truth);
		if(fabs( fabs(angle_fit_truth)-M_PI/2.0) < M_PI/3.0){
			// Fit and truth vector are nearly at 90 degrees. Fit this hit
			// only to wire position
			if(DEBUG_LEVEL>5)_DBG_<<"Downgrading "<<i+1<<"th hit to wire-based (hi.u_err was:"<<hi.u_err<<")"<<endl;
			hi.dist = 0.0;
			hi.err = 0.35;
			hi.u_err = 0.0;
		}else{
			// Fit and truth vector are nearly (anti)parallel. Decide whether to "flip" hit
			if(fabs(angle_fit_truth) > M_PI/2.0){
				if(DEBUG_LEVEL>5)_DBG_<<"Flipping side "<<i+1<<"th hit "<<endl;
				hi.dist = -hi.dist;
			}
		}
	}
}

//------------------
// FillDebugHists
//------------------
void DTrackFitterALT1::FillDebugHists(DReferenceTrajectory *rt, DVector3 &vertex_pos, DVector3 &vertex_mom)
{
	//vertex_mom.SetMagThetaPhi(6.0, 17.2*M_PI/180.0, 90.0*M_PI/180.0);
	//vertex_pos.SetXYZ(0.0,0.0,65.0);
	//rt->Swim(vertex_pos, vertex_mom);
	ptotal->Fill(vertex_mom.Mag());

	// Calculate particle beta
	double beta = 1.0/sqrt(1.0+pow(rt->GetMass(), 2.0)/vertex_mom.Mag2()); // assume this is a pion for now. This should eventually come from outer detectors

	for(unsigned int j=0; j<cdchits.size(); j++){
		const DCDCTrackHit *hit = cdchits[j];
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
	
	for(unsigned int j=0; j<fdchits.size(); j++){
		const DFDCPseudo *hit = fdchits[j];
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
}


