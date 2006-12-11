// $Id$
//
//    File: DTrack_factory.cc
// Created: Sun Apr  3 12:28:45 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include <math.h>

#include <TVector3.h>
#include <TMatrixD.h>

#include <JANA/JEventLoop.h>

#include "GlueX.h"
#include "DANA/DApplication.h"
#include "DMagneticFieldStepper.h"
#include "DTrackCandidate.h"
#include "DTrack_factory.h"
#include "CDC/DCDCTrackHit.h"
#include "FDC/DFDCPseudo.h"
#include "DReferenceTrajectory.h"

#define NaN std::numeric_limits<double>::quiet_NaN()

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


//------------------
// init
//------------------
jerror_t DTrack_factory::init(void)
{
	max_swim_steps = 5000;
	swim_steps = new DReferenceTrajectory::swim_step_t[max_swim_steps];

	max_swim_steps_ls = max_swim_steps;
	swim_steps_ls = new DReferenceTrajectory::swim_step_t[max_swim_steps];

	MAX_HIT_DIST = 10.0; // cm
	DEBUG_HISTS = false;
	USE_CDC = true;
	USE_FDC_ANODE = true;
	USE_FDC_CATHODE = true;
	
	gPARMS->SetDefaultParameter("TRK:MAX_HIT_DIST",	MAX_HIT_DIST);
	gPARMS->SetDefaultParameter("TRK:DEBUG_HISTS",	DEBUG_HISTS);
	gPARMS->SetDefaultParameter("TRK:USE_CDC",			USE_CDC);
	gPARMS->SetDefaultParameter("TRK:USE_FDC_ANODE",	USE_FDC_ANODE);
	gPARMS->SetDefaultParameter("TRK:USE_FDC_CATHODE",USE_FDC_CATHODE);
		
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrack_factory::brun(JEventLoop *loop, int runnumber)
{

	//gPARMS->SetParameter("GEOM:BZ_CONST",  -2.0);	
	bfield = new DMagneticFieldMap(); // temporary until new geometry scheme is worked out
		
	gPARMS->GetParameter("TRK:TRACKHIT_SOURCE",	TRACKHIT_SOURCE);
	
	CDC_Z_MIN = 17.0;
	CDC_Z_MAX = CDC_Z_MIN + 175.0;
	hit_based = false;
	//cout<<__FILE__<<":"<<__LINE__<<"-------------- Helical Fitter TRACKING --------------"<<endl;
	cout<<__FILE__<<":"<<__LINE__<<"-------------- Least Squares TRACKING --------------"<<endl;
	
	if(DEBUG_HISTS){
		cdcdoca_vs_dist = new TH2F("cdcdoca_vs_dist","DOCA vs. DIST",300, 0.0, 1.2, 300, 0.0, 1.2);
		fdcdoca_vs_dist = new TH2F("fdcdoca_vs_dist","DOCA vs. DIST",500, 0.0, 2.0, 500, 0.0, 2.0);
		cdcdoca_vs_dist_vs_ring = new TH3F("cdcdoca_vs_dist_vs_ring","DOCA vs. DIST vs. ring",300, 0.0, 1.2, 300, 0.0, 1.2,23,0.5,23.5);
		dist_axial = new TH1F("dist_axial","Distance from drift time for axial CDC wires",300,0.0,3.0);
		doca_axial = new TH1F("doca_axial","DOCA of track for axial CDC wires",300,0.0,3.0);
		dist_stereo = new TH1F("dist_stereo","Distance from drift time for stereo CDC wires",300,0.0,3.0);
		doca_stereo = new TH1F("doca_stereo","DOCA of track for axial CDC wires",300,0.0,3.0);
		chisq_final_vs_initial = new TH2F("chisq_final_vs_initial","Final vs. initial chi-sq.",200, 0.0, 10.0,50, 0.0, 2.5);
		nhits_final_vs_initial = new TH2F("nhits_final_vs_initial","Final vs. initial nhits used in chi-sq.",30, -0.5, 29.5, 30, -0.5, 29.5);
		residuals = new TH1F("residuals","Residuals",1000,0.0,2.0);
		Npasses = new TH1F("Npasses","Npasses", 21, -0.5, 20.5);
		ptotal = new TH1F("ptotal","ptotal",1000, 0.1, 8.0);
	}

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrack_factory::evnt(JEventLoop *loop, int eventnumber)
{
	// Get the track candidates and hits
	vector<const DTrackCandidate*> trackcandidates;
	loop->Get(trackcandidates);

	cdctrackhits.clear();
	fdctrackhits.clear();
	loop->Get(cdctrackhits);
	loop->Get(fdctrackhits);

	// Loop over track candidates
	for(unsigned int i=0; i<trackcandidates.size(); i++){

		// Fit the track
		DTrack *track = FitTrack(trackcandidates[i]);

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
	delete[] swim_steps;
	delete[] swim_steps_ls;

	return NOERROR;
}

//------------------
// FitTrack
//------------------
DTrack* DTrack_factory::FitTrack(const DTrackCandidate *tc)
{
	/// Fit a track candidate
	
	// Get starting position and momentum from track candidate
	TVector3 mom;
	mom.SetMagThetaPhi(tc->p, tc->theta, tc->phi);
	TVector3 pos(0.0, 0.0, tc->z_vertex);

//mom.SetMagThetaPhi(6.0,6.0*M_PI/180.0,90.0*M_PI/180.0);
//pos.SetXYZ(0.0,0.0,65.0);

	// Generate reference trajectory and use it to find the initial
	// set of hits for this track. Some of these could be dropped by
	// the fitter.
	DReferenceTrajectory *rt = new DReferenceTrajectory(bfield, tc->q, pos, mom, swim_steps, max_swim_steps,0.5);
	GetCDCTrackHits(rt); //Hits are left in private member data "cdchits_on_track"
	GetFDCTrackHits(rt); //Hits are left in private member data "fdchits_on_track"
	if((cdchits_on_track.size() + fdchits_on_track.size())<4)return NULL; // can't fit a track with less than 4 hits!

	// Fit the track and get the results in the form of the
	// position/momentum when the track was closest to the beamline
	TVector3 vertex_pos=pos; // to hold fitted values on return
	TVector3 vertex_mom=mom; // to hold fitted values on return
	double chisq=1.0E6, last_chisq;
	int Niterations;
	for(Niterations=0; Niterations<20; Niterations++){
		last_chisq = chisq;
		chisq = LeastSquares(pos, mom, rt, vertex_pos, vertex_mom);
		if(vertex_pos==pos && vertex_mom==mom)break;
		if(fabs(last_chisq-chisq) < 1.0E-6)break;
		pos = vertex_pos;
		mom = vertex_mom;
	}
	
	if(DEBUG_HISTS)Npasses->Fill(Niterations);

	if(Niterations==0)return NULL;
	
	if(DEBUG_HISTS)FillDebugHists(rt, vertex_pos, vertex_mom, tc);

	// Delete reference trajectory
	delete rt;

	// Create new DTrack object and initialize parameters with those
	// from track candidate
	DTrack *track = new DTrack;
	track->q			= tc->q;
	track->p			= vertex_mom.Mag();
	track->theta	= vertex_mom.Theta();
	track->phi		= vertex_mom.Phi();
	if(track->phi<0.0)track->phi+=2.0*M_PI;
	track->x			= vertex_pos.X();
	track->y			= vertex_pos.Y();
	track->z			= vertex_pos.Z();
	track->candidateid = tc->id;

	return track;
}

//------------------
// GetCDCTrackHits
//------------------
void DTrack_factory::GetCDCTrackHits(DReferenceTrajectory *rt, double max_hit_dist)
{
	/// Determine the distance of each CDC hit to the reference trajectory.
	/// Ones less than TRK:MAX_HIT_DIST are added to cdchits_on_track as
	/// possibly being associated with the track.
	///
	/// We use the wire position here without considering the drift time
	/// (the drift time will be used later when applying the Kalman filter).
	/// The value of TRK:MAX_HIT_DIST should at least be larger than the
	/// straw tube radius and should probably be 2 or 3 times that.
	
	// allow caller to override default value.
	if(max_hit_dist==0.0)max_hit_dist=MAX_HIT_DIST;
	
	cdchits_on_track.clear();
	for(unsigned int j=0; j<cdctrackhits.size(); j++){
		const DCDCTrackHit *hit = cdctrackhits[j];
		
		// To find the closest point in the reference trajectory, we
		// loop over all R.T. points inside the range CDC_Z_MIN to
		// CDC_Z_MAX and find the one closest to the wire. For stereo
		// wires, we adjust the X/Y position using the Z position
		// of the trajectory point. The adjustments use slopes
		// calculated at program start in DCDCTrackHit_factory.cc.
		// See note in that file for more details.
		
		// The adjustment needed for stereo layers is just a linear
		// transformation in both X and Y. According to the HDDS note,
		// the stereo layers should be thought of as though they were
		// axial straws, rotated by the stereo angle about the axis
		// that runs from the beamline through the center of the wire,
		// (perpendicular to both). The slope of the x adjustment
		// as a function of z is proportional to sin(phi) and the
		// y adjustment proportional to cos(phi) where phi is the
		// phi angle of the wire mid-plane position in lab coordinates.
		// The maximum amplitude of the adjustment is L/2*sin(alpha)
		// where L is the CDC length(175cm) and alpha is the stereo
		// angle.		
		cdc_hit_on_track_t closest_hit;
		closest_hit.swim_step = rt->FindClosestSwimStep(hit->wire);
		if(!closest_hit.swim_step)continue;
		closest_hit.cdchit = hit;
		closest_hit.dist = rt->DistToRT(hit->wire, closest_hit.swim_step, &closest_hit.s);

		// Check if this hit is close enough to the RT to be on
		// the track. This should be a function of both the distance
		// along the track and the errors of the measurement in 3D.
		// For now, we use a single distance which may be sufficient
		// for finding hits that belong to the track.
		if(!finite(closest_hit.dist))continue;
		if(closest_hit.dist>max_hit_dist)continue;
		
		cdchits_on_track.push_back(closest_hit);
	}

	// We want the hits to be ordered by the distance along the 
	// reference trajectory starting from the outside.
	sort(cdchits_on_track.begin(), cdchits_on_track.end(), CDCTrkHitSort_C);
}

//------------------
// GetFDCTrackHits
//------------------
void DTrack_factory::GetFDCTrackHits(DReferenceTrajectory *rt, double max_hit_dist)
{
	/// Determine the distance of each CDC hit to the reference trajectory.
	/// Ones less than TRK:MAX_HIT_DIST are added to cdchits_on_track as
	/// possibly being associated with the track.
	///
	/// We use the wire position here without considering the drift time
	/// (the drift time will be used later when applying the Kalman filter).
	/// The value of TRK:MAX_HIT_DIST should at least be larger than the
	/// straw tube radius and should probably be 2 or 3 times that.
	
	// allow caller to override default value.
	if(max_hit_dist==0.0)max_hit_dist=MAX_HIT_DIST;
	
	fdchits_on_track.clear();
	for(unsigned int j=0; j<fdctrackhits.size(); j++){
		const DFDCPseudo *hit = fdctrackhits[j];
		
		// To find the closest point in the reference trajectory, we
		// loop over all R.T. points and find the one closest to the
		// wire.
		fdc_hit_on_track_t closest_hit;
		closest_hit.swim_step = rt->FindClosestSwimStep(hit->wire);
		if(!closest_hit.swim_step)continue;
		closest_hit.fdchit = hit;
		closest_hit.dist = rt->DistToRT(hit->wire, closest_hit.swim_step, &closest_hit.s);

		// Check if this hit is close enough to the RT to be on
		// the track. This should be a function of both the distance
		// along the track and the errors of the measurement in 3D.
		// For now, we use a single distance which may be sufficient
		// for finding hits that belong to the track.
		if(!finite(closest_hit.dist))continue;
		if(closest_hit.dist>max_hit_dist)continue;

		fdchits_on_track.push_back(closest_hit);
	}

	// We want the hits to be ordered by the distance along the 
	// reference trajectory starting from the outside.
	sort(fdchits_on_track.begin(), fdchits_on_track.end(), FDCTrkHitSort_C);
}

//------------------
// ChiSq
//------------------
double DTrack_factory::ChiSq(double q, TMatrixD &state, swim_step_t *start_step, DReferenceTrajectory *rt)
{
	TVector3 pos =   start_step->origin
						+ state[state_x ][0]*start_step->sdir
						+ state[state_y ][0]*start_step->tdir
						+ state[state_z ][0]*start_step->udir;
	TVector3 mom =   state[state_px][0]*start_step->sdir
						+ state[state_py][0]*start_step->tdir
						+ state[state_pz][0]*start_step->udir;

	if(rt)rt->Reswim(pos,mom);

	return ChiSq(q, pos, mom,rt);
}

//------------------
// ChiSq
//------------------
double DTrack_factory::ChiSq(double q, const TVector3 &pos, const TVector3 &mom, DReferenceTrajectory *rt)
{

	// Swim a reference trajectory using the state defined by
	// "state" at "start_step" if one is not provided.
	bool own_rt = false;
	if(!rt){
		rt = new DReferenceTrajectory(bfield, q, pos, mom);
		own_rt = true;
	}

	// Clear chisq and sigma vector
	chisqv.clear();
	sigmav.clear();
	
	// Add CDC hits (if any)
	for(unsigned int i=0; i<cdchits_on_track.size(); i++){
		cdc_hit_on_track_t &hit = cdchits_on_track[i];
		const DCoordinateSystem *wire = hit.cdchit->wire;
		
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
			double beta = 0.8; // use average beta for now. This should eventually come from outer detectors
			double tof = s/(beta*3E10*1E-9);
			dist = (hit.cdchit->tdrift-tof)*22E-4;
			sigma = 0.0200; // 200 um
		}

		// NOTE: Sometimes we push nan or large values on here
		double resi = dist - doca;
		chisqv.push_back(USE_CDC ? resi:NaN);
		sigmav.push_back(sigma);
	}
	
	// Add FDC hits (if any)
	for(unsigned int i=0; i<fdchits_on_track.size(); i++){
		fdc_hit_on_track_t &hit = fdchits_on_track[i];
		const DCoordinateSystem *wire = hit.fdchit->wire;
		
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
			double beta = 0.8;
			double tof = s/(beta*3E10*1E-9);
			dist = (hit.fdchit->time-tof)*22E-4;
			sigma = 0.0200; // 200 um
		}

		// NOTE: Sometimes we push nan or large values on here
		double resi = dist - doca;
		chisqv.push_back(USE_FDC_ANODE ? resi:NaN);
		sigmav.push_back(sigma);
		
		// For the FDC we also have a measurement along the wire
		// which we include as a separate measurement
		double u = rt->GetLastDistAlongWire();
		resi = u - hit.fdchit->s;
		chisqv.push_back(USE_FDC_CATHODE ? resi:NaN);
		sigmav.push_back(0.0200); // 200 um
		
	}

	// Add "good" hits together to get the chi-squared
	double chisq = 0;
	Ngood_chisq_hits=0.0;
	for(unsigned int i=0; i<chisqv.size(); i++){
		// Check if this hit is close enough to the RT to be on
		// the track. This should be a function of both the distance
		// along the track and the errors of the measurement in 3D.
		// For now, we use a single distance which may be sufficient.
		if(!finite(chisqv[i]))continue;
		if(DEBUG_HISTS)residuals->Fill(chisqv[i]);
		if(fabs(chisqv[i]/sigmav[i])>5.0)continue;

		chisq+=pow(chisqv[i]/sigmav[i], 2.0);
		Ngood_chisq_hits += 1.0;
	}
	chisq/=Ngood_chisq_hits;

	// If we created the reference trajectory, then delete
	if(own_rt)delete rt;

	return chisq;
}

//------------------
// LeastSquares
//------------------
double DTrack_factory::LeastSquares(TVector3 &pos, TVector3 &mom, DReferenceTrajectory *rt, TVector3 &vertex_pos, TVector3 &vertex_mom)
{
	// Determine the best fit of the track using the least squares method
	// described by R. Mankel Rep. Prog. Phys. 67 (2004) 553-622 pg 565
	
	// Because we have a non-linear system, we must take the derivatives
	// numerically.
	unsigned int Nhits = cdchits_on_track.size() + 2*fdchits_on_track.size();
	swim_step_t *start_step = &rt->swim_steps[0];
	double deltas[6];

	TVector3 pos_diff = pos-start_step->origin;
	TMatrixD state(6,1);
	state[state_px	][0] = mom.Dot(start_step->sdir);
	state[state_py	][0] = mom.Dot(start_step->tdir);
	state[state_pz	][0] = mom.Dot(start_step->udir);
	state[state_x	][0] = pos_diff.Dot(start_step->sdir);
	state[state_y	][0] = pos_diff.Dot(start_step->tdir);
	state[state_z	][0] = pos_diff.Dot(start_step->udir);
	
	// Create reference trajectory to use in calculating derivatives
	DReferenceTrajectory *myrt = new DReferenceTrajectory(bfield, rt->q, pos, mom, swim_steps_ls, max_swim_steps_ls,0.5);

	// Best-guess
	ChiSq(rt->q, pos, mom, myrt);
	vector<double> resi = chisqv;
	vector<double> errs = sigmav;

	// Note that in the calculations of the deltas below,
	// the change in state should be set first and the value
	// of deltas[...] calculated from that. See Numerical
	// Recipes in C 2nd ed. section 5.7 ppg. 186-189.

	// dpx : tweak by +/- 0.01
	TMatrixD state_dpx = state;
	state_dpx[state_px][0] += 0.001;
	deltas[state_px] = state_dpx[state_px][0] - state[state_px][0];
	ChiSq(rt->q, state_dpx, start_step,myrt);
	vector<double> resi_dpx_hi = chisqv;
	vector<double> &resi_dpx_lo = resi;

	// dpy : tweak by +/- 0.01
	TMatrixD state_dpy = state;
	state_dpy[state_py][0] += 0.001;
	deltas[state_py] = state_dpx[state_px][0] - state[state_px][0];
	ChiSq(rt->q, state_dpy, start_step,myrt);
	vector<double> resi_dpy_hi = chisqv;
	vector<double> &resi_dpy_lo = resi;

	// dpz : tweak by +/- 0.01
	TMatrixD state_dpz = state;
	state_dpz[state_pz][0] += 0.001;
	deltas[state_pz] = state_dpz[state_pz][0] - state[state_pz][0];
	ChiSq(rt->q, state_dpz, start_step,myrt);
	vector<double> resi_dpz_hi = chisqv;
	vector<double> &resi_dpz_lo = resi;

	// Finished with local reference trajectory
	delete myrt;
	
	// Make a list of "clean" hits. Ones with reasonably
	// small, residuals that are not "nan" for the
	// best-guess as well as the tweaked  cases.
	vector<bool> good;
	unsigned int Ngood=0;
	for(unsigned int i=0; i<Nhits; i++){
		double res;
		double max=5.0;
		res =        resi[i]; if(!finite(res) || fabs(res)>max){good.push_back(false); continue;}
		res = resi_dpx_hi[i]; if(!finite(res) || fabs(res)>max){good.push_back(false); continue;}
		//res = resi_dpx_lo[i]; if(!finite(res) || fabs(res)>max){good.push_back(false); continue;}
		res = resi_dpy_hi[i]; if(!finite(res) || fabs(res)>max){good.push_back(false); continue;}
		//res = resi_dpy_lo[i]; if(!finite(res) || fabs(res)>max){good.push_back(false); continue;}
		res = resi_dpz_hi[i]; if(!finite(res) || fabs(res)>max){good.push_back(false); continue;}
		//res = resi_dpz_lo[i]; if(!finite(res) || fabs(res)>max){good.push_back(false); continue;}

		good.push_back(true);
		Ngood++;
	}
	if(Ngood<3){
		//cout<<__FILE__<<":"<<__LINE__<<" Bad number of good distance calculations!"<<endl;
		return 1.0E6;
	}

	// Build "F" matrix of derivatives
	TMatrixD F(Ngood,3);
	for(unsigned int i=0,j=0; j<Nhits; j++){
		if(!good[j])continue; // skip bad hits
		F[i][state_px] = (resi_dpx_hi[j]-resi_dpx_lo[j])/deltas[state_px];
		F[i][state_py] = (resi_dpy_hi[j]-resi_dpy_lo[j])/deltas[state_py];
		F[i][state_pz] = (resi_dpz_hi[j]-resi_dpz_lo[j])/deltas[state_pz];
		i++;
	}
	TMatrixD Ft(TMatrixD::kTransposed, F);

	// V is a diagonal matrix of the measurement errors. In
	// principle, we could fold in the multiple scattering
	// here, but for now, we don't. This is filled in the
	// loop below.
	TMatrixD V(Ngood,Ngood);
	V.Zero();
	
	// Measurement vector. This contains the residuals between
	// DOCAs and DISTs.
	TMatrixD m(Ngood,1);
	for(unsigned int i=0,j=0; j<Nhits; j++){
		if(!good[j])continue; // skip bad hits
		m[i][0] = -resi[j]; // drift time is already subtracted in ChiSq(...)
		V[i][i] = pow(errs[j], 2.0);
		i++;
	}
	TMatrixD Vinv(TMatrixD::kInverted, V);
	TMatrixD B(TMatrixD::kInverted, Ft*Vinv*F);

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
	if(B.E2Norm() < 1.0E6){
		TMatrixD delta_state = B*Ft*Vinv*m;
		for(int i=0; i<3; i++)state[i] += delta_state[i];
	}else{
		cout<<__FILE__<<":"<<__LINE__<<" Fit failed! (B.E2Norm()="<<B.E2Norm()<<")"<<endl;
	}

	// Calculate initial particle position/momentum. We should
	// probably swim this to the DOCA with the beamline
	vertex_pos =     start_step->origin
						+ state[state_x ][0]*start_step->sdir
						+ state[state_y ][0]*start_step->tdir;
						+ state[state_z ][0]*start_step->udir;
	vertex_mom =     state[state_px][0]*start_step->sdir
						+ state[state_py][0]*start_step->tdir
						+ state[state_pz][0]*start_step->udir;

	// Check that the chi-squared has actually decreased!
	rt->Reswim(vertex_pos, vertex_mom);
	double new_chisq = ChiSq(rt->q, vertex_pos, vertex_mom, rt);
	
	return new_chisq;
}

//------------------
// FillDebugHists
//------------------
void DTrack_factory::FillDebugHists(DReferenceTrajectory *rt, TVector3 &vertex_pos, TVector3 &vertex_mom, const DTrackCandidate* tc)
{
	//vertex_mom.SetMagThetaPhi(6.0, 17.2*M_PI/180.0, 90.0*M_PI/180.0);
	//vertex_pos.SetXYZ(0.0,0.0,65.0);
	rt->Reswim(vertex_pos, vertex_mom);
	ptotal->Fill(vertex_mom.Mag());

	for(unsigned int j=0; j<cdctrackhits.size(); j++){
		const DCDCTrackHit *hit = cdctrackhits[j];

		cdc_hit_on_track_t closest_hit;
		closest_hit.swim_step = rt->FindClosestSwimStep(hit->wire);
		if(!closest_hit.swim_step)continue;
		closest_hit.cdchit = hit;
		closest_hit.dist = rt->DistToRT(hit->wire, closest_hit.swim_step, &closest_hit.s);

		// Check if this hit is close enough to the RT to be on
		// the track. This should be a function of both the distance
		// along the track and the errors of the measurement in 3D.
		// For now, we use a single distance which may be sufficient
		// for finding hits that belong to the track.
		if(!finite(closest_hit.dist))continue;
		if(closest_hit.dist>MAX_HIT_DIST)continue;

		// Temporary debugging histos
		cdcdoca_vs_dist->Fill(hit->dist, closest_hit.dist);
		cdcdoca_vs_dist_vs_ring->Fill(hit->dist, closest_hit.dist, hit->wire->ring);
		if(hit->wire->stereo==0.0){
			dist_axial->Fill(hit->dist);
			doca_axial->Fill(closest_hit.dist);
		}else{
			dist_stereo->Fill(hit->dist);
			doca_stereo->Fill(closest_hit.dist);
		}
	}
	for(unsigned int j=0; j<fdctrackhits.size(); j++){
		const DFDCPseudo *hit = fdctrackhits[j];

		fdc_hit_on_track_t closest_hit;
		closest_hit.swim_step = rt->FindClosestSwimStep(hit->wire);
		if(!closest_hit.swim_step)continue;
		closest_hit.fdchit = hit;
		closest_hit.dist = rt->DistToRT(hit->wire, closest_hit.swim_step, &closest_hit.s);

		// Check if this hit is close enough to the RT to be on
		// the track. This should be a function of both the distance
		// along the track and the errors of the measurement in 3D.
		// For now, we use a single distance which may be sufficient
		// for finding hits that belong to the track.
		if(!finite(closest_hit.dist))continue;
		if(closest_hit.dist>MAX_HIT_DIST)continue;

		// Temporary debugging histos
		double beta = 0.8;
		double tof = closest_hit.s/(beta*3E10*1E-9);
		double dist = (hit->time-tof)*22E-4;
		fdcdoca_vs_dist->Fill(dist, closest_hit.dist);
	}
	
	// Find chi-sq for both initial and final values
	double chisq_final = ChiSq(rt->q, vertex_pos, vertex_mom, rt);
	double nhits_final = Ngood_chisq_hits;
	TVector3 mom;
	mom.SetMagThetaPhi(tc->p, tc->theta, tc->phi);
	TVector3 pos(0.0, 0.0, tc->z_vertex);
	rt->Reswim(pos, mom);
	double chisq_initial = ChiSq(rt->q, pos, mom, rt);
	double nhits_initial = Ngood_chisq_hits;
	chisq_final_vs_initial->Fill(chisq_initial,chisq_final);
	nhits_final_vs_initial->Fill(nhits_initial, nhits_final);

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
void DTrack_factory::KalmanFilter(TMatrixD &state
											, TMatrixD &P
											, DReferenceTrajectory *rt
											, TVector3 &vertex_pos
											, TVector3 &vertex_mom)
{
	/// Apply the Kalman filter to the current track starting
	/// with the given initial state and the reference trajectory.

	// State vector for Kalman filter. This is kept as a TMatrixD so
	// we can use the ROOT linear algebra package. We index the
	// elements via enum for readability.
	// We start at the last point on the reference trajectory so we can
	// swim in towards the target to get the state at the vertex.
	swim_step_t *step = &rt->swim_steps[rt->Nswim_steps-1];
	TMatrixD state(5,1);
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
	TMatrixD P(5,5);
	P.UnitMatrix();
	P[state_px ][state_px ] = pow(0.1*fabs(state[state_px][0]) + 0.010 , 2.0);
	P[state_py ][state_py ] = pow(0.1*fabs(state[state_py][0]) + 0.010 , 2.0);
	P[state_pz ][state_pz ] = pow(0.1*fabs(state[state_pz][0]) + 0.010 , 2.0);
	P[state_x  ][state_x  ] = pow(5.0 , 2.0);  // cm
	P[state_y  ][state_y  ] = pow(5.0 , 2.0);  // cm

	// The A matrix propagates the state of the particle from one
	// point to the next. This is recalculated for each.
	TMatrixD A(5, 5);
	
	// The covariance matrix representing the measurement errors.
	// For the CDC we have just one measurement "r". For the FDC
	// we have two: "r" and "w", the distance along the wire.
	TMatrixD R_cdc(1,1);
	TMatrixD R_fdc(2,2);
	
	// The Q matrix represents the "process noise" in our case,
	// this is where multiple scattering comes in. It will need
	// to be calculated at each step to include M.S. for the
	// materials traversed since the last step. To start with,
	// we will set this to zero indicating no M.S., just to keep
	// it as a place holder.
	TMatrixD Q(5,5);
	Q = 0.0*A;
	
	// K is the Kalman "gain matrix"
	TMatrixD K(5,5);
	
	// H is the matrix that converts the state vector values
	// into (predicted) measurement values.
	TMatrixD H_cdc(1,5);
	TMatrixD H_fdc(2,5);

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
		TMatrixD mystate = state;
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
		TMatrixD state_dpx = state;
		deltas[state_px] = 0.001 + 0.10*state_dpx[state_px][0];
		state_dpx[state_px][0] += deltas[state_px];
		dists[state_px] = ProjectStateBackwards(q, state_dpx, &step_prev, &step_docaRT, wire);
		if(!finite(dists[state_px]))continue; // skip this hit if there's a problem

		// dpy : tweak by 10% + 1MeV/c
		TMatrixD state_dpy = state;
		deltas[state_py] = 0.001 + 0.10*state_dpy[state_py][0];
		state_dpy[state_py][0] += deltas[state_py];
		dists[state_py] = ProjectStateBackwards(q, state_dpy, &step_prev, &step_docaRT, wire);
		if(!finite(dists[state_py]))continue; // skip this hit if there's a problem

		// dpz : tweak by 10% + 1MeV/c
		TMatrixD state_dpz = state;
		deltas[state_pz] = 0.001 + 0.10*state_dpz[state_pz][0];
		state_dpz[state_pz][0] += deltas[state_pz];
		dists[state_pz] = ProjectStateBackwards(q, state_dpz, &step_prev, &step_docaRT, wire);
		if(!finite(dists[state_pz]))continue; // skip this hit if there's a problem

		// dx : tweak by 250 microns
		TMatrixD state_dx = state;
		deltas[state_x] = 0.0250;
		state_dx[state_x][0] += deltas[state_x];
		dists[state_x] = ProjectStateBackwards(q, state_dx, &step_prev, &step_docaRT, wire);
		if(!finite(dists[state_x]))continue; // skip this hit if there's a problem

		// dy : tweak by 500 microns
		TMatrixD state_dy = state;
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

		// For Monte Carlo, dist = 22um*tdrift
		TMatrixD z_minus_h(1,1);
		double z = hit.cdchit->tdrift*22.0E-4;
		z=0.0; // Hit-based tracking is equivalent to tdrift=0
		z_minus_h[0][0] = z - dist;
		
		// Measurement error on drift distance
		R_cdc[0][0] = pow(200.0E-4,2.0); // use constant 400um measurement error for now
		
		// Hmmm... it seems like V should be a unit matrix, but then, whay have it?
		TMatrixD V(TMatrixD::kUnit, R_cdc);
		
		// W corresponds to process noise covariance. We set it to the
		// identity matrix, but it is only used to multiply Q which
		// is a NULL matrix so this doesn't really matter.
		TMatrixD W(TMatrixD::kUnit, Q);

		// Apply Kalman filter equations to include this measurement point.
		KalmanStep(mystate, P, z_minus_h, A, H_cdc, Q, R_cdc, W, V);
#endif
		state = mystate;

		// Make current step_docaRT next iterations step_prev.
		step_prev = step_docaRT;
	}
	
	// Swim particle backwards to the beamline
	TVector3 pos = step_prev.pos
						+ state[state_x ][0]*step_prev.xdir
						+ state[state_y ][0]*step_prev.ydir;
	TVector3 mom =   state[state_px][0]*step_prev.xdir
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
void DTrack_factory::KalmanStep(	TMatrixD &x,
											TMatrixD &P,
											TMatrixD &z_minus_h,
											TMatrixD &A,
											TMatrixD &H,
											TMatrixD &Q,
											TMatrixD &R,
											TMatrixD &W,
											TMatrixD &V)
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
	
	TMatrixD At(TMatrixD::kTransposed, A);
	TMatrixD Ht(TMatrixD::kTransposed, H);
	TMatrixD Vt(TMatrixD::kTransposed, V);
	TMatrixD Wt(TMatrixD::kTransposed, W);

_DBG_;	
TMatrixD HPHt = H*P*Ht;
A.Print();
P.Print();
H.Print();
HPHt.Print();

TMatrixD Atmp(TMatrixD::kUnit, A);
A=Atmp;
At=Atmp;

	P = A*P*At + W*Q*Wt;
	
	TMatrixD B(TMatrixD::kInverted, H*P*Ht + V*R*Vt);
cout<<__FILE__<<":"<<__LINE__<<" B[0][0] = "<<B[0][0]<<"  sqrt(1/B[0][0])="<<sqrt(1/B[0][0])<<endl;

	TMatrixD K = P*Ht*B;
K.Print();
	TMatrixD I(TMatrixD::kUnit, P);
x.Print();
	x = x + K*(z_minus_h);
x.Print();
	P = (I - K*H)*P;
}

//------------------
// ProjectStateBackwards
//------------------
double DTrack_factory::ProjectStateBackwards(double q, TMatrixD &state
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
		TVector3 mypos;
		stepper.Step(&mypos);
		TVector3 diff = mypos - stepRT->pos;
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
	const TVector3 &p = stepRT->mom;
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
	TVector3 pos_diff = step.pos - stepRT->pos;
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
TVector3 DTrack_factory::GetDistToRT(TVector3 &hit, swim_step_t *s2)
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
	TVector3 delta_s1 = s1->pos - hit;
	TVector3 delta_s3 = s3->pos - hit;
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
	TVector3 tmpv = s_nn->pos - s2->pos; // shift so origin is that of R.T. at s2
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
	TVector3 hit_rt(tmpv.Dot(s2->xdir), tmpv.Dot(s2->ydir), tmpv.Dot(s2->zdir));
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
	TVector3 doca_point(x,y,z);

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
	TVector3 p;
	p.SetMagThetaPhi(par[0],par[1],par[2]);
	TVector3 pos(par[3],par[4],par[5]);
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
