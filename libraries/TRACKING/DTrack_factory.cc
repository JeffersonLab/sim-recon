// $Id$
//
//    File: DTrack_factory.cc
// Created: Wed Sep  3 09:33:40 EDT 2008
// Creator: davidl (on Darwin harriet.jlab.org 8.11.1 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DTrack_factory.h"
#include <TRACKING/DTrackCandidate.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCPseudo.h>

using namespace jana;

//------------------
// CDCSortByRincreasing
//------------------
bool CDCSortByRincreasing(const DCDCTrackHit* const &hit1, const DCDCTrackHit* const &hit2) {
	// use the ring number to sort by R(decreasing) and then straw(increasing)
	if(hit1->wire->ring == hit2->wire->ring){
		return hit1->wire->straw < hit2->wire->straw;
	}
	return hit1->wire->ring < hit2->wire->ring;
}

//------------------
// FDCSortByZincreasing
//------------------
bool FDCSortByZincreasing(const DFDCPseudo* const &hit1, const DFDCPseudo* const &hit2) {
	// use the layer number to sort by Z(decreasing) and then wire(increasing)
	if(hit1->wire->layer == hit2->wire->layer){
		return hit1->wire->wire < hit2->wire->wire;
	}
	return hit1->wire->layer < hit2->wire->layer;
}


//------------------
// init
//------------------
jerror_t DTrack_factory::init(void)
{
	fitter = NULL;

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrack_factory::brun(jana::JEventLoop *loop, int runnumber)
{
	// Get pointer to TrackFitter object that actually fits a track
	vector<const DTrackFitter *> fitters;
	loop->Get(fitters);
	if(fitters.size()<1){
		_DBG_<<"Unable to get a DTrackFitter object! NO Charged track fitting will be done!"<<endl;
		return RESOURCE_UNAVAILABLE;
	}
	
	// Drop the const qualifier from the DTrackFitter pointer (I'm surely going to hell for this!)
	fitter = const_cast<DTrackFitter*>(fitters[0]);

	// Warn user if something happened that caused us NOT to get a fitter object pointer
	if(!fitter){
		_DBG_<<"Unable to get a DTrackFitter object! NO Charged track fitting will be done!"<<endl;
		return RESOURCE_UNAVAILABLE;
	}

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrack_factory::evnt(JEventLoop *loop, int eventnumber)
{
	if(!fitter)return NOERROR;
	
	// Get candidates and hits
	vector<const DTrackCandidate*> candidates;
	vector<const DCDCTrackHit*> cdctrackhits;
	vector<const DFDCPseudo*> fdcpseudos;
	loop->Get(candidates);
	loop->Get(cdctrackhits);
	loop->Get(fdcpseudos);
	
	// Loop over candidates
	for(unsigned int i=0; i<candidates.size(); i++){
		const DTrackCandidate *candidate = candidates[i];

		// Make sure there are enough DReferenceTrajectory objects
		while(rtv.size()<=_data.size())rtv.push_back(new DReferenceTrajectory(fitter->GetDMagneticFieldMap()));
		DReferenceTrajectory *rt = rtv[_data.size()];
		
		// Swim a reference trajectory with this candidate's parameters
		rt->Swim(candidate->position(), candidate->momentum(), candidate->charge());
		if(rt->Nswim_steps<1)continue;

#if 0
		// Get CDC and FDC hits from candidate
		vector<const DCDCTrackHit *> cdchits;
		vector<const DFDCPseudo *> fdchits;
		candidate->Get(cdchits);
		candidate->Get(fdchits);
		sort(cdchits.begin(), cdchits.end(), CDCSortByRincreasing);
		sort(fdchits.begin(), fdchits.end(), FDCSortByZincreasing);
#endif
		
		// Setup fitter to do fit
		fitter->Reset();
		AddCDCTrackHits(rt, cdctrackhits);
		AddFDCPseudoHits(rt, fdcpseudos);
		fitter->SetFitType(DTrackFitter::kWireBased);
		
		// Do the fit
		DTrackFitter::fit_status_t status = fitter->FitTrack(*candidate);
		switch(status){
			case DTrackFitter::kFitNotDone:
				_DBG_<<"Fitter returned kFitNotDone. This should never happen!!"<<endl;
				break;
			case DTrackFitter::kFitSuccess:
			case DTrackFitter::kFitNoImprovement:
				MakeDTrack(candidate);
				break;
			case DTrackFitter::kFitFailed:
				break;
		}

	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DTrack_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTrack_factory::fini(void)
{
	for(unsigned int i=0; i<rtv.size(); i++)delete rtv[i];
	rtv.clear();

	return NOERROR;
}

//------------------
// AddCDCTrackHits
//------------------
void DTrack_factory::AddCDCTrackHits(DReferenceTrajectory *rt, vector<const DCDCTrackHit*> &cdctrackhits)
{
	/// Determine the probability that for each CDC hit that it came from the track with the given trajectory.
	///
	/// This will calculate a probability for each CDC hit that
	/// it came from the track represented by the given
	/// DReference trajectory. The probability is based on
	/// the residual between the distance of closest approach
	/// of the trajectory to the wire and the drift time.

	// Calculate beta of particle assuming its a pion for now. If the
	// particles is really a proton or an electron, the residual
	// calculated below will only be off by a little.
	double TOF_MASS = 0.13957018;
	double beta = 1.0/sqrt(1.0+TOF_MASS*TOF_MASS/rt->swim_steps[0].mom.Mag2());
	
	// The error on the residual. This should include the
	// error from measurement,track parameters, and multiple 
	// scattering.
	//double sigma = sqrt(pow(SIGMA_CDC,2.0) + pow(0.4000,2.0));
	double sigma = 0.8/sqrt(12.0);
	
	// Minimum probability of hit belonging to wire and still be accepted
	double MIN_HIT_PROB = 0.2;

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
		//double resi = dist - doca;
		double resi = doca - 0.4;
		double chisq = pow(resi/sigma, 2.0);

		// Use chi-sq probaility function with Ndof=1 to calculate probability
		double probability = TMath::Prob(chisq/4.0, 1);
		if(probability>=MIN_HIT_PROB)fitter->AddHit(hit);

		if(debug_level>10)_DBG_<<"s="<<s<<" doca="<<doca<<" dist="<<dist<<" resi="<<resi<<" tof="<<tof<<" prob="<<probability<<endl;
	}
}

//------------------
// AddFDCPseudoHits
//------------------
void DTrack_factory::AddFDCPseudoHits(DReferenceTrajectory *rt, vector<const DFDCPseudo*> &fdcpseudos)
{
	/// Determine the probability that for each FDC hit that it came from the track with the given trajectory.
	///
	/// This will calculate a probability for each FDC hit that
	/// it came from the track represented by the given
	/// DReference trajectory. The probability is based on
	/// the residual between the distance of closest approach
	/// of the trajectory to the wire and the drift time
	/// and the distance along the wire.

	// Calculate beta of particle assuming its a pion for now. If the
	// particles is really a proton or an electron, the residual
	// calculated below will only be off by a little.
	double TOF_MASS = 0.13957018;
	double beta = 1.0/sqrt(1.0+TOF_MASS*TOF_MASS/rt->swim_steps[0].mom.Mag2());
	
	// The error on the residual. This should include the
	// error from measurement,track parameters, and multiple 
	// scattering.
	//double sigma_anode = sqrt(pow(SIGMA_FDC_ANODE,2.0) + pow(1.000,2.0));
	//double sigma_cathode = sqrt(pow(SIGMA_FDC_CATHODE,2.0) + pow(1.000,2.0));
	double sigma_anode = 0.5/sqrt(12.0);
	double sigma_cathode = 0.5/sqrt(12.0);
	
	// Minimum probability of hit belonging to wire and still be accepted
	double MIN_HIT_PROB = 0.2;

	for(unsigned int j=0; j<fdcpseudos.size(); j++){
		const DFDCPseudo *hit = fdcpseudos[j];
		
		// Find the DOCA to this wire
		double s;
		double doca = rt->DistToRT(hit->wire, &s);

		// Distance using drift time
		// NOTE: Right now we assume pions for the TOF
		// and a constant drift velocity of 55um/ns
		double tof = s/(beta*3E10*1E-9);
		double dist = (hit->time - tof)*55E-4;
		
		// Anode Residual
		double resi = dist - doca;		

		// Cathode Residual
		double u=rt->GetLastDistAlongWire();
		double resic = u - hit->s;

		// Probability of this hit being on the track
		double chisq = pow(resi/sigma_anode, 2.0) + pow(resic/sigma_cathode, 2.0);
		double probability = TMath::Prob(chisq/2.0, 2);

		if(probability<=MIN_HIT_PROB)fitter->AddHit(hit);

		if(debug_level>10)_DBG_<<"s="<<s<<" doca="<<doca<<" dist="<<dist<<" resi="<<resi<<" tof="<<tof<<" prob="<<probability<<endl;
	}
}

//------------------
// MakeDTrack
//------------------
void DTrack_factory::MakeDTrack(const DTrackCandidate *candidate)
{
	// Allocate a DReferenceTrajectory object if needed.
	// These each have a large enough memory footprint that
	// it causes noticable performance problems if we allocated
	// and deallocated them every event. Therefore, we allocate
	// when needed, but recycle them on the next event.
	// They are deleted in the fini method.
	while(rtv.size()<=_data.size())rtv.push_back(new DReferenceTrajectory(fitter->GetDMagneticFieldMap()));
	DReferenceTrajectory *rt = rtv[_data.size()];

	DTrack *track = new DTrack;
	
	// Copy over DKinematicData part
	DKinematicData *track_kd = track;
	*track_kd = fitter->GetFitParameters();
	rt->Swim(track->position(), track->momentum(), track->charge());
	
	track->rt = rt;
	track->chisq = fitter->GetChisq();
	track->Ndof = fitter->GetNdof();
	track->candidateid = candidate->id;
	
	// Add hits used as associated objects
	vector<const DCDCTrackHit*> cdchits = fitter->GetCDCFitHits();
	vector<const DFDCPseudo*> fdchits = fitter->GetFDCFitHits();
	sort(cdchits.begin(), cdchits.end(), CDCSortByRincreasing);
	sort(fdchits.begin(), fdchits.end(), FDCSortByZincreasing);
	for(unsigned int i=0; i<cdchits.size(); i++)track->AddAssociatedObject(cdchits[i]);
	for(unsigned int i=0; i<fdchits.size(); i++)track->AddAssociatedObject(fdchits[i]);

	// Add DTrackCandidate as associated object (yes, this is redundant with the candidateid member)
	track->AddAssociatedObject(candidate);
	
	_data.push_back(track);
}

