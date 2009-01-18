// $Id$
//
//    File: DTrack_factory_THROWN.cc
// Created: Mon Sep  3 19:57:11 EDT 2007
// Creator: davidl (on Darwin Amelia.local 8.10.1 i386)
//

#include <cmath>
using namespace std;

#include <DANA/DApplication.h>

#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCPseudo.h>

#include "DTrack_factory_THROWN.h"
#include "DMCThrown.h"
#include "DReferenceTrajectory.h"
#include "DRandom.h"
#include "DMatrix.h"


//------------------
// evnt
//------------------
jerror_t DTrack_factory_THROWN::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DMCThrown*> mcthrowns;
	vector<const DCDCTrackHit*> cdctrackhits;
	vector<const DFDCPseudo*> fdcpseudos;
	loop->Get(mcthrowns);
	loop->Get(cdctrackhits);
	loop->Get(fdcpseudos);

	for(unsigned int i=0; i< mcthrowns.size(); i++){
		const DMCThrown *thrown = mcthrowns[i];
		const DKinematicData *kd_thrown = thrown;

		if(fabs(thrown->charge())<1)continue;

		// First, copy over the DKinematicData part
		DTrack *track = new DTrack;
		DKinematicData *kd_track = track;
		*kd_track = *kd_thrown;

		// Create and fill the covariance matrix for the track.
		// We need to fill this using errors estimated from the thrown
		// momentum and angle. 
		//DMatrixDSym errMatrix(1,7);
		DMatrixDSym errMatrix(7);
		track->setErrorMatrix(errMatrix);

		// Fill in DTrack specific members. (Some of these are redundant)
		DVector3 pos = track->position();
		DVector3 mom = track->momentum();
		track->candidateid = 0;
		track->chisq = 0.0;
		track->Ndof = 0.0;

		// Adapted from Alex's fortran code to paramatrize the resolution 
		// of GlueX. He developed this from a study Dave Lawrence did. 
		//this->SmearMomentum(track);

		// We need to swim a reference trajectory here. To avoid the overhead
		// of allocating/deallocating them every event, we keep a pool and
		// re-use them. If the pool is not big enough, then add one to the
		// pool.
		if(rt.size()<=_data.size()){
			// This is a little ugly, but only gets called a few times throughout the life of the process
			// Note: these never get deleted, even at the end of process.
			DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
			rt.push_back(new DReferenceTrajectory(dapp->GetBfield()));
		}

		rt[_data.size()]->Swim(pos, mom, track->charge());
		track->rt = rt[_data.size()];
		
		// Find hits that look like they belong to this track and add them in as associated objects.
		// This is a bit dangerous since "real" DTrack factories use the associated objects to
		// convey the hits actually used in the fit. The is no fit here so to be consistent, we
		// shouldn't add these. However, other factories that use DTrack (DParticle) rely on the
		// hit list to do their fit so we must add thme here.
		//
		// There is a good probability that this will fall out of alignment with the selection
		// criteria used in DTrack_factory. We should probably have this factory derive from
		// DTrack_factory just so we can guarantee they stay in sync. However, that would take
		// time that I don't have right now.
		cdchits.clear();
		fdchits.clear();
		AddCDCTrackHits(rt[_data.size()], cdctrackhits);
		AddFDCPseudoHits(rt[_data.size()], fdcpseudos);
		for(unsigned int i=0; i<cdchits.size(); i++)track->AddAssociatedObject(cdchits[i]);
		for(unsigned int i=0; i<fdchits.size(); i++)track->AddAssociatedObject(fdchits[i]);

		_data.push_back(track);
	}

	return NOERROR;
}


//==============================================================
// NOTE: The following routines were copied from DTrack_factory
//==============================================================


//------------------
// AddCDCTrackHits
//------------------
void DTrack_factory_THROWN::AddCDCTrackHits(DReferenceTrajectory *rt, vector<const DCDCTrackHit*> &cdctrackhits)
{
	/// Determine the probability that for each CDC hit that it came from the particle with the given trajectory.
	///
	/// This will calculate a probability for each CDC hit that
	/// it came from the particle represented by the given
	/// DReference trajectory. The probability is based on
	/// the residual between the distance of closest approach
	/// of the trajectory to the wire and the drift time.

	// Calculate beta of particle assuming its a pion for now. If the
	// particles is really a proton or an electron, the residual
	// calculated below will only be off by a little.
	double TOF_MASS = 0.13957018;
	double beta = 1.0/sqrt(1.0+TOF_MASS*TOF_MASS/rt->swim_steps[0].mom.Mag2());
	
	// The error on the residual. This should include the
	// error from measurement,particle parameters, and multiple 
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
		if(probability>=MIN_HIT_PROB)cdchits.push_back(hit);

		if(debug_level>10)_DBG_<<"s="<<s<<" doca="<<doca<<" dist="<<dist<<" resi="<<resi<<" tof="<<tof<<" prob="<<probability<<endl;
	}
}

//------------------
// AddFDCPseudoHits
//------------------
void DTrack_factory_THROWN::AddFDCPseudoHits(DReferenceTrajectory *rt, vector<const DFDCPseudo*> &fdcpseudos)
{
	/// Determine the probability that for each FDC hit that it came from the particle with the given trajectory.
	///
	/// This will calculate a probability for each FDC hit that
	/// it came from the particle represented by the given
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
	// error from measurement,particle parameters, and multiple 
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

		if(probability>=MIN_HIT_PROB)fdchits.push_back(hit);

		if(debug_level>10)_DBG_<<"s="<<s<<" doca="<<doca<<" dist="<<dist<<" resi="<<resi<<" tof="<<tof<<" prob="<<probability<<endl;
	}
}
