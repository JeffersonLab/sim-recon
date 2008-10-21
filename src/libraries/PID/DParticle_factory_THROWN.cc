// $Id$
//
//    File: DParticle_factory_THROWN.cc
// Created: Sat Oct  4 22:04:56 EDT 2008
// Creator: davidl (on Darwin Amelia.local 8.11.1 i386)
//

#include <cmath>
using namespace std;

#include <DANA/DApplication.h>

#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCPseudo.h>

#include "DParticle_factory_THROWN.h"
#include <TRACKING/DMCThrown.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <DRandom.h>
#include <DMatrix.h>

//------------------
// Constructor
//------------------
DParticle_factory_THROWN::DParticle_factory_THROWN()
{
	DEBUG_LEVEL = 0;

	// Define target center
	target = new DCoordinateSystem();
	target->origin.SetXYZ(0.0, 0.0, 65.0);
	target->sdir.SetXYZ(1.0, 0.0, 0.0);
	target->tdir.SetXYZ(0.0, 1.0, 0.0);
	target->udir.SetXYZ(0.0, 0.0, 1.0);
	target->L = 30.0;
}

//------------------
// brun
//------------------
jerror_t DParticle_factory_THROWN::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	// Get resolutions of chambers. These should be set somewhere in the TRACKING package
	SIGMA_CDC = 0.0;
	SIGMA_FDC_ANODE = 0.0;
	SIGMA_FDC_CATHODE = 0.0;
	
	JParameter *parm = gPARMS->GetParameter("TRKFIT:SIGMA_CDC");
	if(parm)SIGMA_CDC = parm->d();
	
	parm = gPARMS->GetParameter("TRKFIT:SIGMA_FDC_ANODE");
	if(parm)SIGMA_FDC_ANODE = parm->d();
	
	parm = gPARMS->GetParameter("TRKFIT:SIGMA_FDC_CATHODE");
	if(parm)SIGMA_FDC_CATHODE = parm->d();

	_DBG_<<"        SIGMA_CDC = "<<SIGMA_CDC<<endl;
	_DBG_<<"  SIGMA_FDC_ANODE = "<<SIGMA_FDC_ANODE<<endl;
	_DBG_<<"SIGMA_FDC_CATHODE = "<<SIGMA_FDC_CATHODE<<endl;

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DParticle_factory_THROWN::evnt(JEventLoop *loop, int eventnumber)
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
		DParticle *particle = new DParticle;
		DKinematicData *kd_particle = particle;
		*kd_particle = *kd_thrown;

		// Create and fill the covariance matrix for the particle.
		// We need to fill this using errors estimated from the thrown
		// momentum and angle. 
		//DMatrixDSym errMatrix(1,7);
		DMatrixDSym errMatrix(7);
		particle->setErrorMatrix(errMatrix);

		// Fill in DParticle specific members. (Some of these are redundant)
		DVector3 pos = particle->position();
		DVector3 mom = particle->momentum();
		particle->chisq = 0.0;
		particle->Ndof = 0.0;

		// We need to swim a reference trajectory here. To avoid the overhead
		// of allocating/deallocating them every event, we keep a pool and
		// re-use them. If the pool is not big enough, then add one to the
		// pool.
		if(rt_pool.size()<=_data.size()){
			// This is a little ugly, but only gets called a few times throughout the life of the process
			// Note: these never get deleted, even at the end of process.
			DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
			rt_pool.push_back(new DReferenceTrajectory(dapp->GetBfield()));
		}

		DReferenceTrajectory *rt = rt_pool[_data.size()];
		rt ->Swim(pos, mom, particle->charge());
		particle->rt = rt;
		
		// Find hits that look like they belong to this particle and add them in as associated objects.
		// This is a bit dangerous since "real" DParticle factories use the associated objects to
		// convey the hits actually used in the fit. The is no fit here so to be consistent, we
		// shouldn't add these. However, other factories that use DParticle (DParticle) rely on the
		// hit list to do their fit so we must add thme here.
		//
		// There is a good probability that this will fall out of alignment with the selection
		// criteria used in DParticle_factory. We should probably have this factory derive from
		// DParticle_factory just so we can guarantee they stay in sync. However, that would take
		// time that I don't have right now.
		cdchits.clear();
		fdchits.clear();
		AddCDCTrackHits(rt, cdctrackhits);
		AddFDCPseudoHits(rt, fdcpseudos);
		
		// Add associated objects
		for(unsigned int i=0; i<cdchits.size(); i++)particle->AddAssociatedObject(cdchits[i]);
		for(unsigned int i=0; i<fdchits.size(); i++)particle->AddAssociatedObject(fdchits[i]);
		
		// Calculate chisq and Ndof
		vector<const DCoordinateSystem*> wires;
		vector<DVector3> shifts;
		vector<double> errs;
		GetWiresShiftsErrs(kTimeBased, rt, wires, shifts, errs);
		
		vector<double> chisqv;
		double chisq;
		ChiSq(rt, wires, shifts, errs, chisqv, &chisq, &particle->Ndof);
		particle->chisq = chisq;
		
		_data.push_back(particle);
	}

	return NOERROR;
}


//==============================================================
// NOTE: The following routines were copied from DParticle_factory
//==============================================================


//------------------
// AddCDCTrackHits
//------------------
void DParticle_factory_THROWN::AddCDCTrackHits(DReferenceTrajectory *rt, vector<const DCDCTrackHit*> &cdctrackhits)
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
void DParticle_factory_THROWN::AddFDCPseudoHits(DReferenceTrajectory *rt, vector<const DFDCPseudo*> &fdcpseudos)
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
		
		// Residual
		double resi = dist - doca;		

		// Use chi-sq probaility function with Ndof=1 to calculate probability
		double probability = TMath::Prob(resi/sigma_anode, 1);

		// Cathode
		double u=rt->GetLastDistAlongWire();
		double resic = u - hit->s;

		// Same as for the anode. We multiply the
		// probabilities to get a total probability
		// based on both the anode and cathode hits.
		probability *= TMath::Prob(resic/sigma_cathode, 1);

		if(probability<=MIN_HIT_PROB)fdchits.push_back(hit);

		if(debug_level>10)_DBG_<<"s="<<s<<" doca="<<doca<<" dist="<<dist<<" resi="<<resi<<" tof="<<tof<<" prob="<<probability<<endl;
	}
}


//==============================================================
// NOTE: The following routines were copied from DParticle_factory_THROWN.cc
//==============================================================


//------------------
// ChiSq
//------------------
double DParticle_factory_THROWN::ChiSq(DReferenceTrajectory *rt, vector<const DCoordinateSystem*> &wires, vector<DVector3> &shifts, vector<double> &errs, vector<double> &chisqv, double *chisq_ptr, int *dof_ptr)
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
		for(unsigned int i=0; i<wires.size(); i++)chisqv.push_back(1.0E6);
		return 1.0E6;
	}
	
	// Sometimes, we end up with a zero step rt due to bad parameters.
	// Avoid all the warning messages by detecting this now and bailing early.
	if(rt->Nswim_steps<1){
		for(unsigned int i=0; i<wires.size(); i++)chisqv.push_back(1.0E6);
		return 1.0E6;
	}
	
	// These are for debugging
	const DFDCWire *fdcwire=NULL;
	const DCDCWire *cdcwire=NULL;
	int Npos_shifts=0;
	int Nneg_shifts=0;
	int Nzero_shifts=0;
	
	// Loop over wires
	double chisq = 0.0;
	int dof=0;
	for(unsigned int i=0; i<wires.size(); i++){
		
		if(DEBUG_LEVEL>10){
			fdcwire = dynamic_cast<const DFDCWire*>(wires[i]);
			cdcwire = dynamic_cast<const DCDCWire*>(wires[i]);
		}
	
		DCoordinateSystem wire = *wires[i];
		if(use_shifts){
			wire.origin += shifts[i];
			
			if(DEBUG_LEVEL>10){
				double a = (shifts[i].Cross(wires[i]->origin)).Z();
				if(a==0.0)Nzero_shifts++;
				else if(a<0.0)Nneg_shifts++;
				else Npos_shifts++;
			}
		}
		
		// Get distance of the (shifted) wire from the reference trajectory
		double s;
		double d = rt->DistToRT(&wire, &s);
		
		// It can happen that a track looks enough like a helix that
		// a momentum pointing in the opposite direction of the real track
		// can have a small value of d, even when using a swim step close
		// to the target. In these cases, s will be negative and large.
		// When this happens, set d to R to force a bad chisq.
		if(s<-5.0)d = wire.origin.Perp();
		
		// Add this to the chisq. It may turn out that d is NAN in which case we
		// don't want to include it in the total tally. We still must add it
		// however to the chisqv since that must have an entry for every entry
		// in the wires vector.
		//
		// The value going into chisqv needs to be a signed quantity with positive
		// values meaning the doca is larger than the dist and negative values
		// meaning the doca is smaller than the dist (doca from track, dist
		// from drift time). 
		double c = d/errs[i];
		chisqv.push_back(c);

		if(finite(c)){
			chisq += c*c;
			dof++;
			if(DEBUG_LEVEL>10){
				_DBG_<<"  chisqv["<<dof<<"] = "<<c<<"  d="<<d<<"  Rwire="<<wires[i]->origin.Perp()<<"  Rshifted="<<wire.origin.Perp()<<" s="<<s;
				if(fdcwire)cerr<<" FDC: layer="<<fdcwire->layer<<" wire="<<fdcwire->wire<<" angle="<<fdcwire->angle;
				if(cdcwire)cerr<<" CDC: ring="<<cdcwire->ring<<" straw="<<cdcwire->straw<<" stereo="<<cdcwire->stereo;
				cerr<<endl;
			}
		}else{
			if(DEBUG_LEVEL>10)_DBG_<<"  chisqv[unused] = "<<c<<"  d="<<d<<"  Rwire="<<wires[i]->origin.Perp()<<"  Rshifted="<<wire.origin.Perp()<<" s="<<s<<endl;
		}
	}

	if(DEBUG_LEVEL>10){
		_DBG_<<" Nshifts(+,-,0)=("<<Npos_shifts<<","<<Nneg_shifts<<","<<Nzero_shifts<<")"<<endl;
	}
	
	// If the caller supplied pointers to chisq and dof variables, copy the values into them
	if(chisq_ptr)*chisq_ptr = chisq;
	if(dof_ptr)*dof_ptr = dof - 5; // assume 5 fit parameters

	// If it turns out the dof is zero, return 1.0E6. Otherwise, return
	// the chisq/dof
	return dof<2 ? 1.0E6:chisq/(double)dof;
}

//------------------
// GetWiresShiftsErrs
//------------------
void DParticle_factory_THROWN::GetWiresShiftsErrs(fit_type_t fit_type, DReferenceTrajectory *rt, vector<const DCoordinateSystem*> &wires, vector<DVector3> &shifts, vector<double> &errs)
{
	swim_step_t &start_step = rt->swim_steps[0];
	DVector3 &pos = start_step.origin;

#if 1
	// --- Target ---
	// Define the target as though it were a wire so it is included
	// in the chi-sq
	target->L = 3.0;
//target->L = 0.0;
	target->origin.SetXYZ(0.0, 0.0, pos.Z());
	wires.push_back(target);
	errs.push_back(fit_type==kWireBased ? 0.1:0.1);
//errs.push_back(fit_type==kWireBased ? 0.0001:0.0001);
	if(fit_type!=kWireBased)shifts.push_back(DVector3(0.0, 0.0, 0.0));
#endif

	// --- CDC ---
	for(unsigned int i=0; i<cdchits.size(); i++){
		const DCDCTrackHit *hit = cdchits[i];
		const DCoordinateSystem *wire = hit->wire;
		wires.push_back(wire);

		// Fill in shifts and errs vectors based on whether we're doing
		// hit-based or time-based tracking
		if(fit_type==kWireBased){
			// If we're doing hit-based tracking then only the wire positions
			// are used and the drift time info is ignored. Thus, we don't
			// have to calculate the shift vectors
			//errs.push_back(0.261694); // empirically from single pi+ events
			errs.push_back(0.8/sqrt(12.0)); // variance for evenly illuminated straw
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

	// --- FDC ---
	for(unsigned int i=0; i<fdchits.size(); i++){
		const DFDCPseudo *hit = fdchits[i];
		const DCoordinateSystem *wire = hit->wire;
		wires.push_back(wire);

		// Fill in shifts and errs vectors based on whether we're doing
		// hit-based or time-based tracking
		if(fit_type==kWireBased){
			// If we're doing hit-based tracking then only the wire positions
			// are used and the drift time info is ignored. Thus, we don't
			// have to calculate the shift vectors
			//errs.push_back(0.261694); // empirically from single pi+ events
			errs.push_back(0.5/sqrt(12.0)); // variance for evenly illuminated cell
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
			double dist = hit->dist*((hit->time-tof)/hit->time);
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
			errs.push_back(SIGMA_FDC_ANODE);
		}
	}
}
