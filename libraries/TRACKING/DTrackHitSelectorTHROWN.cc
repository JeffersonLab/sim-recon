// $Id$
//
//    File: DTrackHitSelectorTHROWN.cc
// Created: Mon Mar  9 09:03:03 EDT 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.6.0 i386)
//

#include <TRACKING/DReferenceTrajectory.h>
#include <TRACKING/DMCThrown.h>
#include <TRACKING/DMCTrackHit.h>
#include <GlueX.h>

#include "DTrackHitSelectorTHROWN.h"


//---------------------------------
// DTrackHitSelectorTHROWN    (Constructor)
//---------------------------------
DTrackHitSelectorTHROWN::DTrackHitSelectorTHROWN(jana::JEventLoop *loop):DTrackHitSelector(loop)
{
	HS_DEBUG_LEVEL = 0;
	gPARMS->SetDefaultParameter("TRKFIT:HS_DEBUG_LEVEL", HS_DEBUG_LEVEL);
}

//---------------------------------
// ~DTrackHitSelectorTHROWN    (Destructor)
//---------------------------------
DTrackHitSelectorTHROWN::~DTrackHitSelectorTHROWN()
{

}

//---------------------------------
// GetCDCHits
//---------------------------------
void DTrackHitSelectorTHROWN::GetCDCHits(fit_type_t fit_type, DReferenceTrajectory *rt, const vector<const DCDCTrackHit*> &cdchits_in, vector<const DCDCTrackHit*> &cdchits_out) const
{
	/// Find the hits for the track indicated by the given reference trajectory using truth points.
	///
	/// This will look through the list of charged tracks in DMCThrown and find
	/// the thrown track that most closely matches the starting parameters in 
	/// "rt". Then, it will look through the DCDCTrackHit objects in cdchits_in
	/// and attempt to match them up with CDC hits in the DMCTrackHit objects.
	/// In this way, truth information can be used to generate a list of the
	/// CDC hits that belong with "rt".

	// Get the MC track number
	int track = FindTrackNumber(rt);
	if(track<1)return;
	
	// Get the DMCTrackHit objects
	vector<const DMCTrackHit*> mctrackhits;
	loop->Get(mctrackhits);

	// Here is the hard part. We need to match the hits in cdchits_in with hits
	// in DMCTrackHit objects. We do this by checking on the distance the truth
	// point is from the wire and comparing it to the drift+tof time.
	for(unsigned int i=0; i<cdchits_in.size(); i++){
		const DCDCTrackHit *hit = cdchits_in[i];

		const DMCTrackHit* mchit = GetMCTrackHit(hit->wire, hit->tdrift*55.0E-4, mctrackhits, track);
		if(mchit)cdchits_out.push_back(cdchits_in[i]);		
	}
}

//---------------------------------
// GetFDCHits
//---------------------------------
void DTrackHitSelectorTHROWN::GetFDCHits(fit_type_t fit_type, DReferenceTrajectory *rt, const vector<const DFDCPseudo*> &fdchits_in, vector<const DFDCPseudo*> &fdchits_out) const
{
	/// Find the hits for the track indicated by the given reference trajectory using truth points.
	///
	/// This will look through the list of charged tracks in DMCThrown and find
	/// the thrown track that most closely matches the starting parameters in 
	/// "rt". Then, it will look through the DFDCPseudo objects in cdchits_in
	/// and attempt to match them up with FDC hits in the DMCTrackHit objects.
	/// In this way, truth information can be used to generate a list of the
	/// FDC hits that belong with "rt".

	// Get the MC track number
	int track = FindTrackNumber(rt);
	if(track<1)return;
	
	// Get the DMCTrackHit objects
	vector<const DMCTrackHit*> mctrackhits;
	loop->Get(mctrackhits);

	// Here is the hard part. We need to match the hits in cdchits_in with hits
	// in DMCTrackHit objects. We do this by checking on the distance the truth
	// point is from the wire and comparing it to the drift+tof time.
	for(unsigned int i=0; i<fdchits_in.size(); i++){
		const DFDCPseudo *hit = fdchits_in[i];
		
		const DMCTrackHit* mchit = GetMCTrackHit(hit->wire, hit->time*55.0E-4, mctrackhits, track);
		if(mchit)fdchits_out.push_back(fdchits_in[i]);
	}
}

//---------------------------------
// FindTrackNumber
//---------------------------------
int DTrackHitSelectorTHROWN::FindTrackNumber(DReferenceTrajectory *rt) const
{
	if(rt->Nswim_steps<1)return 0;
	DVector3 &mom = rt->swim_steps[0].mom;
	
	vector<const DMCThrown*> throwns;
	loop->Get(throwns);
	
	int myid = 0;
	int min_chisq=1.0E8;
	for(unsigned int i=0; i<throwns.size(); i++){
		const DMCThrown *thrown = throwns[i];
		DVector3 mc_mom = thrown->momentum();

		if(mom.Mag()<1.0E-9)continue; // Thrown momentum should be > 1eV/c
		
		// Round-off errors can cause cos_phi to be outside the allowed range of -1 to +1
		// Force it to be inside that range in necessary.
		double cos_phi = mom.Dot(mc_mom)/mom.Mag()/mc_mom.Mag();
		if(cos_phi>1.)cos_phi=1.0;
		if(cos_phi<-1.)cos_phi=-1.0;
		
		double delta_pt = (mom.Pt()-mc_mom.Pt())/mc_mom.Pt();
		double delta_theta = (mom.Theta() - mc_mom.Theta())*1000.0; // in milliradians
		double delta_phi = acos(cos_phi)*1000.0; // in milliradians
		double chisq = pow(delta_pt/0.04, 2.0) + pow(delta_theta/20.0, 2.0) + pow(delta_phi/20.0, 2.0);

		if(chisq<min_chisq)myid = thrown->myid;
	}
	
	return myid;
}

//---------------------------------
// GetMCTrackHit
//---------------------------------
const DMCTrackHit* DTrackHitSelectorTHROWN::GetMCTrackHit(const DCoordinateSystem *wire, double rdrift, vector<const DMCTrackHit*> &mctrackhits, int trackno_filter)
{
	/// Identify the DMCTrackHit object (if any) that best matches a hit on 
	/// a given wire with a given drift distance.
	///
	/// This is static function  and so can be used outside of this class.
	/// As such, the list of DMCTrackHits must be supplied by the caller.
	/// Also, the trackno_filter parameter is optional for use on checking
	/// only hits originating from a specific thrown track.

	// We can decide whether this is a CDC or FDC wire based on the z-component.
	// We do this to make it quicker and easier to run through the DMCTrackHits
	// below. Note that wire->udir points along the wire.
	DetectorSystem_t sys = acos(wire->udir.Z())*57.3 > 10.0 ? SYS_FDC:SYS_CDC;

	const DMCTrackHit *best_mchit = NULL;
	double best_resi = 1.0E6;
	double resi_min = 1.0E6, resiu_min=1.0E6;
	for(unsigned int j=0; j<mctrackhits.size(); j++){
		const DMCTrackHit *mchit = mctrackhits[j];
		if(mchit->system!=sys) continue;
		if((trackno_filter>0) && (mchit->track!=trackno_filter)) continue;
		if((trackno_filter>0) && (!mchit->primary)) continue;
		
		// Find the distance from the truth point to the wire. We do this
		// by defining the following 3-vectors:
		//
		// w = hit->wire->origin
		// s = hit->wire->sdir
		// t = hit->wire->tdir
		// u = hit->wire->udir
		//
		// p = mchit truth point (lab coordinates)
		//
		// pw = p - w
		//
		// The matrix "A" is defined from the direction vectors
		// of the wire:
		//
		//       sx  sy  sz
		//  A =  tx  ty  tz
		//       ux  uy  uz
		//
		// The vector abc is defined from the coefficients of each of
		// the wire's direction vectors (s,t,u) such that:
		//
		//  abc*A = pw
		//
		// This lets us solve for abc by inverting A:
		//
		//            -1
		//  abc = pw A
		//
		// Since the "t" direction is along the wire, the "a" and "b"
		// coefficients give the components perpendicular to the wire
		// so adding them in quadrature gives the distance of the
		// point (truth hit) from the wire.
		//
		//  r = a^2 + b^2
		//
		
		DVector3 p(mchit->r*cos(mchit->phi), mchit->r*sin(mchit->phi), mchit->z);
		
		DVector3 w = wire->origin;
		DVector3 sdir = wire->sdir;
		DVector3 tdir = wire->tdir;
		DVector3 udir = wire->udir;
		
		DMatrix pw(3,1);
		pw[0][0] = p.X() - w.X();
		pw[1][0] = p.Y() - w.Y();
		pw[2][0] = p.Z() - w.Z();
		
		DMatrix A(3,3);
		A[0][0]=sdir.X();	A[0][1]=sdir.Y();	A[0][2]=sdir.Z();
		A[1][0]=tdir.X();	A[1][1]=tdir.Y();	A[1][2]=tdir.Z();
		A[2][0]=udir.X();	A[2][1]=udir.Y();	A[2][2]=udir.Z();
		
		DMatrix At(DMatrix::kTransposed, A);
		DMatrix Ainv(DMatrix::kInverted, At);
		DMatrix abc = Ainv*pw;
		
		// Distance of truth hit to wire
		double r = sqrt(pow(abc[0][0], 2.0) + pow(abc[1][0], 2.0));
		
		// Distance of truth hit along the wire
		double u = abc[2][0];
		
		// Now, we need to convert the drift time to a distance and compare it
		// to the distance of the truth hit. Since that includes time-of-flight,
		// the comparison needs to be loose enough to accomodate that.
		// The tof should add at most about 6-7 ns to the time which is less
		// than 400 microns.
		double resi = rdrift - r;

		// Distance of truth point past end of wire in either direction. i.e.
		// a positive value is past the end, a negative value is within the
		// wire length.
		double resiu = fabs(u) - wire->L/2.0;
		
		if(fabs(resi)<fabs(resi_min))resi_min=resi;
		if(resiu<resiu_min)resiu_min=resiu;
		
		if((fabs(resi)<0.3) && (resiu<4.0)){ // Check if truth point is within 2mm of drift time

			// A truth point is roughly the same distance from this wire as indicated
			// by the drift time. Check if this is the best match so far.
			if(fabs(resi)<fabs(best_resi)){
				best_resi = resi;
				best_mchit = mchit;
			}
		}
	}

if(!best_mchit){
//_DBG_<<"resi_min="<<fabs(resi_min)<<"  resiu_min="<<resiu_min<<endl;
}

	return best_mchit;
}

