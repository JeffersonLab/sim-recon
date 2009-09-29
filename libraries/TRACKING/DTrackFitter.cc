// $Id$
//
//    File: DTrackFitter.cc
// Created: Tue Sep  2 09:10:48 EDT 2008
// Creator: davidl (on Darwin Amelia.local 8.11.1 i386)
//

#include <TRACKING/DReferenceTrajectory.h>
#include <TRACKING/DTrackHitSelector.h>

#include "DTrackFitter.h"
using namespace jana;

extern bool FDCSortByZincreasing(const DFDCPseudo* const &hit1, const DFDCPseudo* const &hit2);
extern bool CDCSortByRincreasing(const DCDCTrackHit* const &hit1, const DCDCTrackHit* const &hit2);

//-------------------
// DTrackFitter  (Constructor)
//-------------------
DTrackFitter::DTrackFitter(JEventLoop *loop)
{
	this->loop = loop;
	bfield=NULL;
	fit_status = kFitNotDone;
	
	DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
	if(!dapp){
		_DBG_<<"Cannot get DApplication from JEventLoop! (are you using a JApplication based program?)"<<endl;
		return;
	}
	bfield = dapp->GetBfield(); // this should be run number based!
	lorentz_def=dapp->GetLorentzDeflections();
	geom = dapp->GetDGeometry((loop->GetJEvent()).GetRunNumber());
	RootGeom = dapp->GetRootGeom();
}

//-------------------
// ~DTrackFitter (Destructor)
//-------------------
DTrackFitter::~DTrackFitter()
{
}

//-------------------
// Reset
//-------------------
void DTrackFitter::Reset(void)
{
	cdchits.clear();
	fdchits.clear();
	fit_type = kWireBased;
	chisq = 1.0E6;
	Ndof=0;
	cdchits_used_in_fit.clear();
	fdchits_used_in_fit.clear();
	fit_status = kFitNotDone;
}

//-------------------
// AddHit
//-------------------
void DTrackFitter::AddHit(const DCDCTrackHit* cdchit)
{
	cdchits.push_back(cdchit);
	fit_status = kFitNotDone;
}

//-------------------
// AddHits
//-------------------
void DTrackFitter::AddHits(vector<const DCDCTrackHit*> cdchits)
{
	for(unsigned int i=0; i<cdchits.size(); i++)this->cdchits.push_back(cdchits[i]);
	fit_status = kFitNotDone;
}

//-------------------
// AddHit
//-------------------
void DTrackFitter::AddHit(const DFDCPseudo* fdchit)
{
	fdchits.push_back(fdchit);
	fit_status = kFitNotDone;
}

//-------------------
// AddHits
//-------------------
void DTrackFitter::AddHits(vector<const DFDCPseudo*> fdchits)
{
	for(unsigned int i=0; i<fdchits.size(); i++)this->fdchits.push_back(fdchits[i]);
	fit_status = kFitNotDone;
}

//-------------------
// FitTrack
//-------------------
DTrackFitter::fit_status_t DTrackFitter::FitTrack(const DVector3 &pos, const DVector3 &mom, double q, double mass)
{
	input_params.setPosition(pos);
	input_params.setMomentum(mom);
	input_params.setCharge(q);
	input_params.setMass(mass);
	
	return FitTrack();
}

//-------------------
// FitTrack
//-------------------
DTrackFitter::fit_status_t DTrackFitter::FitTrack(const DKinematicData &starting_params)
{
	SetInputParameters(starting_params);
	
	return FitTrack();
}

//-------------------
// FindHitsAndFitTrack
//-------------------
DTrackFitter::fit_status_t DTrackFitter::FindHitsAndFitTrack(const DKinematicData &starting_params, DReferenceTrajectory *rt, JEventLoop *loop, double mass)
{
	/// Fit a DTrackCandidate using a given mass hypothesis.
	///
	/// This will perform a full wire-based and time-based
	/// fit using the given mass and starting from the given
	/// candidate. The given DReferenceTrajectory is used to
	/// swim the track numerous times during the various stages
	/// but will be left with the final time-based fit result.
	/// The JEventLoop given will be used to get the hits (CDC
	/// and FDC) and default DTrackHitSelector to use for the
	/// fit.

	// Reset fitter saving the type of fit we're doing
	fit_type_t save_type = fit_type;
	Reset();
	fit_type = save_type;
	
	// If a mass<0 is passed in, get it from starting_params instead
	if(mass<0.0)mass = starting_params.mass();

	// Correct for energy loss in target etc. based on particle mass in starting_params
	DVector3 pos, mom; // (holds parameters at vertex after correction)
	if(fit_type==kTimeBased){
		CorrectForELoss(starting_params, rt, pos, mom, mass);
	}else{
		pos = starting_params.position();
		mom = starting_params.momentum();
	}

	// Swim a reference trajectory with this candidate's parameters
	rt->Swim(pos, mom, starting_params.charge());
	if(rt->Nswim_steps<1)return fit_status = kFitFailed;

	// Get pointer to DTrackHitSelector object
	vector<const DTrackHitSelector *> hitselectors;
	loop->Get(hitselectors);
	if(hitselectors.size()<1){
		_DBG_<<"Unable to get a DTrackHitSelector object! NO Charged track fitting will be done!"<<endl;
		return fit_status = kFitNotDone;
	}
	const DTrackHitSelector * hitselector = hitselectors[0];

	// Get hits to be used for the fit
	vector<const DCDCTrackHit*> cdctrackhits;
	vector<const DFDCPseudo*> fdcpseudos;
	loop->Get(cdctrackhits);
	loop->Get(fdcpseudos);
	DTrackHitSelector::fit_type_t input_type = fit_type==kTimeBased ? DTrackHitSelector::kWireBased:DTrackHitSelector::kHelical;
	hitselector->GetAllHits(input_type, rt, cdctrackhits, fdcpseudos, this);

	// In case the subclass doesn't actually set the mass ....
	fit_params.setMass(starting_params.mass());

	// Do the fit
	return fit_status = FitTrack(pos, mom, starting_params.charge(), mass);	
}

//------------------
// CorrectForELoss
//------------------
jerror_t DTrackFitter::CorrectForELoss(const DKinematicData &starting_params, DReferenceTrajectory *rt, DVector3 &pos, DVector3 &mom, double mass)
{
	// Find first wire hit by this track
	const DCoordinateSystem *first_wire = NULL;
	vector<const DCDCTrackHit*> cdchits;
	starting_params.Get(cdchits);
	if(cdchits.size()>0){
		sort(cdchits.begin(), cdchits.end(), CDCSortByRincreasing);
		first_wire = cdchits[0]->wire;
	}else{
		vector<const DFDCPseudo*> fdchits;
		starting_params.Get(fdchits);
		if(fdchits.size()!=0){
			sort(fdchits.begin(), fdchits.end(), FDCSortByZincreasing);
			first_wire = fdchits[0]->wire;
		}
	}
	if(!first_wire){
		//_DBG_<<"NO WIRES IN CANDIDATE!! (event "<<eventnumber<<")"<<endl;
		return RESOURCE_UNAVAILABLE;
	}

	// Swim from vertex to first wire hit. Disable momentum loss.
	rt->SetDGeometry(geom);
	rt->SetMass(0.0);
	rt->Swim(starting_params.position(), starting_params.momentum(), starting_params.charge(), 1000.0, first_wire);
	rt->DistToRT(first_wire);
	rt->GetLastDOCAPoint(pos, mom);

	// Define target center
	DCoordinateSystem target;
	target.origin.SetXYZ(0.0, 0.0, 65.0);
	target.sdir.SetXYZ(1.0, 0.0, 0.0);
	target.tdir.SetXYZ(0.0, 1.0, 0.0);
	target.udir.SetXYZ(0.0, 0.0, 1.0);
	target.L = 30.0;

	// Swim backwards to target, setting momentum to *increase* due to material
	rt->SetMass(mass);
	rt->SetPLossDirection(DReferenceTrajectory::kBackward);
	rt->Swim(pos, -mom, -starting_params.charge(), 1000.0, &target);
	rt->SetPLossDirection(DReferenceTrajectory::kForward);
	rt->DistToRT(&target);
	rt->GetLastDOCAPoint(pos, mom);
	
	// Reverse momentum
	mom = -mom;

	return NOERROR;
}

// Calculate the path length for a single hit in a straw and return ds and the 
// energy deposition in the straw.  It returns dE as the first element in the 
// dedx pair and ds as the second element in the dedx pair.
jerror_t DTrackFitter::CalcdEdxHit(const DVector3 &mom,
				   const DVector3 &pos,
				   const DCDCTrackHit *hit,
				   pair <double,double> &dedx){
  if (hit==NULL || hit->wire==NULL) return RESOURCE_UNAVAILABLE;

  // Track direction parameters
  double phi=mom.Phi();
  double lambda=M_PI_2-mom.Theta();
  double cosphi=cos(phi);
  double sinphi=sin(phi);
  double tanl=tan(lambda);
  
  //Position relative to wire origin
  double dz=pos.z()-hit->wire->origin.z();
  double dx=pos.x()-hit->wire->origin.x();
  double dy=pos.y()-hit->wire->origin.y();
  
  // square of straw radius
  double rs2=0.8*0.8;
  
  // Useful temporary variables related to the direction of the wire
  double A=1.-hit->wire->udir.x()*hit->wire->udir.x();
  double B=-2.*hit->wire->udir.x()*hit->wire->udir.y();
  double C=-2.*hit->wire->udir.x()*hit->wire->udir.z();
  double D=-2.*hit->wire->udir.y()*hit->wire->udir.z();
  double E=1.-hit->wire->udir.y()*hit->wire->udir.y();
  double F=1.-hit->wire->udir.z()*hit->wire->udir.z();

  // The path length in the straw is given by  s=sqrt(b*b-4*a*c)/a/cosl.
  // a, b, and c follow.
  double a=A*cosphi*cosphi+B*cosphi*sinphi+C*cosphi*tanl+D*sinphi*tanl
    +E*sinphi*sinphi+F*tanl*tanl;
  double b=2.*A*dx*cosphi+B*dx*sinphi+B*dy*cosphi+C*dx*tanl+C*cosphi*dz
    +D*dy*tanl+D*sinphi*dz+2.*E*dy*sinphi+2.*F*dz*tanl;
  double c=A*dx*dx+B*dx*dy+C*dx*dz+D*dy*dz+E*dy*dy+F*dz*dz-rs2;
  
  // Check for valid arc length and compute dEdx
  double temp=b*b-4.*a*c;
  if (temp>0){
    double cosl=fabs(cos(lambda));
    double gas_density=0.0018;

    // arc length and energy deposition
    dedx.second=gas_density*sqrt(temp)/a/cosl; // g/cm^2
    dedx.first=1000.*hit->dE; //MeV

    return NOERROR;
  }
  
  return VALUE_OUT_OF_RANGE;
}
