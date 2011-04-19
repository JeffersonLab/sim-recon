// $Id$
//
//    File: DTrackFitter.cc
// Created: Tue Sep  2 09:10:48 EDT 2008
// Creator: davidl (on Darwin Amelia.local 8.11.1 i386)
//

#include <TRACKING/DReferenceTrajectory.h>
#include <TRACKING/DTrackHitSelector.h>

#include "DTrackFitter.h"
#include "HDGEOMETRY/DRootGeom.h"
using namespace jana;

// hi-res timers for profiling
#include <prof_time.h>
static map<string, prof_time::time_diffs> prof_times;

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
	DEBUG_LEVEL=0;

	CORRECT_FOR_ELOSS=true;
	gPARMS->SetDefaultParameter("TRKFIT:CORRECT_FOR_ELOSS",CORRECT_FOR_ELOSS);

	DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
	if(!dapp){
		_DBG_<<"Cannot get DApplication from JEventLoop! (are you using a JApplication based program?)"<<endl;
		return;
	}
	bfield = dapp->GetBfield(); // this should be run number based!
	lorentz_def=dapp->GetLorentzDeflections();
	geom = dapp->GetDGeometry((loop->GetJEvent()).GetRunNumber());
	RootGeom = dapp->GetRootGeom();

	// Use a special entry to hold how many tracks we fit
	prof_time::time_diffs tdiff_zero;
	prof_times["Ntracks"] = tdiff_zero;
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
	prof_time start_time;
	
	cdchits.clear();
	fdchits.clear();
	fit_type = kWireBased;
	chisq = 1.0E6;
	Ndof=0;
	cdchits_used_in_fit.clear();
	fdchits_used_in_fit.clear();
	fit_status = kFitNotDone;
	
	start_time.TimeDiffNow(prof_times, "Reset");
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
DTrackFitter::fit_status_t DTrackFitter::FitTrack(const DVector3 &pos, const DVector3 &mom, double q, double mass,double t0)
{
	prof_time start_time;

	input_params.setPosition(pos);
	input_params.setMomentum(mom);
	input_params.setCharge(q);
	input_params.setMass(mass);
	input_params.setT0(t0,0.,SYS_NULL);

	DTrackFitter::fit_status_t status = FitTrack();

	start_time.TimeDiffNow(prof_times, "Fit Track 1");

	return status;
}

//-------------------
// FitTrack
//-------------------
DTrackFitter::fit_status_t DTrackFitter::FitTrack(const DKinematicData &starting_params)
{
	prof_time start_time;

	SetInputParameters(starting_params);
	DTrackFitter::fit_status_t status = FitTrack();

	start_time.TimeDiffNow(prof_times, "Fit Track 2");
	
	return status;
}

//-------------------
// FindHitsAndFitTrack
//-------------------
DTrackFitter::fit_status_t DTrackFitter::FindHitsAndFitTrack(const DKinematicData &starting_params, DReferenceTrajectory *rt, JEventLoop *loop, double mass,
							     double t0)
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
	prof_times["Ntracks"].real += 1.0; // keep count of the number of tracks we fit
	prof_time start_time;

	// Reset fitter saving the type of fit we're doing
	fit_type_t save_type = fit_type;
	Reset();
	fit_type = save_type;
	
	// If a mass<0 is passed in, get it from starting_params instead
	if(mass<0.0)mass = starting_params.mass();

	// Correct for energy loss in target etc. based on particle mass in starting_params
	DVector3 pos, mom; // (holds parameters at vertex after correction)
	if(CORRECT_FOR_ELOSS && fit_type==kWireBased){
		jerror_t err = CorrectForELoss(starting_params, rt, pos, mom, mass);
		if(err != NOERROR){
			pos = starting_params.position();
			mom = starting_params.momentum();
		}
	}else{
	  pos = starting_params.position();
	  mom = starting_params.momentum();
	}
	double q=starting_params.charge();

	// Swim a reference trajectory with this candidate's parameters
	rt->Swim(pos, mom, q);
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

	// If the hit selector found no hits at all on the track, the most
	// likely explanation is that the charge of the candidate was wrong,
	// especially for stiff forward-going tracks.  Try rotating phi by 
	// 180 degrees, switching the charge, and trying the hit selector again.
	if (fdchits.size()+cdchits.size()==0){
	  // Swim a reference trajectory with this candidate's parameters
	  mom.SetPhi(mom.Phi()+M_PI);
	  q*=-1.;
	  rt->Swim(pos, mom,q);
	  if(rt->Nswim_steps<1)return fit_status = kFitFailed;
	  hitselector->GetAllHits(input_type, rt, cdctrackhits, fdcpseudos, this);
       	
	  if (fdchits.size()+cdchits.size()!=0){
	    if (DEBUG_LEVEL>0)
	      _DBG_ << "Switching the charge and phi of the track..." <<endl;
	  }
	}
	if (fdchits.size()+cdchits.size()==0) return fit_status=kFitFailed;
	
	

	// In case the subclass doesn't actually set the mass ....
	//fit_params.setMass(starting_params.mass());
	fit_params.setMass(mass);

	start_time.TimeDiffNow(prof_times, "Find Hits");

	// Do the fit
	fit_status = FitTrack(pos, mom,q, mass,t0);

	start_time.TimeDiffNow(prof_times, "Find Hits and Fit Track");
	
	return fit_status;
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
		first_wire = cdchits[cdchits.size()/2]->wire;
	}else{
		vector<const DFDCPseudo*> fdchits;
		starting_params.Get(fdchits);
		if(fdchits.size()!=0){
			sort(fdchits.begin(), fdchits.end(), FDCSortByZincreasing);
			first_wire = fdchits[fdchits.size()/2]->wire;
		}
	}
	if(!first_wire){
		//_DBG_<<"NO WIRES IN CANDIDATE!! (event "<<eventnumber<<")"<<endl;
		return RESOURCE_UNAVAILABLE;
	}

	// Swim from vertex to first wire hit. Disable momentum loss.
	rt->SetDGeometry(geom);
	rt->SetMass(0.0);
	rt->FastSwim(starting_params.position(),starting_params.momentum(),
		     pos,mom,starting_params.charge(), 1000.0, first_wire);

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
	rt->FastSwim(pos,-mom,pos,mom,-starting_params.charge(), 1000.0, 
		     &target);
	rt->SetPLossDirection(DReferenceTrajectory::kForward);
	
	// Reverse momentum
	mom = -mom;

	return NOERROR;
}

// Calculate the density effect correction to the energy loss formulae
double DTrackFitter::CalcDensityEffect(double p,double mass,double density,
				       double Z_over_A,double I){
  double rho_Z_over_A=density*Z_over_A;
  double LnI=log(I);
  return CalcDensityEffect(p,mass,rho_Z_over_A,LnI);
}



// Calculate the density effect correction to the energy loss formulae
double DTrackFitter::CalcDensityEffect(double p,double mass,
				       double rho_Z_over_A,double LnI){

  double betagamma=p/mass;
  return CalcDensityEffect(betagamma,rho_Z_over_A,LnI);
}

// Calculate the density effect correction to the energy loss formulae
double DTrackFitter::CalcDensityEffect(double betagamma,
				       double rho_Z_over_A,double LnI){
  double X=log10(betagamma);
  double X0,X1;
  double Cbar=2.*(LnI-log(28.816e-9*sqrt(rho_Z_over_A)))+1.;
  if (rho_Z_over_A>0.01){ // not a gas
    if (LnI<-1.6118){ // I<100
      if (Cbar<=3.681) X0=0.2;
      else X0=0.326*Cbar-1.;
      X1=2.;
    }
    else{
      if (Cbar<=5.215) X0=0.2;
      else X0=0.326*Cbar-1.5;
      X1=3.;
    }
  }
  else { // gases
    X1=4.;
    if (Cbar<=9.5) X0=1.6;
    else if (Cbar>9.5 && Cbar<=10.) X0=1.7;
    else if (Cbar>10 && Cbar<=10.5) X0=1.8;    
    else if (Cbar>10.5 && Cbar<=11.) X0=1.9;
    else if (Cbar>11.0 && Cbar<=12.25) X0=2.;
    else if (Cbar>12.25 && Cbar<=13.804){
      X0=2.;
      X1=5.;
    }
    else {
      X0=0.326*Cbar-2.5;
      X1=5.;
    } 
  }
  double delta=0;
  if (X>=X0 && X<X1)
    delta=4.606*X-Cbar+(Cbar-4.606*X0)*pow((X1-X)/(X1-X0),3.);
  else if (X>=X1)
    delta= 4.606*X-Cbar;  
  return delta;
}

//------------------
// GetProfilingTimes
//------------------
void DTrackFitter::GetProfilingTimes(map<string, prof_time::time_diffs> &my_prof_times) const
{
	my_prof_times = prof_times;
}
