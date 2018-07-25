// $Id$
//
//    File: DTrackFitter.cc
// Created: Tue Sep  2 09:10:48 EDT 2008
// Creator: davidl (on Darwin Amelia.local 8.11.1 i386)
//

#include <TRACKING/DReferenceTrajectory.h>
#include <TRACKING/DTrackHitSelector.h>


#include "DTrackFitter.h"
#include "PID/DParticleID.h"
#include "HDGEOMETRY/DRootGeom.h"
using namespace jana;

// hi-res timers for profiling
// The PROFILE_TRK_TIMES compile time option is intended to enable keeping track
// of the time needed to do the tracking using three system-level interval timers.  
#ifdef PROFILE_TRK_TIMES
#include <prof_time.h>
static map<string, prof_time::time_diffs> prof_times;
#endif

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
	unsigned int run_number = (loop->GetJEvent()).GetRunNumber();
	DEBUG_LEVEL=0;

	CORRECT_FOR_ELOSS=true;
	gPARMS->SetDefaultParameter("TRKFIT:CORRECT_FOR_ELOSS",CORRECT_FOR_ELOSS);


	DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
	if(!dapp){
		_DBG_<<"Cannot get DApplication from JEventLoop! (are you using a JApplication based program?)"<<endl;
		return;
	}
	bfield = dapp->GetBfield(run_number); 
	geom = dapp->GetDGeometry(run_number);

	RootGeom=NULL;
	MATERIAL_MAP_MODEL = "DGeometry";
	gPARMS->SetDefaultParameter("TRKFIT:MATERIAL_MAP_MODEL",MATERIAL_MAP_MODEL);
	if(MATERIAL_MAP_MODEL=="DRootGeom"){
	  RootGeom = dapp->GetRootGeom(run_number);
	}
	// Create the extrapolation vectors
	vector<Extrapolation_t>myvector;
	extrapolations.emplace(SYS_BCAL,myvector);
	extrapolations.emplace(SYS_TOF,myvector);
	extrapolations.emplace(SYS_FCAL,myvector);
	extrapolations.emplace(SYS_FDC,myvector);
	extrapolations.emplace(SYS_CDC,myvector);
	extrapolations.emplace(SYS_START,myvector);
	extrapolations.emplace(SYS_DIRC,myvector);	

	extrapolations[SYS_TOF].reserve(1);
	extrapolations[SYS_BCAL].reserve(300);
	extrapolations[SYS_FCAL].reserve(2);
	extrapolations[SYS_FDC].reserve(24);
	extrapolations[SYS_CDC].reserve(200);
	extrapolations[SYS_START].reserve(1);
	extrapolations[SYS_DIRC].reserve(1);
	
	pulls.reserve(30);

#ifdef PROFILE_TRK_TIMES
	// Use a special entry to hold how many tracks we fit
	prof_time::time_diffs tdiff_zero;
	prof_times["Ntracks"] = tdiff_zero;
#endif	
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
#ifdef PROFILE_TRK_TIMES
    prof_time start_time;
#endif

	cdchits.clear();
	fdchits.clear();
	fit_type = kWireBased;
	chisq = 1.0E6;
	Ndof=0;
	cdchits_used_in_fit.clear();
	fdchits_used_in_fit.clear();
	ClearExtrapolations();
	pulls.clear();
      
	fit_status = kFitNotDone;
	
#ifdef PROFILE_TRK_TIMES
	start_time.TimeDiffNow(prof_times, "Reset");
#endif	
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
DTrackFitter::fit_status_t DTrackFitter::FitTrack(const DVector3 &pos, const DVector3 &mom, double q, double mass,double t0,DetectorSystem_t t0_det)
{
#ifdef PROFILE_TRK_TIMES
    prof_time start_time;
#endif	

	input_params.setPosition(pos);
	input_params.setMomentum(mom);
	input_params.setPID(IDTrack(q, mass));
	input_params.setTime(t0);
	input_params.setT0(t0,0.,t0_det);

	DTrackFitter::fit_status_t status = FitTrack();

#ifdef PROFILE_TRK_TIMES
	start_time.TimeDiffNow(prof_times, "Fit Track 1");
#endif

	return status;
}

//-------------------
// FitTrack
//-------------------
DTrackFitter::fit_status_t DTrackFitter::FitTrack(const DTrackingData &starting_params)
{
#ifdef PROFILE_TRK_TIMES
  prof_time start_time;
#endif

	SetInputParameters(starting_params);
	DTrackFitter::fit_status_t status = FitTrack();

#ifdef PROFILE_TRK_TIMES
	start_time.TimeDiffNow(prof_times, "Fit Track 2");
#endif

	return status;
}
//-------------------
// FindHitsAndFitTrack
//-------------------
DTrackFitter::fit_status_t 
DTrackFitter::FindHitsAndFitTrack(const DKinematicData &starting_params, 
				  const map<DetectorSystem_t,vector<DTrackFitter::Extrapolation_t> >&extrapolations,
				  JEventLoop *loop, 
				  double mass,int N,double t0,
				  DetectorSystem_t t0_det){
  // Reset fitter saving the type of fit we're doing
  fit_type_t save_type = fit_type;
  Reset();
  fit_type = save_type;
	
  // If a mass<0 is passed in, get it from starting_params instead
  if(mass<0.0)mass = starting_params.mass();
  // charge of the track
  double q=starting_params.charge();

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

  // Get Bfield at the position at the middle of the extrapolations, i.e. the 
  // region where we actually have measurements...
  bool got_hits=false;
  if (extrapolations.at(SYS_CDC).size()>0){
    vector<Extrapolation_t>extraps=extrapolations.at(SYS_CDC);
    DVector3 mypos=extraps[extraps.size()/2].position;
    double Bz=GetDMagneticFieldMap()->GetBz(mypos.x(),mypos.y(),mypos.z());
    hitselector->GetCDCHits(Bz,q,extraps,cdctrackhits,this,N);
    got_hits=true;
  }
  if (extrapolations.at(SYS_FDC).size()>0){
    vector<Extrapolation_t>extraps=extrapolations.at(SYS_FDC);
    DVector3 mypos=extraps[extraps.size()/2].position;
    double Bz=GetDMagneticFieldMap()->GetBz(mypos.x(),mypos.y(),mypos.z());
    hitselector->GetFDCHits(Bz,q,extraps,fdcpseudos,this,N);	
    got_hits=true;
  }
  if (got_hits==false){
    return fit_status = kFitNotDone;
  }

  // In case the subclass doesn't actually set the mass ....
  fit_params.setPID(IDTrack(q, mass));
  
#ifdef PROFILE_TRK_TIMES
  start_time.TimeDiffNow(prof_times, "Find Hits");
#endif
  
  // Do the fit 
  DVector3 pos = starting_params.position();
  DVector3 mom = starting_params.momentum();
  fit_status = FitTrack(pos, mom,q, mass,t0,t0_det);
  
#ifdef PROFILE_TRK_TIMES
  start_time.TimeDiffNow(prof_times, "Find Hits and Fit Track");
#endif
  
  return fit_status;
}
//-------------------
// FindHitsAndFitTrack
//-------------------
DTrackFitter::fit_status_t 
DTrackFitter::FindHitsAndFitTrack(const DKinematicData &starting_params,
				  const DReferenceTrajectory *rt, JEventLoop *loop, 
				  double mass,int N,double t0,
				  DetectorSystem_t t0_det)
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
#ifdef PROFILE_TRK_TIMES
  prof_times["Ntracks"].real += 1.0; // keep count of the number of tracks we fit

  prof_time start_time;
#endif

	// Reset fitter saving the type of fit we're doing
	fit_type_t save_type = fit_type;
	Reset();
	fit_type = save_type;
	
	// If a mass<0 is passed in, get it from starting_params instead
	if(mass<0.0)mass = starting_params.mass();

	// Correct for energy loss in target etc. based on particle mass in starting_params
	DVector3 pos, mom; // (holds parameters at vertex after correction)
	/*
	if(CORRECT_FOR_ELOSS && fit_type==kWireBased){
		jerror_t err = CorrectForELoss(starting_params, rt, pos, mom, mass);
		if(err != NOERROR){
			pos = starting_params.position();
			mom = starting_params.momentum();
		}
	}else{
	*/
	  pos = starting_params.position();
	  mom = starting_params.momentum();
	  //}
	double q=starting_params.charge();
       
	// Swim a reference trajectory with this candidate's parameters
	//rt->Swim(pos, mom, q);
	//if(rt->Nswim_steps<1)return fit_status = kFitFailed;

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
	hitselector->GetAllHits(input_type, rt, cdctrackhits, fdcpseudos, this,N);

	// If the condition below is met, it seems that the track parameters 
	// are inconsistent with the hits used to create the track candidate, 
	// in which case we have to bail...
	if (fdchits.size()+cdchits.size()==0) return fit_status=kFitNotDone;
	
	// In case the subclass doesn't actually set the mass ....
	fit_params.setPID(IDTrack(q, mass));

#ifdef PROFILE_TRK_TIMES
	start_time.TimeDiffNow(prof_times, "Find Hits");
#endif

	// Do the fit
	fit_status = FitTrack(pos, mom,q, mass,t0,t0_det);

#ifdef PROFILE_TRK_TIMES
	start_time.TimeDiffNow(prof_times, "Find Hits and Fit Track");
#endif

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


// Extrapolate to a radius R given two extrapolation points before and after 
bool DTrackFitter::ExtrapolateToRadius(double R,
				       const vector<Extrapolation_t>&extraps,
				       DVector3 &pos,DVector3 &mom,double &t,
				       double &s) const{
  if (extraps.size()<2) return false;

  for (unsigned int j=1;j<extraps.size();j++){
    if (extraps[j].position.Perp()>R){
      // At this point, the location where the track intersects the cyclinder 
      // is somewhere between extrapolated point and the previous one.  For
      // simplicity, we're going to just find the intersection of the cylinder 
      // with the line that joins the 2 positions. We do this by working in 
      // the X/Y plane only and finding the value of "alpha" which is the 
      // fractional distance the intersection point is between our two 
      // extrapolations.  We'll then apply the alpha found in the 2D X/Y space 
      // to the 3D x/y/Z space to find the actual intersection point.
      Extrapolation_t extrap=extraps[j];
      Extrapolation_t prev=extraps[j-1];
      pos=extrap.position;
      mom=extrap.momentum;
      t=extrap.t;
      s=extrap.s;
      // The next part of the code refines the extrapolation
      DVector3 prevpos=prev.position;
      DVector2 x1(pos.X(),pos.Y());
      DVector2 x2(prevpos.X(),prevpos.Y());
      DVector2 dx = x2-x1;
      double A = dx.Mod2();
      double B = 2.0*(x1.X()*dx.X() + x1.Y()*dx.Y());
      double C = x1.Mod2() - R*R;
      
      double sqrt_D=sqrt(B*B-4.0*A*C);
      double one_over_denom=0.5/A;
      double alpha1 = (-B + sqrt_D)*one_over_denom;
      double alpha2 = (-B - sqrt_D)*one_over_denom;
      double alpha = alpha1;
      if(alpha1<0.0 || alpha1>1.0)alpha=alpha2;
      if (isfinite(alpha)){
	DVector3 delta = prevpos - pos;
	pos+=alpha*delta;
	// Flight time and path length (approximate)
	double ds=(1.-alpha)*delta.Mag();
	double v=(extrap.s-prev.s)/(extrap.t-prev.t);
	s-=ds;
	t-=ds/v;

	return true;
      }
      break;
    }
  }
  return false;
}
			
bool DTrackFitter::ExtrapolateToRadius(double R,
				       const vector<Extrapolation_t>&extraps,
				       DVector3 &pos) const{
  double s=0,t=0;
  DVector3 mom;
  return ExtrapolateToRadius(R,extraps,pos,mom,s,t);
}

// Loop through extrapolations to find the distance of closest approach to a 
// wire.
double DTrackFitter::DistToWire(const DCoordinateSystem *wire,
				const vector<Extrapolation_t>&extrapolations,
				DVector3 *pos, DVector3 *mom,
				DVector3 *position_along_wire) const{
  if (extrapolations.size()<3) return 1000.;

  // Wire info
  double z0w=wire->origin.z();
  double ux=wire->udir.x();
  double uy=wire->udir.y();
  double uz=wire->udir.z();

  double doca_old=1000.,doca=1000.;
  for (unsigned int i=1;i<extrapolations.size();i++){
    DVector3 trackpos=extrapolations[i].position;
    double z=trackpos.z();
    double dzw=z-z0w;
    DVector3 wirepos=wire->origin+(dzw/uz)*wire->udir;    
    doca=(wirepos-trackpos).Perp();
    if (doca>doca_old){
      DVector3 trackdir=extrapolations[i-1].momentum;
      trackdir.SetMag(1.);
      // Variables relating wire direction and track direction
      double lambda=M_PI_2-trackdir.Theta();
      double sinl=sin(lambda);
      double cosl=cos(lambda);
      double phi=trackdir.Phi();
      double sinphi=sin(phi);
      double cosphi=cos(phi);
      double my_ux=ux*sinl-cosl*cosphi;
      double my_uy=uy*sinl-cosl*sinphi;
      double denom=my_ux*my_ux+my_uy*my_uy;
      double ds=((trackpos.X()-wire->origin.X()-ux*dzw)*my_ux
		 +(trackpos.Y()-wire->origin.Y()-uy*dzw)*my_uy)/denom;
      if (fabs(ds)<2.*fabs(extrapolations[i].s-extrapolations[i-1].s)){
	trackpos+=ds*trackdir;
	wirepos=wire->origin+((trackpos.z()-z0w)/uz)*wire->udir;
	double cosstereo=wire->udir.Dot(DVector3(0.,0.,1.));
	doca=(wirepos-trackpos).Perp()*cosstereo;

	if (pos!=NULL) *pos=trackpos;
	if (mom!=NULL) *mom=extrapolations[i-1].momentum;
	if (position_along_wire!=NULL) *position_along_wire=wirepos;
      }
      break;
    }
    
    doca_old=doca;
  }

  return doca;
}



//------------------
// GetProfilingTimes
//------------------
#ifdef PROFILE_TRK_TIMES
void DTrackFitter::GetProfilingTimes(map<string, prof_time::time_diffs> &my_prof_times) const
{
	my_prof_times = prof_times;
}
#endif
