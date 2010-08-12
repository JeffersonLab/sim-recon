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

#define TCUT 100.e-6 // energy cut for bethe-bloch 

extern bool FDCSortByZincreasing(const DFDCPseudo* const &hit1, const DFDCPseudo* const &hit2);
extern bool CDCSortByRincreasing(const DCDCTrackHit* const &hit1, const DCDCTrackHit* const &hit2);

// Routine for sorting dEdx data
bool static DTrackFitter_dedx_cmp(pair<double,double>a,pair<double,double>b){
  double dEdx1=a.first/a.second;
  double dEdx2=b.first/b.second;
  return dEdx1<dEdx2;  
}

bool static DTrackFitter_dedx_cmp2(DTrackFitter::dedx_t a,DTrackFitter::dedx_t b){
  double dEdx1=a.dE/a.dx;
  double dEdx2=b.dE/b.dx;
  return dEdx1<dEdx2;
}

//-------------------
// DTrackFitter  (Constructor)
//-------------------
DTrackFitter::DTrackFitter(JEventLoop *loop)
{
	this->loop = loop;
	bfield=NULL;
	fit_status = kFitNotDone;
	DEBUG_LEVEL=0;

	DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
	if(!dapp){
		_DBG_<<"Cannot get DApplication from JEventLoop! (are you using a JApplication based program?)"<<endl;
		return;
	}
	bfield = dapp->GetBfield(); // this should be run number based!
	lorentz_def=dapp->GetLorentzDeflections();
	geom = dapp->GetDGeometry((loop->GetJEvent()).GetRunNumber());
	RootGeom = dapp->GetRootGeom();

	// Get material properties for chamber gas
	double rho_Z_over_A_LnI=0,radlen=0;
	RootGeom->FindMat("CDchamberGas",mRhoZoverAGas,rho_Z_over_A_LnI,
			  radlen);
	mLnIGas=rho_Z_over_A_LnI/mRhoZoverAGas;
	mKRhoZoverAGas=0.1535*mRhoZoverAGas;
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
DTrackFitter::fit_status_t DTrackFitter::FitTrack(const DVector3 &pos, const DVector3 &mom, double q, double mass,double t0)
{
	input_params.setPosition(pos);
	input_params.setMomentum(mom);
	input_params.setCharge(q);
	input_params.setMass(mass);
	input_params.setT0(t0,0.,SYS_NULL);
	
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

	// Reset fitter saving the type of fit we're doing
	fit_type_t save_type = fit_type;
	Reset();
	fit_type = save_type;
	
	// If a mass<0 is passed in, get it from starting_params instead
	if(mass<0.0)mass = starting_params.mass();

	// Correct for energy loss in target etc. based on particle mass in starting_params
	DVector3 pos, mom; // (holds parameters at vertex after correction)
	if(fit_type==kWireBased){
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

	// Do the fit
	return fit_status = FitTrack(pos, mom,q, mass,t0);
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
	rt->Swim(starting_params.position(), starting_params.momentum(), starting_params.charge(), 1000.0, first_wire);
	if(rt->Nswim_steps<2)return NOERROR; // no enough swim steps to make this correction worthwile
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
	if(rt->Nswim_steps<2)return NOERROR; // no enough swim steps to make this correction worthwile
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
  double rs2=0.776*0.776;
  
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
    //    double gas_density=0.0018;

    // arc length and energy deposition
    //dedx.second=gas_density*sqrt(temp)/a/cosl; // g/cm^2
    dedx.second=sqrt(temp)/a/cosl;
    dedx.first=1000.*hit->dE; //MeV

    return NOERROR;
  }
  
  return VALUE_OUT_OF_RANGE;
}

// Calculate the variance in the energy loss in a Gaussian approximation
// in units of MeV^2/cm^2
// The full width at half maximum of the energy loss distribution is
// approximated by Gamma=4.018 Xi, where
//      Xi=0.1535*density*(Z/A)*x/beta^2  [MeV]
// To convert that to the sigma of a Gaussian, use Gamma=2.354*sigma.
double DTrackFitter::GetdEdxVariance(double p,double mass_hyp,double dx,
				     DVector3 pos){
  if (p<0.001) return 1000.;  // try to avoid division by zero errors
   
  double beta2=1./(1.+mass_hyp*mass_hyp/p/p);
 
  // Get material properties
  double rho_Z_over_A=0.,rho_Z_over_A_LnI=0.,radlen=0.;
  geom->FindMat(pos,rho_Z_over_A,rho_Z_over_A_LnI,radlen);

  // Ignore errors in dx for now 
  const double full_width_half_maximum=4.018;
  double sigma=0.1535*full_width_half_maximum/2.354*rho_Z_over_A/beta2;
  return sigma*sigma;
}

// Calculate the most probable energy loss per unit length in units of 
// MeV/cm in a material with a certain density, mean excitation energy 
// I (in eV), and Z/A for a particle of momentum p and mass mass_hyp
double DTrackFitter::GetdEdx(double p,double mass_hyp,double dx,
			     DVector3 pos){
  if (p<0.001) return 0.;  // try to avoid division by zero errors
  double betagamma=p/mass_hyp;
  double beta2=1./(1.+1./betagamma/betagamma);

  // Electron mass 
  double Me=0.000511; //GeV
  
   // Get material properties
  double rho_Z_over_A=0.,rho_Z_over_A_LnI=0.,radlen=0.;
  geom->FindMat(pos,rho_Z_over_A,rho_Z_over_A_LnI,radlen);

  // First (non-logarithmic) term in Bethe-Bloch formula
  double mean_dedx=0.1535*rho_Z_over_A/beta2;

  // density effect
  double LnI=rho_Z_over_A_LnI/rho_Z_over_A;
  double delta=CalcDensityEffect(betagamma,rho_Z_over_A,LnI);
  
  // Most probable energy loss from Landau theory (see Leo, pp. 51-52)
  return mean_dedx*(log(mean_dedx*dx/1000.)-log((1.-beta2)/2./Me/beta2)
		    -2.*LnI-beta2+0.198-delta);
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

// Calculate the most probable energy loss per unit length in units of 
// MeV/cm in the FDC or CDC gas for a particle of momentum p and mass mass_hyp
double DTrackFitter::GetdEdx(double p,double mass_hyp,double mean_path_length){
  double betagamma=p/mass_hyp;
  double beta2=1./(1.+1./betagamma/betagamma);
  if (beta2<1e-6) beta2=1e-6;
  
  // Electron mass 
  double Me=0.000511; //GeV

  // First (non-logarithmic) term in Bethe-Bloch formula
  double mean_dedx=mKRhoZoverAGas/beta2;
 
  // density effect
  double delta=CalcDensityEffect(betagamma,mRhoZoverAGas,mLnIGas);

  // Most probable energy loss from Landau theory (see Leo, pp. 51-52)
  return mean_dedx*(log(mean_dedx*mean_path_length/1000.)
			     -log((1.-beta2)/2./Me/beta2)
			     -2.*mLnIGas-beta2+0.198-delta);
}

// Empirical form for sigma/mean for gaseous detectors with num_dedx 
// samples and sampling thickness path_length.  Assumes that the number of 
// hits has already been converted from an (unsigned) int to a double.
double DTrackFitter::GetdEdxSigma(double num_hits,double p,double mass,
				  double mean_path_length){
  // kinematic quantities
  double betagamma=p/mass;
  double betagamma2=betagamma*betagamma;
  double beta2=1./(1.+1./betagamma2);
  if (beta2<1e-6) beta2=1e-6;

  double Me=0.000511; //GeV
  double m_ratio=Me/mass;
  double two_Me_betagamma_sq=2.*Me*betagamma2;
  double Tmax
    =two_Me_betagamma_sq/(1.+2.*sqrt(1.+betagamma2)*m_ratio+m_ratio*m_ratio);
  // Energy truncation for knock-on electrons
  double T0=(Tmax>TCUT)?TCUT:Tmax;

  //density effect   
  double delta=CalcDensityEffect(p,mass,mRhoZoverAGas,mLnIGas); 
  
  // Bethe-Bloch
  double mean_dedx=mKRhoZoverAGas/beta2
    *(log(two_Me_betagamma_sq*T0)-2.*mLnIGas-beta2*(1.+T0/Tmax)-delta);
  
  return 0.41*mean_dedx*pow(num_hits,-0.43)*pow(mean_path_length,-0.32);
  //return 0.41*mean_dedx*pow(double(num_hits),-0.5)*pow(mean_path_length,-0.32);
}
// Overload to allow the first argument to be an integer
double DTrackFitter::GetdEdxSigma(unsigned int num_hits,double p,
				  double mass,double mean_path_length){
  double N=double(num_hits);
  return GetdEdxSigma(N,p,mass,mean_path_length);
}



// Get the most probable value and the standard deviation of the energy loss
jerror_t DTrackFitter::GetdEdxMPandSigma(unsigned int num_hits,double p,
					 double mass,double mean_path_length,
					 double &dedx_mp,double &dedx_sigma){
   // kinematic quantities
  double betagamma=p/mass;
  double betagamma2=betagamma*betagamma;
  double beta2=1./(1.+1./betagamma2);
  if (beta2<1e-6) beta2=1e-6;

  double Me=0.000511; //GeV
  double m_ratio=Me/mass;
  double two_Me_betagamma_sq=2.*Me*betagamma2;
  double Tmax
    =two_Me_betagamma_sq/(1.+2.*sqrt(1.+betagamma2)*m_ratio+m_ratio*m_ratio);
  // Energy truncation for knock-on electrons
  double T0=(Tmax>TCUT)?TCUT:Tmax;

  //density effect   
  double delta=CalcDensityEffect(betagamma,mRhoZoverAGas,mLnIGas); 
  
  // First (non-logarithmic) term in Bethe-Bloch formula
  double Xi=mKRhoZoverAGas/beta2;

  // Most probable energy loss from Landau theory (see Leo, pp. 51-52)
  dedx_mp=Xi*(log(Xi*mean_path_length/1000.)-log((1.-beta2)/2./Me/beta2)
	      -2.*mLnIGas-beta2+0.198-delta);

  // Bethe-Bloch: mean energy loss
  double mean_dedx
    =Xi*(log(two_Me_betagamma_sq*T0)-2.*mLnIGas-beta2*(1.+T0/Tmax)-delta);

  dedx_sigma=0.41*mean_dedx*pow(double(num_hits),-0.43)*pow(mean_path_length,-0.32);
  
  return NOERROR;
}


// Compute the dEdx for a track described by the reference trajectory rt using
// a subset of all the FDC and CDC hits on the track. Also returns the mean
// path length per hit in the active volume of the detector and the average
// measured momentum within the active region.
jerror_t DTrackFitter::GetdEdx(const DReferenceTrajectory *rt, double &dedx,
			       double &mean_path_length, double &p_avg,
			       unsigned int &num_hits){
  DVector3 pos,mom;
  //Get the list of cdc hits used in the fit
  vector<const DCDCTrackHit*>cdchits;
  cdchits=GetCDCFitHits();

  //Vector of dE and dx pairs
  vector<pair<double,double> >dEdx_list;
  pair<double,double>de_and_dx;

  // Average measured momentum
  p_avg=0.;

  // Initialize other output variables
  dedx=mean_path_length=0.;
  num_hits=0;

  // We cast away the const-ness of the reference trajectory so that we can use the DisToRT method
  DReferenceTrajectory *my_rt=const_cast<DReferenceTrajectory*>(rt);

  // Loop over cdc hits
  for (unsigned int i=0;i<cdchits.size();i++){
    my_rt->DistToRT(cdchits[i]->wire);
    my_rt->GetLastDOCAPoint(pos, mom);

    // Create the dE,dx pair from the position and momentum using a helical approximation for the path 
    // in the straw and keep track of the momentum in the active region of the detector
    if (CalcdEdxHit(mom,pos,cdchits[i],de_and_dx)==NOERROR){
      dEdx_list.push_back(de_and_dx);
      
      p_avg+=mom.Mag();
    }
  }
  
  //Get the list of fdc hits used in the fit
  vector<const DFDCPseudo*>fdchits;
  fdchits=GetFDCFitHits();

  // loop over fdc hits
  for (unsigned int i=0;i<fdchits.size();i++){
    my_rt->DistToRT(fdchits[i]->wire);
    my_rt->GetLastDOCAPoint(pos, mom);
   
    pair<double,double>de_and_dx;
    de_and_dx.first=1000.*fdchits[i]->dE; // MeV
    //double gas_density=0.0018; // g/cm^3
    double gas_thickness=1.0; // cm
    //de_and_dx.second=gas_density*gas_thickness/cos(mom.Theta());// g/cm^2
    de_and_dx.second=gas_thickness/cos(mom.Theta());
  }
    
  // Sort the dEdx entries from smallest to largest
  sort(dEdx_list.begin(),dEdx_list.end(),DTrackFitter_dedx_cmp);  

  // Compute the dEdx in the active volume for the track using a truncated 
  // mean to minimize the effect of the long Landau tail
  //num_hits=(dEdx_list.size()>5)?int(0.6*dEdx_list.size()):dEdx_list.size();
  num_hits=int(0.5*dEdx_list.size());
  if (num_hits>0){    
    double sum_dE=0.;
    double sum_ds=0.;
    for (unsigned int i=0;i<num_hits;i++){
      sum_ds+=dEdx_list[i].second;
      sum_dE+=dEdx_list[i].first; 
    }
    dedx=sum_dE/sum_ds;// MeV/cm
    mean_path_length=sum_ds/double(num_hits);    
    p_avg/=double(dEdx_list.size());

    return NOERROR;
  }
 
  return RESOURCE_UNAVAILABLE;
}

// Compute the dEdx for a track described by the reference trajectory rt 
// using a subset of all the FDC and CDC hits on the track. Also returns 
// the mean path length per hit in the active volume of the detector and 
// the average measured momentum within the active region.
jerror_t DTrackFitter::GetdEdx(const DReferenceTrajectory *rt,
			       vector<dedx_t>&dEdx_list){
  DVector3 pos,mom;
  //Get the list of cdc hits used in the fit
  vector<const DCDCTrackHit*>cdchits;
  cdchits=GetCDCFitHits();

  //dE and dx pairs
  pair<double,double>de_and_dx;

  // We cast away the const-ness of the reference trajectory so that we can use the DisToRT method
  DReferenceTrajectory *my_rt=const_cast<DReferenceTrajectory*>(rt);

  // Loop over cdc hits
  for (unsigned int i=0;i<cdchits.size();i++){
    my_rt->DistToRT(cdchits[i]->wire);
    my_rt->GetLastDOCAPoint(pos, mom);

    // Create the dE,dx pair from the position and momentum using a helical approximation for the path 
    // in the straw and keep track of the momentum in the active region of the detector
    if (CalcdEdxHit(mom,pos,cdchits[i],de_and_dx)==NOERROR){
      dEdx_list.push_back(dedx_t(de_and_dx.first,de_and_dx.second,
				 mom.Mag()));
    }
  }
  
  //Get the list of fdc hits used in the fit
  vector<const DFDCPseudo*>fdchits;
  fdchits=GetFDCFitHits();

  // loop over fdc hits
  for (unsigned int i=0;i<fdchits.size();i++){
    my_rt->DistToRT(fdchits[i]->wire);
    my_rt->GetLastDOCAPoint(pos, mom);
   
    double gas_thickness=1.0; // cm
    dEdx_list.push_back(dedx_t(1000.*fdchits[i]->dE,
			       gas_thickness/cos(mom.Theta()),
			       mom.Mag()));
  }
    
  // Sort the dEdx entries from smallest to largest
  sort(dEdx_list.begin(),dEdx_list.end(),DTrackFitter_dedx_cmp2);  
 
  return NOERROR;
}
