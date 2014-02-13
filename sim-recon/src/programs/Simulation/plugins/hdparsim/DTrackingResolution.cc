// $Id$
//
//    File: DTrackingResolution.cc
// Created: Mon Feb 25 15:06:17 EST 2008
// Creator: davidl (on Darwin fwing-dhcp13.jlab.org 8.11.1 i386)
//

#include <cmath>
#include <iostream>
using namespace std;

#include "DTrackingResolution.h"
#include "DFactoryGeneratorHDParSim.h"

// Routine used If we're a plugin
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddFactoryGenerator(new DFactoryGeneratorHDParSim());
}
} // "C"


//---------------------------------
// DTrackingResolution    (Constructor)
//---------------------------------
DTrackingResolution::DTrackingResolution()
{
	scale_err_pt = 1.0;
	scale_err_theta = 1.0;
	scale_err_phi = 1.0;
}

//---------------------------------
// ~DTrackingResolution    (Destructor)
//---------------------------------
DTrackingResolution::~DTrackingResolution()
{

}

//----------------
// SetErrorScaleFactors
//----------------
void DTrackingResolution::SetErrorScaleFactors(double scale_err_pt, double scale_err_theta, double scale_err_phi)
{
	this->scale_err_pt = scale_err_pt;
	this->scale_err_theta =scale_err_theta;
	this->scale_err_phi = scale_err_phi;
}

//----------------
// GetErrorScaleFactors
//----------------
void DTrackingResolution::GetErrorScaleFactors(double &scale_err_pt, double &scale_err_theta, double &scale_err_phi)
{
	scale_err_pt = this->scale_err_pt;
	scale_err_theta = this->scale_err_theta;
	scale_err_phi = this->scale_err_phi;
}

//----------------
// Smear
//----------------
bool DTrackingResolution::Smear(int geanttype, TVector3 &mom)
{
	/// Smear the momentum vector of a charged particle
	/// based on reolutions obtained from the GetResolution
	/// method.
	///
	/// The value of geanttype should specify the particle
	/// type using the GEANT particle ids.
	///
	/// The units of mom should be GeV/c

	// Efficiency should be based on input values. Calculate that
	// now before they are smeared.
	bool is_reconstructed = Efficiency(geanttype, mom);

	// Get resolutions
	double pt_res, theta_res, phi_res;
	GetResolution(geanttype, mom, pt_res, theta_res, phi_res);
	
	// Scale resolutions (default scale factors are all 1.0)
	pt_res *= scale_err_pt;
	theta_res *= scale_err_theta;
	phi_res *= scale_err_phi;
	
	// Calculate new values. For each, we make a check that it is
	// still in a valid range (e.g. theta is not less than 0).
	// In reality, the probablity functions near these hard limits
	// would not be gaussian in shape and so should be handled
	// quite differently.
	double theta_new=-1.0, phi_new=-1.0;
	
	if(theta_res>0.0){
		while(theta_new<=0.0 || theta_new>M_PI)theta_new = mom.Theta() + rnd.Gaus(0.0, theta_res)/1000.0;
	}else{
		theta_new = mom.Theta();
	}
	phi_new = mom.Phi() + rnd.Gaus(0.0, phi_res)/1000.0;
	while(phi_new<-M_PI)phi_new+=M_PI;
	while(phi_new>=M_PI)phi_new-=M_PI;
	
	// Overwrite input vector with new values.
	// For photons, the value returned in "pt_res" is actually the
	// total energy resolution.
	if(geanttype==1){
		// photons
		double E_res = pt_res, E_new=-1.0;
		if(mom.Mag()>0.0){
			while(E_new<=0.0) E_new = mom.Mag()*(1.0 + rnd.Gaus(0.0, E_res));
		}else{
			E_new = 0.0;
		}
		mom.SetMagThetaPhi(E_new, theta_new, phi_new);
	}else{
		// not photons
		double pt_new=-1.0;
		if(mom.Perp()>0.0){
			while(pt_new<=0.0) pt_new = mom.Perp()*(1.0 + rnd.Gaus(0.0, pt_res));
		}else{
			pt_new = 0.0;
		}
		mom.SetMagThetaPhi(pt_new/sin(theta_new), theta_new, phi_new);
	}
	return is_reconstructed;
}

//----------------
// Efficiency
//----------------
bool DTrackingResolution::Efficiency(int geanttype, const TVector3 &mom)
{
	/// Return a boolean saying whether this event would be reconstructed
	/// or not. The value returned will vary in that if mom points to
	/// an area of the detector with 80% efficienct, then this will return
	/// "true" 80% of the time and "false" 20% of the time.
	///
	/// Geometric acceptance is also included so values of mom pointing
	/// away from the detector will always returen "false".
	///
	/// This works by calling the GetEfficiency method and then picking a 
	/// random number between 0 and 1. If the random number is less than or
	/// equal to the efficiency value, then true is returned . Otherwise,
	/// false is returned.

	double eff = GetEfficiency(geanttype, mom);
	
	double s = rnd.Rndm();
	
	return s<=eff;
}



