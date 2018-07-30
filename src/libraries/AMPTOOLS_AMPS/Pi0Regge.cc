
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/Pi0Regge.h"

#include <AMPTOOLS_MCGEN/CobremsGeneration.hh>

Pi0Regge::Pi0Regge( const vector< string >& args ) :
UserAmplitude< Pi0Regge >( args )
{
	assert( args.size() == 1 );
	// Polarization plane angle (PARA = 0 and PERP = PI/2)
	PolPlane = atof( args[0].c_str() );

	// Initialize coherent brem table
	// Do this over the full range since we will be using this as a lookup
	float Emax  = 12.0;
	float Epeak = 9.0;
	float Elow  = 0.135;
	float Ehigh = 12.0;
	
	int doPolFlux=0;  // want total flux (1 for polarized flux)
	float emitmr=10.e-9; // electron beam emittance
	float radt=50.e-6; // radiator thickness in m
	float collDiam=0.005; // meters
	float Dist = 76.0; // meters
	CobremsGeneration cobrems(Emax, Epeak);
	cobrems.setBeamEmittance(emitmr);
	cobrems.setTargetThickness(radt);
	cobrems.setCollimatorDistance(Dist);
	cobrems.setCollimatorDiameter(collDiam);
	cobrems.setCollimatedFlag(true);
	cobrems.setPolarizedFlag(doPolFlux);
	
	// Create histogram
	totalFlux_vs_E = new TH1D("totalFlux_vs_E", "Total Flux vs. E_{#gamma}", 1000, Elow, Ehigh);
	polFlux_vs_E   = new TH1D("polFlux_vs_E", "Polarized Flux vs. E_{#gamma}", 1000, Elow, Ehigh);
	polFrac_vs_E   = new TH1D("polFrac_vs_E", "Polarization Fraction vs. E_{#gamma}", 1000, Elow, Ehigh);

	// Fill totalFlux
	for(int i=1;i<=totalFlux_vs_E->GetNbinsX(); i++){
		double x = totalFlux_vs_E->GetBinCenter(i)/Emax;
		double y = 0;
		//if(Epeak<Elow) y = cobrems.Rate_dNidx(x);
		y = cobrems.Rate_dNtdx(x);
		totalFlux_vs_E->SetBinContent(i, y);
	}

	doPolFlux=1;
	cobrems.setPolarizedFlag(doPolFlux);
	// Fill totalFlux
	for(int i=1;i<=polFlux_vs_E->GetNbinsX(); i++){
		double x = polFlux_vs_E->GetBinCenter(i)/Emax;
		double y = 0;
		//if(Epeak<Elow) y = cobrems.Rate_dNidx(x);
		y = cobrems.Rate_dNcdx(x);
		polFlux_vs_E->SetBinContent(i, y);
	}
	
	polFrac_vs_E->Divide(polFlux_vs_E, totalFlux_vs_E);
}


complex< GDouble >
Pi0Regge::calcAmplitude( GDouble** pKin ) const {
  
	TLorentzVector target  ( 0., 0., 0., 0.938);
	TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] ); 
	TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 
	TLorentzVector p1     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] ); 
	
	TLorentzVector cm = recoil + p1;
	TLorentzRotation cmBoost( -cm.BoostVector() );
	
	TLorentzVector beam_cm = cmBoost * beam;
	TLorentzVector target_cm = cmBoost * target;
	TLorentzVector recoil_cm = cmBoost * recoil;
	
	// phi dependence needed for polarized distribution
	TLorentzVector p1_cm = cmBoost * p1;
	GDouble phi = p1_cm.Phi() + PolPlane*TMath::Pi()/180.;
	GDouble cos2Phi = cos(2.*phi);
	
	// polarization from cobrem.F
	int bin = polFrac_vs_E->GetXaxis()->FindBin(beam.E());
	GDouble Pgamma;
	if (bin == 0 || bin > polFrac_vs_E->GetXaxis()->GetNbins()){
		Pgamma = 0.;
	}
	else Pgamma = polFrac_vs_E->GetBinContent(bin);

	// factors needed to calculate amplitude in c++ code
	GDouble Ecom = cm.M();
	GDouble theta = p1_cm.Theta();	

	// amplitude coded in c++ (include calculation of beam asymmetry)
	GDouble BeamSigma = 0.;
	GDouble W = Pi0PhotCS_S(Ecom, theta, BeamSigma);
	W *= (1 - Pgamma * BeamSigma * cos2Phi);

	return complex< GDouble > ( sqrt( fabs(W) ) );
}

