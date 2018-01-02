/*
 *  BeamProperties.cc
 *
 *  Contains histograms for beam properties to be used in event generation and fitting.  Source
 *  of beam properties is from CombremsGeneration, external ROOT file or CCDB (to be implemented). 
 *
 *  Created by Justin Stevens on 12/29/17
 */

#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "TROOT.h"
#include "TFile.h"

#include "BeamProperties.h"
#include "CobremsGeneration.hh"

using namespace std;

BeamProperties::BeamProperties( TString configFile ) {

	// check if histograms already exist before re-creating
	gDirectory->cd("/");
	fluxVsEgamma = (TH1D*)gDirectory->Get("BeamPoperties_FluxVsEgamma");
	polFracVsEgamma = (TH1D*)gDirectory->Get("BeamPoperties_PolFracVsEgamma");
	if(!fluxVsEgamma || !polFracVsEgamma) 
		createHistograms(configFile);
}

void BeamProperties::createHistograms( TString configFile ) {

	// First parse configuration file
	mConfigFile = configFile;
	parseConfig();

	// Fill beam property histograms based on config file
	if(mIsROOT) 
		fillFromROOT();
	else if(mIsCCDB) 
		fillFromCCDB();
	else 
		generateCobrems(); // default to CobremsGeneration
		
}

void BeamProperties::parseConfig(){

	cout<<endl<<"Parsing BeamProperty config file: "<<mConfigFile.Data()<<endl<<endl;

	// start assuming parameters for CobremsGeneration
	mIsROOT = false;
	mIsCCDB = false;

	ifstream inputFile;
	inputFile.open(mConfigFile.Data());

	while(true) { 
		if(!inputFile.good()) break;
		
		string line_string;
		getline(inputFile,line_string);

		// skip blank lines and comments
		if(line_string.length() == 0 || line_string[0] == '/') continue;

		TString line = line_string;
		int firstEquals = line.First("=");
		TString varName = TString(line(0,firstEquals));
		int firstDash = line.First("/");
		TString varValue = TString(line(firstEquals+1,firstDash-(firstEquals+1)));
		varValue.ReplaceAll(" ","");

		if(varName.Contains("ROOT")) {
			mIsROOT = true;
			mBeamHistNameMap.insert( std::pair<std::string,std::string>( varName.Data(), varValue.Data() ) );
		}
		else 
			mBeamParametersMap.insert( std::pair<std::string,double>( varName.Data(), atof(varValue.Data()) ) );
	}
	inputFile.close();

	return;
}

// create histograms for flux and polarization fraction using CobremsGeneration
void BeamProperties::generateCobrems(){

	// Set parameters from config file
	double Emax  = mBeamParametersMap.at("ElectronBeamEnergy");
	double Epeak = mBeamParametersMap.at("CoherentPeakEnergy");
	double Elow  = mBeamParametersMap.at("PhotonBeamLowEnergy");
	double Ehigh = mBeamParametersMap.at("PhotonBeamHighEnergy");

	// Create histograms
	int nBinsEgamma = 1000;
	fluxVsEgamma = new TH1D("BeamPoperties_FluxVsEgamma", "Flux vs. E_{#gamma}", nBinsEgamma, Elow, Ehigh);
	polFracVsEgamma = new TH1D("BeamPoperties_PolFracVsEgamma", "Polarization Fraction vs. E_{#gamma}", nBinsEgamma, Elow, Ehigh);
	TH1D *polFluxVsEgamma   = new TH1D("BeamPoperties_PolFluxVsEgamma", "Polarized Flux vs. E_{#gamma}", nBinsEgamma, Elow, Ehigh);
	
	// Setup cobrems
	CobremsGeneration cobrems(Emax, Epeak);
	cobrems.setBeamEmittance(mBeamParametersMap.at("Emittance"));
	cobrems.setTargetThickness(mBeamParametersMap.at("RadiatorThickness"));
	cobrems.setCollimatorDistance(mBeamParametersMap.at("CollimatorDiameter"));
	cobrems.setCollimatorDiameter(mBeamParametersMap.at("CollimatorDistance"));
	cobrems.setCollimatedFlag(true);
	
	// Fill flux
	cobrems.setPolarizedFlag(0); // 0=total flux
	for(int i=1; i<=nBinsEgamma; i++){
		double x = fluxVsEgamma->GetBinCenter(i)/Emax;
		double y = 0;
		if(Epeak<Elow) y = cobrems.Rate_dNidx(x);
		else y = cobrems.Rate_dNtdx(x);
		fluxVsEgamma->SetBinContent(i, y);
	}

	// Fill polarized flux
	cobrems.setPolarizedFlag(1); // 1=polarized flux
	for(int i=1;i<=nBinsEgamma; i++){
		double x = fluxVsEgamma->GetBinCenter(i)/Emax;
		double y = 0;
		if(Epeak<Elow) y = 0.;
		else y = cobrems.Rate_dNcdx(x);
		polFluxVsEgamma->SetBinContent(i, y);
	}
	
	// Polarization fraction from ratio
	polFracVsEgamma->Divide(polFluxVsEgamma, fluxVsEgamma);

	return;
}

// load ROOT histograms for flux and polarization fraction from external file
void BeamProperties::fillFromROOT() {

	TFile *fFlux = TFile::Open(mBeamHistNameMap.at("ROOTFluxFile").data());
	fluxVsEgamma = (TH1D*)fFlux->Get(mBeamHistNameMap.at("ROOTFluxName").data())->Clone("BeamProperties_FluxVsEgamma");
	// set energy range for event generation
	fluxVsEgamma->GetXaxis()->SetRangeUser(mBeamParametersMap.at("PhotonBeamLowEnergy"), mBeamParametersMap.at("PhotonBeamHighEnergy")); 

	TFile *fPol = TFile::Open(mBeamHistNameMap.at("ROOTPolFile").data());
	polFracVsEgamma = (TH1D*)fPol->Get(mBeamHistNameMap.at("ROOTPolName").data())->Clone("BeamProperties_PolFracVsEgamma");

	// keep in memory after file is closed
	fluxVsEgamma->SetDirectory(gROOT);
	polFracVsEgamma->SetDirectory(gROOT);

	fFlux->Close();
	fPol->Close();

        return;
}

// create histograms for flux and polarization fraction from external file for specified run number (not implemented!)
void BeamProperties::fillFromCCDB() {

        return;
}

double BeamProperties::GetPolAngle() {
	if(mBeamParametersMap.count("PolarizationAngle"))
		return mBeamParametersMap.at("PolarizationAngle");
	else 
		cout<<endl<<"Warning: Polarization angle requeseted by generator/fitter but not found in beam configuration file.  Using PolAngle = 0 degrees by default"<<endl<<endl;

	return 0.;
}
