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
	bool isParsed = parseConfig();
	if(!isParsed) exit(1);

	// Fill beam property histograms based on config file
	if(mIsROOT) 
		fillFromROOT();
	else if(mIsCCDB) 
		fillFromCCDB();
	else 
		generateCobrems(); // default to CobremsGeneration
		
}

bool BeamProperties::parseConfig(){

	//cout<<endl<<"Parsing BeamProperty config file: "<<mConfigFile.Data()<<endl<<endl;

	// start assuming parameters for CobremsGeneration
	mIsROOT = false;
	mIsCCDB = false;

	ifstream inputFile;
	inputFile.open(mConfigFile.Data());
	if (!inputFile.is_open()){
		cout << "BeamProperties ERROR:  Could not open configuration file: " << mConfigFile << endl;
		return false;
	}

	while(inputFile) { 
		
		string line;
		getline(inputFile,line);

		// parse the line into words (from AmpTools ConfigFileParser)
		vector<string> words;
		string word("");
		for (unsigned int j = 0; j < line.size(); j++){
			if (!isspace(line[j])){
				word += line[j];
				if ((j == (line.size()-1))&&(!word.empty())){
					words.push_back(word);
					word = "";
				}
			}
			else if (!word.empty()){
				words.push_back(word);
				word = "";
			}
		}	

		// skip blank or incomplete lines and comments
		if(words.size() < 2 || words[0][0] == '#') 
			continue;

		// 1st is variable name and 2nd word is value
		if(words[0].find("ROOT") != std::string::npos) {
			mIsROOT = true;
			mBeamHistNameMap.insert( std::pair<std::string,std::string>( words[0].data(), words[1].data() ) );
		}
		else 
			mBeamParametersMap.insert( std::pair<std::string,double>( words[0].data(), atof(words[1].data()) ));
	}
	inputFile.close();

	return true;
}

// create histograms for flux and polarization fraction using CobremsGeneration
void BeamProperties::generateCobrems(){

	// check that required parameters are provided
	const int nParameters = 8;
	string parameterNames[nParameters] = {"ElectronBeamEnergy", "CoherentPeakEnergy", "PhotonBeamLowEnergy", "PhotonBeamHighEnergy",  "Emittance", "RadiatorThickness", "CollimatorDiameter", "CollimatorDistance"};
	for(int i=0; i<nParameters; i++) {
		if(!mBeamParametersMap.count(parameterNames[i].data())) {
			cout << "BeamProperties ERROR:  generateCombrems parameter " << parameterNames[i] << " missing" << endl;
			exit(1);
		}
	}

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

	// open files and check that they exist
	TFile *fFlux, *fPol;
	if(mBeamHistNameMap.count("ROOTFluxFile") && mBeamHistNameMap.count("ROOTPolFile")) {
		fFlux = TFile::Open(mBeamHistNameMap.at("ROOTFluxFile").data());
		fPol = TFile::Open(mBeamHistNameMap.at("ROOTPolFile").data());
	}
	else {
		cout << "BeamProperties ERROR:  ROOT flux or polarization file name not defined in configuration file" << endl;
		exit(1);
	}
	if(!fFlux->IsOpen() || !fPol->IsOpen()) {
		cout << "BeamProperties ERROR:  Could not open ROOT flux " << mBeamHistNameMap.at("ROOTFluxFile").data() << " or polarization " << mBeamHistNameMap.at("ROOTPolFile").data() << " files don't exist" << endl;
		exit(1);
	}
	
	// open histograms and check that they exist
	if(mBeamHistNameMap.count("ROOTFluxName") && mBeamHistNameMap.count("ROOTPolName")) {
		fluxVsEgamma = (TH1D*)fFlux->Get(mBeamHistNameMap.at("ROOTFluxName").data())->Clone("BeamProperties_FluxVsEgamma");
		polFracVsEgamma = (TH1D*)fPol->Get(mBeamHistNameMap.at("ROOTPolName").data())->Clone("BeamProperties_PolFracVsEgamma");
	}
	else {
		cout << "BeamProperties ERROR:  ROOT flux or polarization histogram name not defined in configuration file" << endl;
		exit(1);
	}
	if(!fluxVsEgamma || !polFracVsEgamma) {
		cout << "BeamProperties ERROR:  ROOT flux " << mBeamHistNameMap.at("ROOTFluxFile").data() << " or polarization " << mBeamHistNameMap.at("ROOTPolFile").data() << " histograms don't exist" << endl;
		exit(1);
	}
	
	// set energy range for event generation
	if(!mBeamParametersMap.count("PhotonBeamLowEnergy") || !mBeamParametersMap.count("PhotonBeamHighEnergy")) {
		cout << "BeamProperties ERROR:  PhotonBeamLowEnergy or PhotonBeamHighEnergy not specified for event generation" << endl;
		exit(1);
	}
	fluxVsEgamma->GetXaxis()->SetRangeUser(mBeamParametersMap.at("PhotonBeamLowEnergy"), mBeamParametersMap.at("PhotonBeamHighEnergy")); 

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
		cout<<endl<<"WARNING: Polarization angle requeseted by generator/fitter but not found in beam configuration file.  Using PolAngle = 0 degrees by default"<<endl<<endl;

	return 0.;
}
