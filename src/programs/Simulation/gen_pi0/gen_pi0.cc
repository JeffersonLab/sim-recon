
#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include <cassert>
#include <cstdlib>

#include "particleType.h"

#include "AMPTOOLS_DATAIO/ROOTDataWriter.h"
#include "AMPTOOLS_DATAIO/HDDMDataWriter.h"
#include "AMPTOOLS_DATAIO/Pi0PlotGenerator.h"

#include "AMPTOOLS_AMPS/Pi0Regge.h"
#include "AMPTOOLS_AMPS/Pi0SAID.h"

#include "AMPTOOLS_MCGEN/GammaPToXP.h"

#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/PlotGenerator.h"
#include "IUAmpTools/FitResults.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TRandom3.h"

using std::complex;
using namespace std;

int main( int argc, char* argv[] ){
  
	string  configfile("");
	string  outname("");
	string  hddmname("");
	
	bool diag = false;
	bool genFlat = false;
	
	int runNum = 9001;
	int seed = 0;

	int nEvents = 10000;
	int batchSize = 10000;
	
	//parse command line:
	for (int i = 1; i < argc; i++){
		
		string arg(argv[i]);
		
		if (arg == "-c"){  
			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
			else  configfile = argv[++i]; }
		if (arg == "-o"){  
			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
			else  outname = argv[++i]; }
		if (arg == "-hd"){
			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
			else  hddmname = argv[++i]; }
		if (arg == "-n"){  
			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
			else  nEvents = atoi( argv[++i] ); }
		if (arg == "-r"){
                        if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
                        else  runNum = atoi( argv[++i] ); }
		if (arg == "-s"){
                        if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
                        else  seed = atoi( argv[++i] ); }
		if (arg == "-d"){
			diag = true; }
		if (arg == "-f"){
			genFlat = true; }
		if (arg == "-h"){
			cout << endl << " Usage for: " << argv[0] << endl << endl;
			cout << "\t -c  <file>\t Config file" << endl;
			cout << "\t -o  <name>\t ROOT file output name" << endl;
			cout << "\t -hd <name>\t HDDM file output name [optional]" << endl;
			cout << "\t -n  <value>\t Minimum number of events to generate [optional]" << endl;
			cout << "\t -r  <value>\t Run number assigned to generated events [optional]" << endl;
			cout << "\t -s  <value>\t Random number seed initialization [optional]" << endl;
			cout << "\t -f \t\t Generate flat in M(X) (no physics) [optional]" << endl;
			cout << "\t -d \t\t Plot only diagnostic histograms [optional]" << endl << endl;
			exit(1);
		}
	}
	
	if( configfile.size() == 0 || outname.size() == 0 ){
		cout << "No config file or output specificed:  run gen_pi0 -h for help" << endl;
		exit(1);
	}
	
	// open config file and be sure only one reaction is specified
	ConfigFileParser parser( configfile );
	ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
	assert( cfgInfo->reactionList().size() == 1 );
	ReactionInfo* reaction = cfgInfo->reactionList()[0];
	
	// random number initialization (set to 0 by default)
	gRandom->SetSeed(seed);

	// setup AmpToolsInterface
	AmpToolsInterface::registerAmplitude( Pi0Regge() );
	AmpToolsInterface::registerAmplitude( Pi0SAID() );
	AmpToolsInterface ati( cfgInfo, AmpToolsInterface::kMCGeneration );
	
	// loop to look for beam configuration file
	TString beamConfigFile;
	const vector<ConfigFileLine> configFileLines = parser.getConfigFileLines();
	for (vector<ConfigFileLine>::const_iterator it=configFileLines.begin(); it!=configFileLines.end(); it++) {
		if ((*it).keyword() == "define") {
			TString beamArgument =  (*it).arguments()[0].c_str();
			if(beamArgument.Contains("beam")) {
				beamConfigFile = (*it).arguments()[1].c_str();
				//cout<<beamConfigFile.Data()<<endl;
			}
		}
	}

	// generate single pi0 production
	GammaPToXP phasespace( 0.135, beamConfigFile);
	
	vector< int > pTypes;
	pTypes.push_back( Gamma );
	pTypes.push_back( Proton );
	pTypes.push_back( Pi0 );
	
	HDDMDataWriter* hddmOut = NULL;
	if( hddmname.size() != 0 ) hddmOut = new HDDMDataWriter( hddmname, runNum );
	ROOTDataWriter rootOut( outname );
	
	// use plot generator for diagnostic histograms
	TFile* diagOut = new TFile( "gen_diagnostic.root", "recreate" );
	Pi0PlotGenerator plotGen;
	
	int eventCounter = 0;
	while( eventCounter < nEvents ){
		
		if( batchSize < 1E4 ){
			
			cout << "WARNING:  small batches could have batch-to-batch variations\n"
			     << "          due to different maximum intensities!" << endl;
		}
		
		cout << "Generating four-vectors..." << endl;
		
		ati.clearEvents();
		for( int i = 0; i < batchSize; ++i ){
			
			Kinematics* kin = phasespace.generate();
			ati.loadEvent( kin, i, batchSize );
			delete kin;
		}
		
		cout << "Processing events..." << endl;
		
		// include factor of 1.5 to be safe in case we miss peak -- avoid
		// intensity calculation of we are generating flat data
		double maxInten = ( genFlat ? 1 : 1.5 * ati.processEvents( reaction->reactionName() ) );
		
		for( int i = 0; i < batchSize; ++i ){
			
			Kinematics* evt = ati.kinematics( i );
			
			// cannot ask for the intensity if we haven't called process events above
			double weightedInten = ( genFlat ? 1 : ati.intensity( i ) );
			
			if( !diag ){
				
				// obtain this by looking at the maximum value of intensity * genWeight
				double rand = gRandom->Uniform() * maxInten;
				
				if( weightedInten > rand || genFlat ){
					
					// fill PlotGenerator histograms
					plotGen.projectEvent(evt);
					
					// we want to save events with weight 1
					evt->setWeight( 1.0 );
					
					if( hddmOut ) hddmOut->writeEvent( *evt, pTypes );
					rootOut.writeEvent( *evt );
					++eventCounter;
				}
			}
			else{
				
				++eventCounter;
			}
      
			delete evt;
		}
		
		cout << eventCounter << " events were processed." << endl;
	}
	
	// write PlotGenerator histograms to file
	for(int iHist=0; iHist<plotGen.kNumHists; iHist++) {
		Histogram* hist = plotGen.getHistogram(iHist);
		if(hist) {
			TH1* thist = hist->toRoot();
			thist->Write();
		}
	}
	diagOut->Close();
	
	if( hddmOut ) delete hddmOut;
	
	return 0;
}


