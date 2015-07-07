
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

#include "AMPTOOLS_AMPS/Pi0Regge.h"
#include "AMPTOOLS_AMPS/Pi0SAID.h"

#include "AMPTOOLS_MCGEN/GammaPToXP.h"

#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/ConfigFileParser.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

using std::complex;
using namespace std;

int main( int argc, char* argv[] ){
  
	// random number initialization - this is not GlueX standard and
	// should be standardized in the future
	
	srand48( time( NULL ) );
	
	string  configfile("");
	string  outname("");
	string  hddmname("");
	
	bool diag = false;
	bool genFlat = false;
	
	double beamMaxE   = 5.5;
	double beamPeakE  = 3.0;
	double beamLowE   = 1.5;
	double beamHighE  = 5.5;
	
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
		if (arg == "-m"){  
			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
			else  beamMaxE = atof( argv[++i] ); }
		if (arg == "-p"){  
			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
			else  beamPeakE = atof( argv[++i] ); }
		if (arg == "-a"){  
			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
			else  beamLowE = atof( argv[++i] ); }
		if (arg == "-b"){  
			if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
			else  beamHighE = atof( argv[++i] ); }
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
			cout << "\t -m  <value>\t Electron beam energy (or photon energy endpoint) [optional]" << endl;
			cout << "\t -p  <value>\t Coherent peak photon energy [optional]" << endl;
			cout << "\t -a  <value>\t Minimum photon energy to simulate events [optional]" << endl;
			cout << "\t -b  <value>\t Maximum photon energy to simulate events [optional]" << endl;
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
	
	// setup AmpToolsInterface
	AmpToolsInterface::registerAmplitude( Pi0Regge() );
	AmpToolsInterface::registerAmplitude( Pi0SAID() );
	AmpToolsInterface ati( cfgInfo, AmpToolsInterface::kMCGeneration );
	
	// generate single pi0 production
	GammaPToXP phasespace( 0.135, beamMaxE, beamPeakE, beamLowE, beamHighE); 
	
	vector< int > pTypes;
	pTypes.push_back( Gamma );
	pTypes.push_back( Proton );
	pTypes.push_back( Pi0 );
	
	HDDMDataWriter* hddmOut = NULL;
	if( hddmname.size() != 0 ) hddmOut = new HDDMDataWriter( hddmname, 9306 );
	ROOTDataWriter rootOut( outname );
	
	TFile* diagOut = new TFile( "gen_ppi0_diagnostic.root", "recreate" );
	TH2F* hCosTheta_phi = new TH2F( "CosTheta_phi", "cos#theta vs. #phi; #phi; cos#theta", 180, -3.14, 3.14, 100, -1, 1);
	TH2F* ht_phi = new TH2F( "t_phi", "-t vs. #phi; #phi; -t (GeV^{2})", 100, -3.14, 3.14, 100, 0, 2);
	
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
				double rand = drand48() * maxInten;
				
				if( weightedInten > rand || genFlat ){
					
					// calculate angular variables
					TLorentzVector target  ( 0., 0., 0., 0.938);	
					TLorentzVector beam = evt->particle ( 0 );
					TLorentzVector recoil = evt->particle ( 1 );
					TLorentzVector p1 = evt->particle ( 2 );
					
					TLorentzVector cm = recoil + p1;
					TLorentzRotation cmBoost( -cm.BoostVector() );
					
					TLorentzVector recoil_cm = cmBoost * recoil;
					TLorentzVector p1_cm = cmBoost * p1;
					
					GDouble t = (target - recoil).M2();
					GDouble CosTheta = p1_cm.CosTheta();
					GDouble phi = p1_cm.Phi();
					if(phi < -1*PI) phi += 2*PI;
					if(phi > PI) phi -= 2*PI;
					
					hCosTheta_phi->Fill( phi, CosTheta);
					ht_phi->Fill( phi, -1.*t);
					
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
	
	hCosTheta_phi->Write();
	ht_phi->Write();
	diagOut->Close();
	
	if( hddmOut ) delete hddmOut;
	
	return 0;
}


