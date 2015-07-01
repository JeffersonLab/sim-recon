
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

#include "AMPTOOLS_AMPS/TwoPiAngles.h"
#include "AMPTOOLS_AMPS/BreitWigner.h"

#include "AMPTOOLS_MCGEN/ProductionMechanism.h"
#include "AMPTOOLS_MCGEN/GammaPToXYP.h"

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
  
  // default upper and lower bounds 
  double lowMass = 0.2;
  double highMass = 2.0;

  double beamE = 3.0;
  
  int nEvents = 1000;
  int batchSize = 1000;
    
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
    if (arg == "-l"){
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  lowMass = atof( argv[++i] ); }
    if (arg == "-u"){
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  highMass = atof( argv[++i] ); }
    if (arg == "-n"){  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  nEvents = atoi( argv[++i] ); }
    if (arg == "-d"){
      diag = true; }
    if (arg == "-f"){
      genFlat = true; }
    if (arg == "-h"){
      cout << endl << " Usage for: " << argv[0] << endl << endl;
      cout << "\t -c  <file>\t Config file" << endl;
      cout << "\t -o  <name>\t ROOT file output name" << endl;
      cout << "\t -hd <name>\t HDDM file output name [optional]" << endl;
      cout << "\t -l  <value>\t Low edge of mass range (GeV) [optional]" << endl;
      cout << "\t -u  <value>\t Upper edge of mass range (GeV) [optional]" << endl;
      cout << "\t -n  <value>\t Minimum number of events to generate [optional]" << endl;
      cout << "\t -f \t\t Generate flat in M(X) (no physics) [optional]" << endl;
      cout << "\t -d \t\t Plot only diagnostic histograms [optional]" << endl << endl;
      exit(1);
    }
  }
  
  if( configfile.size() == 0 || outname.size() == 0 ){
    cout << "No config file or output specificed:  run gen_2pi -h for help" << endl;
    exit(1);
  }
  
  // open config file and be sure only one reaction is specified
  ConfigFileParser parser( configfile );
  ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
  assert( cfgInfo->reactionList().size() == 1 );
  ReactionInfo* reaction = cfgInfo->reactionList()[0];
  
  // setup AmpToolsInterface
  AmpToolsInterface::registerAmplitude( TwoPiAngles() );
  AmpToolsInterface::registerAmplitude( BreitWigner() );
  AmpToolsInterface ati( cfgInfo, AmpToolsInterface::kMCGeneration );
  
  ProductionMechanism::Type type =
    ( genFlat ? ProductionMechanism::kFlat : ProductionMechanism::kResonant );
   
  // generate over a range of mass -- the daughters are two charged pions
  GammaPToXYP resProd( lowMass, highMass, 0.140, 0.140, beamE, type );
  
  // seed the distribution with a sum of noninterfering Breit-Wigners
  // we can easily compute the PDF for this and divide by that when
  // doing accept/reject -- improves efficiency if seeds are picked well

  if( !genFlat ){
     
    // the lines below should be tailored by the user for the particular desired
    // set of amplitudes -- doing so will improve efficiency.  Leaving as is
    // won't make MC incorrect, it just won't be as fast as it could be
    
    resProd.addResonance( 0.775, 0.146,  1.0 );
  }
  
  vector< int > pTypes;
  pTypes.push_back( Gamma );
  pTypes.push_back( Proton );
  pTypes.push_back( PiPlus );
  pTypes.push_back( PiMinus );
  
  HDDMDataWriter* hddmOut = NULL;
  if( hddmname.size() != 0 ) hddmOut = new HDDMDataWriter( hddmname, 9306 );
  ROOTDataWriter rootOut( outname );
  
  TFile* diagOut = new TFile( "gen_2pi_diagnostic.root", "recreate" );
  
  TH1F* mass = new TH1F( "M", "Resonance Mass", 180, lowMass, highMass );
  TH1F* massW = new TH1F( "M_W", "Weighted Resonance Mass", 180, lowMass, highMass );
  massW->Sumw2();
  TH1F* intenW = new TH1F( "intenW", "True PDF / Gen. PDF", 1000, 0, 100 );
  TH2F* intenWVsM = new TH2F( "intenWVsM", "Ratio vs. M", 100, lowMass, highMass, 1000, 0, 10 );
 
  TH2F* CosTheta_psi = new TH2F( "CosTheta_psi", "cos#theta vs. #psi", 180, -3.14, 3.14, 100, -1, 1);

  //TH2F* dalitz = new TH2F( "dalitz", "Dalitz Plot", 100, 0, 3.0, 100, 0, 3.0 );
  
  int eventCounter = 0;
  while( eventCounter < nEvents ){
    
    if( batchSize < 1E4 ){
      
      cout << "WARNING:  small batches could have batch-to-batch variations\n"
      << "          due to different maximum intensities!" << endl;
    }
    
    cout << "Generating four-vectors..." << endl;
    
    ati.clearEvents();
    for( int i = 0; i < batchSize; ++i ){
      
      Kinematics* kin = resProd.generate();
      ati.loadEvent( kin, i, batchSize );
      delete kin;
    }
    
    cout << "Processing events..." << endl;

    // include factor of 1.5 to be safe in case we miss peak -- avoid
    // intensity calculation of we are generating flat data
    double maxInten = ( genFlat ? 1 : 1.5 * ati.processEvents( reaction->reactionName() ) );
    
    
    for( int i = 0; i < batchSize; ++i ){
  		
      Kinematics* evt = ati.kinematics( i );
      TLorentzVector resonance( evt->particle( 2 ) + 
                                  evt->particle( 3 ) );

      double genWeight = evt->weight();
      
      // cannot ask for the intensity if we haven't called process events above
      double weightedInten = ( genFlat ? 1 : ati.intensity( i ) );
      
      if( !diag ){
        
        // obtain this by looking at the maximum value of intensity * genWeight
        double rand = drand48() * maxInten;
        
        if( weightedInten > rand || genFlat ){
          
          mass->Fill( resonance.M() );
          massW->Fill( resonance.M(), genWeight );
          
          intenW->Fill( weightedInten );
          intenWVsM->Fill( resonance.M(), weightedInten );

	  // calculate angular variables
	  TLorentzVector beam = evt->particle ( 0 );
	  TLorentzVector recoil = evt->particle ( 1 );
	  TLorentzVector p1 = evt->particle ( 2 );

	  TLorentzRotation resonanceBoost( -resonance.BoostVector() );

	  TLorentzVector beam_res = resonanceBoost * beam;
	  TLorentzVector recoil_res = resonanceBoost * recoil;
	  TLorentzVector p1_res = resonanceBoost * p1;
	  
	  TVector3 z = -recoil_res.Vect().Unit();
	  TVector3 y = beam_res.Vect().Cross(z).Unit();
	  TVector3 x = y.Cross(z).Unit();
	  
	  TVector3 angles(   (p1_res.Vect()).Dot(x),
			       (p1_res.Vect()).Dot(y),
			       (p1_res.Vect()).Dot(z) );
	  
	  GDouble CosTheta = angles.CosTheta();
	 
	  GDouble phi = angles.Phi();
	  GDouble Phi = -1.*resonance.Vect().Phi();
	  
	  GDouble psi = phi - Phi;
	  if(psi < -1*PI) psi += 2*PI;
	  if(psi > PI) psi -= 2*PI;

	  CosTheta_psi->Fill( psi, CosTheta);
           
          // we want to save events with weight 1
          evt->setWeight( 1.0 );
          
          if( hddmOut ) hddmOut->writeEvent( *evt, pTypes );
          rootOut.writeEvent( *evt );
          ++eventCounter;
        }
      }
      else{
        
        mass->Fill( resonance.M() );
        massW->Fill( resonance.M(), genWeight );

        intenW->Fill( weightedInten );
        intenWVsM->Fill( resonance.M(), weightedInten );
	TLorentzVector recoil = evt->particle ( 1 );
 
        ++eventCounter;
      }
      
      delete evt;
    }
    
    cout << eventCounter << " events were processed." << endl;
  }
  
  mass->Write();
  massW->Write();
  intenW->Write();
  intenWVsM->Write();
  CosTheta_psi->Write();
  diagOut->Close();
  
  if( hddmOut ) delete hddmOut;
  
	return 0;
}


