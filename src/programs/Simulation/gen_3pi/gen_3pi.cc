
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

#include "AMPTOOLS_AMPS/ThreePiAngles.h"
#include "AMPTOOLS_AMPS/BreitWigner.h"

#include "AMPTOOLS_MCGEN/ProductionMechanism.h"
#include "AMPTOOLS_MCGEN/GammaPToXYZP.h"

#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "CLHEP/Vector/LorentzVector.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

using std::complex;
using namespace std;
using namespace CLHEP;

int main( int argc, char* argv[] ){
  
  string  configfile("");
  string  outname("");
  string  hddmname("");
  
  bool diag = false;
  bool genFlat = false;
  
  // default upper and lower bounds 
  double lowMass = 0.7;
  double highMass = 2.0;
  
  int nEvents = 100000;
  int batchSize = 100000;
    
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
    cout << "No config file or output specificed:  run gen_3pi -h for help" << endl;
    exit(1);
  }
  
  // open config file and be sure only one reaction is specified
  ConfigFileParser parser( configfile );
  ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
  assert( cfgInfo->reactionList().size() == 1 );
  ReactionInfo* reaction = cfgInfo->reactionList()[0];
  
  // setup AmpToolsInterface
  AmpToolsInterface::registerAmplitude( ThreePiAngles() );
  AmpToolsInterface::registerAmplitude( BreitWigner() );
  AmpToolsInterface ati( cfgInfo, AmpToolsInterface::kMCGeneration );
  
  ProductionMechanism::Type type =
    ( genFlat ? ProductionMechanism::kFlat : ProductionMechanism::kResonant );
   
  // generate over a range of mass -- the daughters are three charged pions
  GammaPToXYZP resProd( lowMass, highMass, 0.140, 0.140, 0.140, type );
  
  // seed the distribution with a sum of noninterfering Breit-Wigners
  // we can easily compute the PDF for this and divide by that when
  // doing accept/reject -- improves efficiency if seeds are picked well

  if( !genFlat ){
     
    // the lines below should be tailored by the user for the particular desired
    // set of amplitudes -- doing so will improve efficiency.  Leaving as is
    // won't make MC incorrect, it just won't be as fast as it could be
    
    resProd.addResonance( 1.230, 0.400,  0.4 );
    resProd.addResonance( 1.318, 0.105,  0.3 );
    resProd.addResonance( 1.600, 0.200,  0.2 );
    resProd.addResonance( 1.670, 0.260,  0.4 );
  }

  vector< int > pTypes;
  pTypes.push_back( Gamma );
  pTypes.push_back( Neutron );
  pTypes.push_back( PiPlus );
  pTypes.push_back( PiMinus );
  pTypes.push_back( PiPlus );
  
  HDDMDataWriter* hddmOut = NULL;
  if( hddmname.size() != 0 ) hddmOut = new HDDMDataWriter( hddmname );
  ROOTDataWriter rootOut( outname );
  
  TFile* diagOut = new TFile( "gen_3pi_diagnostic.root", "recreate" );
  
  TH1F* mass = new TH1F( "M", "Resonance Mass", 180, lowMass, highMass );
  TH1F* massW = new TH1F( "M_W", "Weighted Resonance Mass", 180, lowMass, highMass );
  massW->Sumw2();
  TH1F* intenW = new TH1F( "intenW", "True PDF / Gen. PDF", 1000, 0, 100 );
  TH2F* intenWVsM = new TH2F( "intenWVsM", "Ratio vs. M", 100, lowMass, highMass, 1000, 0, 10 );
  
  TH2F* dalitz = new TH2F( "dalitz", "Dalitz Plot", 100, 0, 3.0, 100, 0, 3.0 );
  
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
      HepLorentzVector resonance( evt->particle( 2 ) + 
                                  evt->particle( 3 ) + 
                                  evt->particle( 4 ) );

      double genWeight = evt->weight();
      
      // cannot ask for the intensity if we haven't called process events above
      double weightedInten = ( genFlat ? 1 : ati.intensity( i ) );
      
      if( !diag ){
        
        // obtain this by looking at the maximum value of intensity * genWeight
        double rand = drand48() * maxInten;
        
        if( weightedInten > rand || genFlat ){
          
          mass->Fill( resonance.m() );
          massW->Fill( resonance.m(), genWeight );
          
          intenW->Fill( weightedInten );
          intenWVsM->Fill( resonance.m(), weightedInten );
          
          dalitz->Fill( ( evt->particle( 2 ) + evt->particle( 3 ) ).m2(),
                        ( evt->particle( 3 ) + evt->particle( 4 ) ).m2() );
           
          // we want to save events with weight 1
          evt->setWeight( 1.0 );
          
          if( hddmOut ) hddmOut->writeEvent( *evt, pTypes, true );
          rootOut.writeEvent( *evt );
          ++eventCounter;
        }
      }
      else{
        
        mass->Fill( resonance.m() );
        massW->Fill( resonance.m(), genWeight );

        dalitz->Fill( ( evt->particle( 2 ) + evt->particle( 3 ) ).m2(),
                      ( evt->particle( 3 ) + evt->particle( 4 ) ).m2() );

        intenW->Fill( weightedInten );
        intenWVsM->Fill( resonance.m(), weightedInten );
        
        ++eventCounter;
      }
      
      delete evt;
    }
    
    cout << eventCounter << " events were processed." << endl;
  }
  
  mass->Write();
  massW->Write();
  dalitz->Write();
  intenW->Write();
  intenWVsM->Write();
  diagOut->Close();
  
  if( hddmOut ) delete hddmOut;
  
	return 0;
}


