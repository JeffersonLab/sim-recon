
#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include <cassert>
#include <cstdlib>

#include "AMPTOOLS_DATAIO/ROOTDataWriter.h"

#include "AMPTOOLS_AMPS/ThreePiAngles.h"
#include "AMPTOOLS_AMPS/BreitWigner.h"

#include "AMPTOOLS_MCGEN/ProductionMechanism.h"
#include "AMPTOOLS_MCGEN/GammaPToXYZP.h"

#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "CLHEP/Vector/LorentzVector.h"

#include "TH1F.h"
#include "TH2F.h"

using std::complex;
using namespace std;
using namespace CLHEP;

int main( int argc, char* argv[] ){
  
  string  configfile("");
  string  outname("");
  
  bool diag = false;
  bool genFlat = false;
  
  // default upper and lower bounds 
  double lowMass = 0.7;
  double highMass = 2.0;
  
  int nEvents = 1E5;
  int batchSize = 1E5;
    
	//parse command line:
  for (int i = 1; i < argc; i++){
    
    string arg(argv[i]);
    
    if (arg == "-c"){  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  configfile = argv[++i]; }
    if (arg == "-o"){  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  outname = argv[++i]; }
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
      cout << "\t -c <file>\t Config file" << endl;
      cout << "\t -o <name>\t Output name" << endl;
      cout << "\t -l <value>\t Low edge of mass range (GeV) [optional]" << endl;
      cout << "\t -u <value>\t Upper edge of mass range (GeV) [optional]" << endl;
      cout << "\t -n <value>\t Minimum number of events to generate [optional]" << endl;
      cout << "\t -f \t\t Generate flat in M(X) (no physics) [optional]" << endl;
      cout << "\t -d \t\t Plot only diagnostic histograms [optional]" << endl << endl;
      exit(1);
    }
  }
  
  if( configfile.size() == 0 || outname.size() == 0 ){
    cout << "No config file or output specificed:  run gen_3pi -h for help" << endl;
    exit(1);
  }
  
  // open config file
  ConfigFileParser parser( configfile );
  ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
  vector< ReactionInfo* > reactions = cfgInfo->reactionList();
  assert( reactions.size() == 1 );
  
  ReactionInfo* info = reactions[0];
  
  // setup amplitude manager
  AmplitudeManager ampManager( info->particleList(), info->reactionName() );	
  ampManager.registerAmplitudeFactor( ThreePiAngles() );
  ampManager.registerAmplitudeFactor( BreitWigner() );
  ampManager.setupFromConfigurationInfo( cfgInfo );  
  
  
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
  
  // open output file
	ROOTDataWriter rootOut( outname );
  
  TH1F* mass = new TH1F( "M", "Resonance Mass", 180, 0.7, 2.5 );
  TH1F* massW = new TH1F( "M_W", "Weighted Resonance Mass", 180, 0.7, 2.5 );
  massW->Sumw2();
  TH1F* intenW = new TH1F( "intenW", "True PDF / Gen. PDF", 1000, 0, 100 );
  TH2F* intenWVsM = new TH2F( "intenWVsM", "Ratio vs. M", 100, .7, 2.5, 1000, 0, 10 );
  
  TH2F* dalitz = new TH2F( "dalitz", "Dalitz Plot", 100, 0, 3.0, 100, 0, 3.0 );
  
  int eventCounter = 0;
  while( eventCounter < nEvents ){
    
    if( batchSize < 1E4 ){
      
      cout << "WARNING:  small batches could have batch-to-batch variations\n"
      << "          due to different maximum intensities!" << endl;
    }
    
    cout << "Generating four-vectors..." << endl;
    AmpVecs* aVecs = resProd.generateMany( batchSize );
    cout << "Calculating amplitudes..." << endl;
    aVecs->allocateAmps( ampManager, true );
    
    // include factor of 1.5 to be safe in case we miss peak
    double maxInten = 1.5 * ampManager.calcIntensities( *aVecs );
    
    cout << "Processing events.." << endl;
    
    for( int i = 0; i < batchSize; ++i ){
  		
      Kinematics* evt = aVecs->getEvent( i );
      HepLorentzVector resonance( evt->particle( 2 ) + 
                                  evt->particle( 3 ) + 
                                  evt->particle( 4 ) );

      double genWeight = evt->weight();
      
      double weightedInten = aVecs->m_pdIntensity[i];
      
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
    delete aVecs;
  }
  
  mass->Write();
  massW->Write();
  dalitz->Write();
  intenW->Write();
  intenWVsM->Write();
  
	return 0;
}


