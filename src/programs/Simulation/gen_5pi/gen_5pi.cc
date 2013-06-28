
#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include <cassert>
#include <cstdlib>

using std::complex;
using namespace std;

#include "particleType.h"
#include "AMPTOOLS_DATAIO/ROOTDataReader.h"
#include "AMPTOOLS_DATAIO/ROOTDataWriter.h"
#include "AMPTOOLS_DATAIO/ASCIIDataWriter.h"
#include "AMPTOOLS_DATAIO/HDDMDataWriter.h"

#include "AMPTOOLS_AMPS/b1piAngAmp.h"
#include "AMPTOOLS_AMPS/BreitWigner.h"

#include "AMPTOOLS_MCGEN/ProductionMechanism.h"
#include "AMPTOOLS_MCGEN/GammaPToNPartP.h"
//#include "GammaPTob1piP.h"

#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/LorentzRotation.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

#include "b1piAmpCheck.h"

using namespace CLHEP;

void Usage(char* progName )
{
  cout << endl << " Usage for: " << progName << endl << endl;
  cout << "    -c <file>\t Config file [required]" << endl;
  cout << "    -o <name>\t Output name base (ascii, hddm and root files generated)" << endl;
  cout << "    -a <file>\t File to record events before accept/reject" << endl;
  cout << "    -e <file>\t Input file with MC sample to use instead of generating" << endl;
  cout << "    -l <value>\t Low edge of mass range (GeV)" << endl;
  cout << "    -u <value>\t Upper edge of mass range (GeV)" << endl;
  cout << "    -n <value>\t Minimum number of events to generate (allows spill-over)" << endl;
  cout << "    -N <value>\t Exact number of events to generate" << endl;
  cout << "    -f \t\t Generate flat in M(X) (no physics)" << endl;
  cout << "    -d \t\t Compute diagnostic histograms" << endl ;
  cout << "    -s <value>\t Specify random number generator seed" << endl;
  cout << "    -b <value>\t Batch size for intensities in accept/reject alg. (def. 200k)" << endl; 
  cout << "    -i <value>\t Specify maximum intensity (accept/reject range)" << endl << endl;
  
}


bool fileGood(string fname)
{
  ifstream file(fname.c_str());
  return file.good();
}


int main( int argc, char* argv[] ){
  
  string  configfile("");
  string  outname("gen_5pi"), allGenFName, inMCFName;
  b1piAmpCheck AmpCheck;
  bool diag = false, genFlat = false, StrictEvtLimit=false;
  bool saveAll=false, readInEvents=false;
  
  // default upper and lower bounds 
  double lowMass = 0.7, highMass = 3.0, Mpipm,Mpi0;
  Mpipm=ParticleMass(PiPlus);
  Mpi0=ParticleMass(Pi0);

  //Exprected particle list: 
  // pi- b1(pi+ omega(pi0 "rho"(pi- pi+)))
  //  2      3         4         5   6
  int par_types_list[]={1,14,9,8,7,9,8};
  vector<int> part_types(par_types_list,par_types_list+7);
  float part_masses_list[]={Mpipm, Mpipm, Mpi0, Mpipm, Mpipm};
  vector<double> part_masses(part_masses_list,part_masses_list+5);

  
  double CustMaxInten=-1, maxInten=0, runs_maxInten=0.0;
  long int Seed=0;

  int nEvents = 30, batchSize = 200000;

  //Parse command line: -----------------------------------------------
  for (int i = 1; i < argc; i++){
    
    string arg(argv[i]);
    
    if (arg == "-c"){  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  configfile = argv[++i]; }
    if (arg == "-o"){  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  outname = argv[++i]; }
    if (arg == "-a"){  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  allGenFName = argv[++i]; saveAll=true;}
    if (arg == "-e"){  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  inMCFName = argv[++i]; readInEvents=true;}
    if (arg == "-l"){  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  lowMass = atof( argv[++i] ); }
    if (arg == "-u"){  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  highMass = atof( argv[++i] ); }
    if (arg == "-n"){  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  nEvents = atoi( argv[++i] ); }
    if (arg == "-N"){  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else { nEvents = atoi( argv[++i] ); StrictEvtLimit=true;} }
    if (arg == "-b"){  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  batchSize = atoi( argv[++i] ); }
    if (arg == "-s"){  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else Seed = atoi( argv[++i] ); }
    if (arg == "-i"){  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else CustMaxInten = runs_maxInten = atof( argv[++i] ); }
    
    if (arg == "-d"){  
      diag = true; }
    if (arg == "-f"){  
      genFlat = true; }
    if (arg == "-h"){
      Usage(argv[0]);
      exit(1);
    }
  }
  
  if( configfile.size() == 0 ){
    cerr << "No config file specified." << endl;
    Usage(argv[0]);
    exit(1);
  }
  if(!fileGood(configfile)){
    cerr << "Invalid config file!" << endl;
    exit(1);
  }    
     

  // END OF ARGUMENT PARSING/CHECKING /////////////////////////////////////////

  
  srand48(Seed);

  // prepare output pipes -----------------------------
  string rootfname=outname + ".root";
  string asciifname=outname + ".ascii";
  string hddmfname=outname + ".hddm";
  
  // open output file
  ROOTDataWriter *rootAllOut=NULL;
  if(saveAll) rootAllOut= new ROOTDataWriter( allGenFName, "kin", true );
  HDDMDataWriter hddmOut( hddmfname );
  ROOTDataWriter rootOut( rootfname);
  ASCIIDataWriter asciiOut( asciifname ); 
  //----------------------------------------------------

  vector< string > readerArgs;
  readerArgs.push_back( inMCFName );
  readerArgs.push_back( "kin" );
  
  if(readInEvents){
    cout << "Performing accept/reject on pre-generated sample from: " << inMCFName << endl;
    ROOTDataReader rootIn( readerArgs );
    Kinematics *evt;

    if(CustMaxInten>0){
      cout << "  using the specified intensity ceiling: " << CustMaxInten << endl;
      maxInten=CustMaxInten;
    }else{
      cout << "Looking for peak intensity...\n";
      int i=0;
      while((evt=rootIn.getEvent())!=NULL){
	if(evt->weight() > maxInten) maxInten = evt->weight();
	delete evt;
	if(++i % 50000==0) cout << i << " events searched\n";
      }

      cout << "Maximum intensity found: " << maxInten << endl;	  
      rootIn.resetSource();
    }

    for(int i=0 ; (evt=rootIn.getEvent())!=NULL && 
	  ((StrictEvtLimit && i<nEvents) || !StrictEvtLimit) ; ++i){
      if(evt->weight() > drand48() * maxInten * 1.5 ) {
	evt->setWeight( 1.0 );
	rootOut.writeEvent( *evt );
	asciiOut.writeEvent( *evt, part_types );
	hddmOut.writeEvent( *evt, part_types );
      }
      delete evt;
      if(i%9999==0) cout << i+1 << " events written\n";
    }
    
    return 0 ;
  }
  



  // prepare dyagnostic histograms ---------------------
  TH1F* mass = new TH1F( "M", "Resonance Mass", 180, lowMass, highMass );
  TH1F* massW = new TH1F( "M_W", "Weighted Resonance Mass", 180, lowMass, highMass );
  massW->Sumw2();
  TH1F* intenW = new TH1F( "intenW", "True PDF / Gen. PDF", 1000, 0, 100 );
  TH2F* intenWVsM = new TH2F( "intenWVsM", "Ratio vs. M", 100, lowMass, highMass, 1000, 0, 10 );
  
  TH2F* dalitz = new TH2F( "dalitz", "Dalitz Plot", 100, 0, 3.0, 100, 0, 3.0 );
  
  TH1F* prod_ang = new TH1F( "alpha", "Production angle #alpha", 100, -PI, PI);
  
  TH1F* M_rho = new TH1F( "M_rho", "#rho Mass", 200, 0.27, .9 );
  TH1F* M_omega = new TH1F( "M_omega", "#omega Mass", 200, 0.6, 1.1 );
  TH1F* M_b1 = new TH1F( "M_b1", "Isobar Mass", 200, 0.5, 3.0 );
  
  TH2F *XAng = new TH2F("XAng","Angular distribution of X decay",
			50,-M_PI,M_PI,50,-1,1);
  TH2F *b1Ang = new TH2F("b1Ang","Angular distribution of b_{1} decay",
			 50,-M_PI,M_PI,50,-1,1);
  TH2F *OmegaAng = new TH2F("OmegaAng","Angular distribution of #omega decay",
			    50,-M_PI,M_PI,50,-1,1);
  TH2F *RhoAng = new TH2F("RhoAng","Angular distribution of #rho decay",
			  50,-M_PI,M_PI,50,-1,1);
  // -----------------------------------------------------

  

  //FILE *Ifid;

  
  // open config file and be sure only one reaction is specified
  ConfigFileParser parser( configfile );
  ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
  assert( cfgInfo->reactionList().size() == 1 );
  ReactionInfo* reaction = cfgInfo->reactionList()[0];
  
  AmpToolsInterface::registerAmplitude( b1piAngAmp() );
  AmpToolsInterface::registerAmplitude( BreitWigner() );
  AmpToolsInterface ati( cfgInfo );
  
  ProductionMechanism::Type type =
    ( genFlat ? ProductionMechanism::kFlat : ProductionMechanism::kResonant );
  
  //generate over a range mass -- the daughters are pi-,pi+,omega
  GammaPToNPartP resProd( lowMass, highMass, part_masses,type );
  //GammaPTob1piP resProd( lowMass, highMass, type );
  
  // seed the distribution with a sum of noninterfering Breit-Wigners
  // we can easily compute the PDF for this and divide by that when
  // doing accept/reject -- improves efficiency if seeds are picked well
  if( !genFlat ){
    // the lines below should be tailored by the user for the particular desired
    // set of amplitudes -- doing so will improve efficiency.  Leaving as is
    // won't make MC incorrect, it just won't be as fast as it could be
    resProd.addResonance( 1.89, 0.16,  1.0 );
    resProd.addResonance( 2.0, 0.25,  0.75 );
  }

  
  //Ifid=fopen("Idump.dat","w");
  
  int eventCounter = 0;
  while( eventCounter < nEvents ){
    
    if( batchSize < 1E4 ){  
      cout << "WARNING:  small batches could have batch-to-batch variations\n"
	   << "          due to different maximum intensities!" << endl;
    }
    
    cout << "Number to generate: " << batchSize << endl;
    cout << "Generating four-vectors..." << endl;
    
    ati.clearEvents();
    for( int i = 0; i < batchSize; ++i ){
      
      Kinematics* kin = resProd.generate();
      ati.loadEvent( kin, i, batchSize );
      delete kin;
    }

    cout << "Calculating intensities..." << endl;
    
    maxInten=0;
    if(!genFlat) // Signal intensity calculation
      maxInten = ati.processEvents( reaction->reactionName() );
    else
      for(int i = 0; i < batchSize; ++i ){ // Get the max. inten. for PS gen
        Kinematics* kin = ati.kinematics(i);
        if(kin->weight() > maxInten)
          maxInten = kin->weight();
        delete kin;
      }
    
    cout << "Beginning accept/reject..." << endl;
    
    printf("MAXINTEN: %25.20f\n",maxInten);
    if(runs_maxInten < maxInten) runs_maxInten=maxInten;
    
    //override the max intensity found with that passed in through cmd line args
    if( CustMaxInten > 0 ){
      if(maxInten>CustMaxInten){
	printf("WARNING: Event found with intensity greater than custom-specified maximum\n");
	CustMaxInten=maxInten;
      }else maxInten=CustMaxInten;
    }else maxInten*=1.5;
    
    
    cout << "Processing events..." << endl;
    
    
    //double IbatchSum=0;
    for( int i = 0; i < batchSize ; ++i ){
      
      Kinematics* evt = ati.kinematics( i );
      

      double genWeight = evt->weight();
      double weightedInten = ati.intensity( i );
      


      // obtain this by looking at the maximum value of intensity * genWeight
      if((!genFlat && weightedInten > drand48() * maxInten) ||
	 (genFlat && genWeight > drand48() * maxInten) ){
	
	double histWeight = 1.0;//genFlat ? genWeight : 1.0; 
	
	//Fill some useful histograms
	if(diag){
	  AmpCheck.SetEvent(*evt);
	  mass->Fill( AmpCheck.GetMX(), histWeight );
	  massW->Fill( AmpCheck.GetMX(), genWeight );
	  
	  intenW->Fill( weightedInten, histWeight );
	  intenWVsM->Fill( AmpCheck.GetMX(), weightedInten );
	  
	  dalitz->Fill( (AmpCheck.GetXsPi() + AmpCheck.Getb1sPi()).M2(),
			(AmpCheck.Getb1sPi() + AmpCheck.GetOmega()).M2(),histWeight);
	  M_rho->Fill(AmpCheck.GetMrho(), histWeight);
	  M_omega->Fill(AmpCheck.GetMomega(), histWeight);
	  M_b1->Fill(AmpCheck.GetMb1(), histWeight);
	  
	  // orientation of production plane in lab
	  prod_ang->Fill(AmpCheck.GetAlpha(), histWeight);
	  XAng->Fill(AmpCheck.GetXPhi(), AmpCheck.GetXCosTheta(), histWeight);
	  b1Ang->Fill(AmpCheck.Getb1Phi(), AmpCheck.Getb1CosTheta(), histWeight);
	  OmegaAng->Fill(AmpCheck.GetOmegaPhi(),AmpCheck.GetOmegaCosTheta(), histWeight);
	  RhoAng->Fill(AmpCheck.GetRhoPhi(), AmpCheck.GetRhoCosTheta(), histWeight);
	}
	// we want to save events with weight 1
	evt->setWeight( 1.0 );
	
	rootOut.writeEvent( *evt ); 	  
	asciiOut.writeEvent( *evt, part_types );
	hddmOut.writeEvent( *evt, part_types );
	++eventCounter;
	if(StrictEvtLimit && eventCounter>=nEvents) break;
      }
      
      if(saveAll) {
	evt->setWeight( weightedInten );
	rootAllOut->writeEvent( *evt );
      }      
      delete evt;
    }
    //printf("BATCH AVERAGE:  %f\n",IbatchSum/batchSize);
    
    cout << eventCounter << " events were processed." << endl;
  }

  if(!genFlat) printf("RUN_MAX_INTEN: %25.20f\n",runs_maxInten);
  
  mass->Write();
  massW->Write();
  dalitz->Write();
  intenW->Write();
  intenWVsM->Write();
  prod_ang->Write();
  M_b1->Write();
  M_omega->Write();
  M_rho->Write();
  XAng->Write();
  b1Ang->Write();
  OmegaAng->Write();
  RhoAng->Write();
  if(saveAll) delete rootAllOut;
  //fclose(Ifid);

  return 0;
}


