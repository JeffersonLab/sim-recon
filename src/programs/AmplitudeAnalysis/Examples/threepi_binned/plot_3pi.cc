#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <unistd.h>

#include "IUAmpTools/FitResults.h"

using namespace std;

int main( int argc, char* argv[] ){
  
  // these params should probably come in on the command line
  double lowMass = 0.7;
	double highMass = 2.0;
	enum{ kNumBins = 65 };
  string fitDir( "threepi_fit" );
  
  // set default parameters
  
  string outfileName("");
  
  // parse command line
  
  for (int i = 1; i < argc; i++){
    
    string arg(argv[i]);
    
    if (arg == "-o"){  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  outfileName = argv[++i]; }
    if (arg == "-h"){
      cout << endl << " Usage for: " << argv[0] << endl << endl;
      cout << "\t -o <file>\t Ouput text file" << endl;
      exit(1);}
  }
  
  if (outfileName.size() == 0){
    cout << "No output file specified" << endl;
    exit(1);
  }
  
	double step = ( highMass - lowMass ) / kNumBins;
	
  ofstream outfile;
  outfile.open( outfileName.c_str() );
  
  // descend into the directory that contains the bins
  chdir( fitDir.c_str() );
  
	for( int i = 0; i < kNumBins; ++i ){
		
    ostringstream dir;
    dir << "bin_" << i;
    chdir( dir.str().c_str() );

    ostringstream resultsFile;
    resultsFile << "bin_" << i << ".fit";
    
    FitResults results( resultsFile.str() );
    if( !results.valid() ){
      
      chdir( ".." );
      continue;
    }
    
  	// print out the bin center
		outfile << lowMass + step * i + step / 2. << "\t";
		  
    bool yPol = false;
    
    vector< string > rhoPiS;
		rhoPiS.push_back( "Pi+Pi-Pi+::xpol::J1_rhopi_S" );
    if( yPol ) rhoPiS.push_back( "Pi+Pi-Pi+::ypol::J1_rhopi_S" );
		pair< double, double > rhoPiSInt = results.intensity( rhoPiS );
		outfile << rhoPiSInt.first << "\t" << rhoPiSInt.second << "\t";
    
    vector< string > rhoPiD;
    rhoPiD.push_back( "Pi+Pi-Pi+::xpol::J2_rhopi_D" );
    if( yPol ) rhoPiD.push_back( "Pi+Pi-Pi+::ypol::J2_rhopi_D" );
    pair< double, double > rhoPiDInt = results.intensity( rhoPiD );
    outfile << rhoPiDInt.first << "\t" << rhoPiDInt.second << "\t";
    
    vector< string > rhoPiPX;
		rhoPiPX.push_back( "Pi+Pi-Pi+::xpol::J1_rhopi_P" );
    if( yPol ) rhoPiPX.push_back( "Pi+Pi-Pi+::ypol::J1_rhopi_P" );
    pair< double, double > rhoPiPXInt = results.intensity( rhoPiPX );
		outfile << rhoPiPXInt.first << "\t" << rhoPiPXInt.second << "\t";

    vector< string > f2PiS;
		f2PiS.push_back( "Pi+Pi-Pi+::xpol::J2_f2pi_S" );
    if( yPol ) f2PiS.push_back( "Pi+Pi-Pi+::ypol::J2_f2pi_S" );
		pair< double, double > f2PiSInt = results.intensity( f2PiS );
		outfile << f2PiSInt.first << "\t" << f2PiSInt.second << "\t";

    vector< string > rhoPiP;
		rhoPiP.push_back( "Pi+Pi-Pi+::xpol::J2_rhopi_P" );
    if( yPol ) rhoPiP.push_back( "Pi+Pi-Pi+::ypol::J2_rhopi_P" );
		pair< double, double > rhoPiPInt = results.intensity( rhoPiP );
		outfile << rhoPiPInt.first << "\t" << rhoPiPInt.second << "\t";
    
		vector< string > all;
		all.push_back( "Pi+Pi-Pi+::xpol::J1_rhopi_S" );
    all.push_back( "Pi+Pi-Pi+::xpol::J2_rhopi_D" );
    all.push_back( "Pi+Pi-Pi+::xpol::J1_rhopi_P" );
    all.push_back( "Pi+Pi-Pi+::xpol::J2_f2pi_S" );
    all.push_back( "Pi+Pi-Pi+::xpol::J2_rhopi_P" );
    if( yPol ) all.push_back( "Pi+Pi-Pi+::ypol::J1_rhopi_S" );
    if( yPol ) all.push_back( "Pi+Pi-Pi+::ypol::J2_rhopi_D" );
    if( yPol ) all.push_back( "Pi+Pi-Pi+::ypol::J1_rhopi_P" );
    if( yPol ) all.push_back( "Pi+Pi-Pi+::ypol::J2_f2pi_S" );
    if( yPol ) all.push_back( "Pi+Pi-Pi+::ypol::J2_rhopi_P" );
    pair< double, double > allInt = results.intensity( all );
		outfile << allInt.first << "\t" << allInt.second << "\t";
    
    pair< double, double > phaseDP =
    results.phaseDiff( "Pi+Pi-Pi+::xpol::J2_rhopi_D",
                       "Pi+Pi-Pi+::xpol::J1_rhopi_P" );
  
    outfile << phaseDP.first << "\t" << phaseDP.second << "\t";

    pair< double, double > phaseDS =
    results.phaseDiff( "Pi+Pi-Pi+::xpol::J2_rhopi_D",
                      "Pi+Pi-Pi+::xpol::J2_f2pi_S" );
    
    outfile << phaseDS.first << "\t" << phaseDS.second << "\t";

    outfile << endl;
    
    chdir( ".." );
	}
  
	return 0;
}
