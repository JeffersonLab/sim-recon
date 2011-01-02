
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cassert>
#include <utility>
#include <cstdlib>

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_DATAIO/ROOTDataReader.h"
#include "AMPTOOLS_DATAIO/ROOTDataWriter.h"

#include "CLHEP/Vector/LorentzVector.h"

#include "TRandom.h"

using namespace std;
using namespace CLHEP;

int main( int argc, char* argv[] ){
  
  string treeName( "kin" );
  
	if( argc < 3 ){
		
		cout << "Usage:  toy_detector <infile> <outfile>" << endl;
		return 1;
	}
	
	ROOTDataReader in(  argv[1], treeName );
  ROOTDataWriter out( argv[2] );
	
  Kinematics* event;
  while( ( event = in.getEvent() ) != NULL ){
    
    vector< HepLorentzVector > fs = event->particleList();

    HepLorentzVector x;
    // the first two entries in this list are the beam and the recoil
    // skip them in computing the mass
		for( vector< HepLorentzVector >::iterator particle = fs.begin() + 2;
        particle != fs.end(); ++particle ){
			
			x += *particle;
		}

    // an acceptance that is linearly rising with mass
    if( x.m() > drand48() * 3 )
      out.writeEvent( *event );
    
    delete event;
  }
  
	return 0;
}
