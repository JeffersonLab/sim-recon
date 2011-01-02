#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cassert>
#include <cstdlib>

#include "AMPTOOLS_DATAIO/ROOTDataReader.h"
#include "AMPTOOLS_DATAIO/ROOTDataWriter.h"

#include "CLHEP/Vector/LorentzVector.h"

#include "TH1F.h"
using namespace std;
using namespace CLHEP;

int main( int argc, char* argv[] ){
  
  int maxEvents = 1E10;
  
  string treeName( "kin" );
  
	if( argc < 6 ){
		
		cout << "Usage:  split_mass <infile> <outputBase> <lowMass> <highMass> <nBins> [maxEvents]" << endl;
		return 1;
	}
	
	ROOTDataReader in( argv[1], treeName );
	string outBase( argv[2] );
	
	double lowMass = atof( argv[3] );
	double highMass = atof( argv[4] );
	int numBins = atoi( argv[5] );
	
  if( argc == 7 ) maxEvents = atoi( argv[6] );
  
  enum { kMaxBins = 1000 };
  assert( numBins < kMaxBins );
  
	double step = ( highMass - lowMass ) / numBins;
	
	ROOTDataWriter* outFile[kMaxBins];
	
	for( int i = 0; i < numBins; ++i ){
		
		ostringstream outName;
		outName << outBase << "_" << i << ".root";
		outFile[i] = new ROOTDataWriter( outName.str() );
	}
	
  int eventCount = 0;
  
	Kinematics* event;
	while( ( event = in.getEvent() ) != NULL && eventCount++ < maxEvents ){
		
		vector< HepLorentzVector > fs = event->particleList();
		
		HepLorentzVector x;
    // the first two entries in this list are the beam and the recoil
    // skip them in computing the mass
		for( vector< HepLorentzVector >::iterator particle = fs.begin() + 2;
        particle != fs.end(); ++particle ){
			
			x += *particle;
		}
		
		int bin = static_cast< int >( floor( ( x.m() - lowMass ) / step ) );
		if( ( bin < numBins ) && ( bin >= 0 ) ){
			
			outFile[bin]->writeEvent( *event );
		}
	}
	
	for( int i = 0; i < numBins; ++i ){
		
		delete outFile[i];
	}
	
	return 0;
}
