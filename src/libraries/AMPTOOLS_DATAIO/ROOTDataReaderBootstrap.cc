
#include <vector>
#include <cassert>
#include <iostream>
#include <string>
#include <cmath>

#include "TLorentzVector.h"

#include "AMPTOOLS_DATAIO/ROOTDataReaderBootstrap.h"
#include "IUAmpTools/Kinematics.h"

#include "TH1.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

ROOTDataReaderBootstrap::ROOTDataReaderBootstrap( const vector< string >& args ):
UserDataReader< ROOTDataReaderBootstrap >( args ),
m_eventCounter( 0 ),
m_useWeight( false )
{
  
  // arguments:
  // 0:  file name
  // 1:  random seed
  // 2:  tree name (optional; deafult: "kin")
  
  assert( args.size() == 3 || args.size() == 2 );
  
  TH1::AddDirectory( kFALSE );
  
  //this way of opening files works with URLs of the form
  // root://xrootdserver/path/to/myfile.root
  m_inFile = TFile::Open( args[0].c_str() );
  
  int seed = stoi( args[1] );
  m_randGenerator = new TRandom2( seed );
  
  cout << "******************** WARNING ***********************" << endl;
  cout << "*  You are using the boostrap data reader, which   *" << endl;
  cout << "*  should only be used for evaluating errors.      *" << endl;
  cout << "*  The results with different seeds will be random *" << endl;
  cout << "*  due to random oversampling of the input file.   *" << endl;
  cout << "****************************************************" << endl;
  cout << endl;
  cout << "   Random Seed:  " << seed << endl << endl;
  
  // default to tree name of "kin" if none is provided
  if( args.size() == 2 ){
    
    m_inTree = dynamic_cast<TTree*>( m_inFile->Get( "kin" ) );
  }
  else{
    
    m_inTree = dynamic_cast<TTree*>( m_inFile->Get( args[2].c_str() ) );
  }
  
  m_inTree->SetBranchAddress( "NumFinalState", &m_nPart );
  m_inTree->SetBranchAddress( "E_FinalState", m_e );
  m_inTree->SetBranchAddress( "Px_FinalState", m_px );
  m_inTree->SetBranchAddress( "Py_FinalState", m_py );
  m_inTree->SetBranchAddress( "Pz_FinalState", m_pz );
  m_inTree->SetBranchAddress( "E_Beam", &m_eBeam );
  m_inTree->SetBranchAddress( "Px_Beam", &m_pxBeam );
  m_inTree->SetBranchAddress( "Py_Beam", &m_pyBeam );
  m_inTree->SetBranchAddress( "Pz_Beam", &m_pzBeam );
  
  if(m_inTree->GetBranch("Weight") != NULL) {
    
    m_useWeight = true;
    m_inTree->SetBranchAddress( "Weight", &m_weight );
  }
  else{
    
    m_useWeight = false;
  }
}

ROOTDataReaderBootstrap::~ROOTDataReaderBootstrap()
{
  if( m_inFile != NULL ) m_inFile->Close();
  if( m_randGenerator ) delete m_randGenerator;
}

void
ROOTDataReaderBootstrap::resetSource()
{
  
  cout << "Resetting source " << m_inTree->GetName()
       << " in " << m_inFile->GetName() << endl;
  
  // this will cause the read to start back at event 0
  m_eventCounter = 0;
}

Kinematics*
ROOTDataReaderBootstrap::getEvent()
{
  if( m_eventCounter++ < numEvents() ){

    int thisEntry = (int)floor( m_randGenerator->Rndm()*numEvents() );
    
    m_inTree->GetEntry( thisEntry );
    assert( m_nPart < Kinematics::kMaxParticles );
    
    vector< TLorentzVector > particleList;
    
    particleList.
    push_back( TLorentzVector( m_pxBeam, m_pyBeam, m_pzBeam, m_eBeam ) );
    
    for( int i = 0; i < m_nPart; ++i ){
      
      particleList.push_back( TLorentzVector( m_px[i], m_py[i], m_pz[i], m_e[i] ) );
    }
    
    return new Kinematics( particleList, m_useWeight ? m_weight : 1.0 );
  }
  else{
    
    return NULL;
  }
}

unsigned int
ROOTDataReaderBootstrap::numEvents() const
{
  return static_cast< unsigned int >( m_inTree->GetEntries() );
}
