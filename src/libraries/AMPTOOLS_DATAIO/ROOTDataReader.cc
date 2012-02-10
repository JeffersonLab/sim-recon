
#include <vector>
#include <cassert>

#include "AMPTOOLS_DATAIO/ROOTDataReader.h"
#include "IUAmpTools/Kinematics.h"
#include "CLHEP/Vector/LorentzVector.h"

#include "TH1.h"
#include "TFile.h"
#include "TTree.h"

ROOTDataReader::ROOTDataReader( const string& inFileName,
				const string& inTreeName, bool useWeight) :
  m_eventCounter( 0 ),
  m_useWeight(useWeight)
{
  
	TH1::AddDirectory( kFALSE );

	//this way of opening files works with URLs of the form
	// root://xrootdserver/path/to/myfile.root
	m_inFile = TFile::Open( inFileName.c_str() ); 
	//m_inFile = new TFile( inFileName.c_str() );

	m_inTree = static_cast<TTree*>( m_inFile->Get( inTreeName.c_str() ) );
  
	m_inTree->SetBranchAddress( "nPart", &m_nPart );
	m_inTree->SetBranchAddress( "e", m_e );
	m_inTree->SetBranchAddress( "px", m_px );
	m_inTree->SetBranchAddress( "py", m_py );
	m_inTree->SetBranchAddress( "pz", m_pz );
	m_inTree->SetBranchAddress( "eBeam", &m_eBeam );
	m_inTree->SetBranchAddress( "pxBeam", &m_pxBeam );
	m_inTree->SetBranchAddress( "pyBeam", &m_pyBeam );
	m_inTree->SetBranchAddress( "pzBeam", &m_pzBeam );
	m_inTree->SetBranchAddress( "eRecoil", &m_eRecoil );
	m_inTree->SetBranchAddress( "pxRecoil", &m_pxRecoil );
	m_inTree->SetBranchAddress( "pyRecoil", &m_pyRecoil );
	m_inTree->SetBranchAddress( "pzRecoil", &m_pzRecoil );
	if(m_inTree->GetBranch("weight") != NULL)
	  m_inTree->SetBranchAddress( "weight", &m_weight );
	else
	  m_useWeight=false;
}

void
ROOTDataReader::resetSource()
{
	
	cout << "Resetting source " << m_inTree->GetName() 
	<< " in " << m_inFile->GetName() << endl;
  
	// this will cause the read to start back at event 0
	m_eventCounter = 0;
}

Kinematics*
ROOTDataReader::getEvent()
{
  if( m_eventCounter < static_cast< unsigned int >( m_inTree->GetEntries() ) ){
    //  if( m_eventCounter < 10 ){
    
    m_inTree->GetEntry( m_eventCounter++ );
    assert( m_nPart < Kinematics::kMaxParticles );
    
    vector< HepLorentzVector > particleList;
    
    particleList.
      push_back( HepLorentzVector( m_pxBeam, m_pyBeam, m_pzBeam, m_eBeam ) );
    particleList.
      push_back( HepLorentzVector( m_pxRecoil, m_pyRecoil, m_pzRecoil, m_eRecoil ) );
    
    for( int i = 0; i < m_nPart; ++i ){
      
      particleList.push_back( HepLorentzVector( m_px[i], m_py[i], m_pz[i], m_e[i] ) );
    }
    
    return new Kinematics( particleList, m_useWeight ? m_weight : 1.0 );
  }
  else{
    
    return NULL;
  }
}

unsigned int
ROOTDataReader::numEvents() const
{	
  return static_cast< unsigned int >( m_inTree->GetEntries() );
}
