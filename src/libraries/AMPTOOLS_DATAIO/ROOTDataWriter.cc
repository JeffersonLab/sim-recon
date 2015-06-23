#include <iostream>
#include <vector>
#include <cassert>

#include "AMPTOOLS_DATAIO/ROOTDataWriter.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"


void ROOTDataWriter::IOinit( const string& outFile,
                             const string& outTreeName,
                             bool overwrite, bool writeWeight)
{

  TH1::AddDirectory( kFALSE );
  string writeMode="recreate";
  if(!overwrite) writeMode="update";

  m_outFile = new TFile( outFile.c_str(), writeMode.c_str() );
  m_outTree = new TTree( outTreeName.c_str(), "Kinematics" );

  m_outTree->Branch( "NumFinalState", &m_nPart, "NumFinalState/I" );
  m_outTree->Branch( "E_FinalState", m_e, "E_FinalState[NumFinalState]/F" );
  m_outTree->Branch( "Px_FinalState", m_px, "Px_FinalState[NumFinalState]/F" );
  m_outTree->Branch( "Py_FinalState", m_py, "Py_FinalState[NumFinalState]/F" );
  m_outTree->Branch( "Pz_FinalState", m_pz, "Pz_FinalState[NumFinalState]/F" );
  m_outTree->Branch( "E_Beam", &m_eBeam, "E_Beam/F" );
  m_outTree->Branch( "Px_Beam", &m_pxBeam, "Px_Beam/F" );
  m_outTree->Branch( "Py_Beam", &m_pyBeam, "Py_Beam/F" );
  m_outTree->Branch( "Pz_Beam", &m_pzBeam, "Pz_Beam/F" );
  if(writeWeight)
    m_outTree->Branch( "Weight", &m_weight, "Weight/F" );  
  
  m_eventCounter = 0;
}

ROOTDataWriter::~ROOTDataWriter()
{
	m_outFile->cd();
	m_outTree->Write();
	m_outFile->Close();
}

void
ROOTDataWriter::writeEvent( const Kinematics& kin )
{
  vector< TLorentzVector > particleList = kin.particleList();
  
  m_nPart = particleList.size() - 1;
  
  assert( particleList.size() <= Kinematics::kMaxParticles );
  
  m_eBeam = particleList[0].E();
  m_pxBeam = particleList[0].Px();
  m_pyBeam = particleList[0].Py();
  m_pzBeam = particleList[0].Pz();
    
  for( int i = 0; i < m_nPart; ++i ){
    
    m_e[i] = particleList[i+1].E();
    m_px[i] = particleList[i+1].Px();
    m_py[i] = particleList[i+1].Py();
    m_pz[i] = particleList[i+1].Pz();
  }

  m_weight = kin.weight(); //will not get saved if branch not added in IOinit()
  
  m_outTree->Fill();
  
  m_eventCounter++;
  
}
