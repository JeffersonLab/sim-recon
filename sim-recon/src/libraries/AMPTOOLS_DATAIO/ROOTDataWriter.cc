#include <iostream>
#include <vector>
#include <cassert>

#include "AMPTOOLS_DATAIO/ROOTDataWriter.h"
#include "CLHEP/Vector/LorentzVector.h"

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

  m_outTree->Branch( "nPart", &m_nPart, "nPart/I" );
  m_outTree->Branch( "e", m_e, "e[nPart]/F" );
  m_outTree->Branch( "px", m_px, "px[nPart]/F" );
  m_outTree->Branch( "py", m_py, "py[nPart]/F" );
  m_outTree->Branch( "pz", m_pz, "pz[nPart]/F" );
  m_outTree->Branch( "eBeam", &m_eBeam, "eBeam/F" );
  m_outTree->Branch( "pxBeam", &m_pxBeam, "pxBeam/F" );
  m_outTree->Branch( "pyBeam", &m_pyBeam, "pyBeam/F" );
  m_outTree->Branch( "pzBeam", &m_pzBeam, "pzBeam/F" );
  m_outTree->Branch( "eRecoil", &m_eRecoil, "eRecoil/F" );
  m_outTree->Branch( "pxRecoil", &m_pxRecoil, "pxRecoil/F" );
  m_outTree->Branch( "pyRecoil", &m_pyRecoil, "pyRecoil/F" );
  m_outTree->Branch( "pzRecoil", &m_pzRecoil, "pzRecoil/F" );
  if(writeWeight) 
    m_outTree->Branch( "weight", &m_weight, "weight/F" );  
  
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
  vector< HepLorentzVector > particleList = kin.particleList();
  
  m_nPart = particleList.size() - 2;
  
  assert( particleList.size() <= Kinematics::kMaxParticles );
  
  m_eBeam = particleList[0].t();
  m_pxBeam = particleList[0].px();
  m_pyBeam = particleList[0].py();
  m_pzBeam = particleList[0].pz();
  
  m_eRecoil = particleList[1].t();
  m_pxRecoil = particleList[1].px();
  m_pyRecoil = particleList[1].py();
  m_pzRecoil = particleList[1].pz();
  
  for( int i = 0; i < m_nPart; ++i ){
    
    m_e[i] = particleList[i+2].t();
    m_px[i] = particleList[i+2].x();
    m_py[i] = particleList[i+2].y();
    m_pz[i] = particleList[i+2].z();
  }

  m_weight = kin.weight();//will not get saved if branch not added in IOinit()
  
  m_outTree->Fill();
  
  m_eventCounter++;
  
}
