
#include <cassert>
#include <cstdio>

#include "ASCIIDataWriter.h"
#include "CLHEP/Vector/LorentzVector.h"


ASCIIDataWriter::ASCIIDataWriter( const string& outFile )
{
  
  // Open output file
  fid=fopen((char *)(outFile.c_str()),"w");

  m_eventCounter = 0;
}

ASCIIDataWriter::~ASCIIDataWriter()
{
  fclose(fid);
}


/**
 * This function writes one event. It is presumed that the first 
 * two particles are the beam photon and the recoiling proton
 * in that order.
 * \param[in] kin - kinematics of the event
 * \param[in] types - geant particles types corresponding to the particles specified in kin 
 */

void
ASCIIDataWriter::writeEvent( const Kinematics& kin, vector<int> &types)
{
  vector< HepLorentzVector > particleList = kin.particleList();
	
  m_nPart = particleList.size() - 2;
  
  assert( particleList.size() <= Kinematics::kMaxParticles );
  
  
  // Start a new event
  fprintf(fid,"9000 %d %d\n",m_eventCounter+1,m_nPart+1);
  
  
  fprintf(fid,"1 %d %f\n",types[1],particleList[1].m());
  fprintf(fid,"   1 %f %f %f %f\n",
	  particleList[1].px(), particleList[1].py(),
	  particleList[1].pz(), particleList[1].t());
  
  for( int i = 0; i < m_nPart; i++ ){

    int charge=0;
    switch (types[i+2]){
    case 8: //pi plus
      charge=1; break;
    case 9: //pi minus
      charge=-1; break;
    case 11: //K plus
      charge=1; break;
    case 12: //K minus
      charge=-1; break;
    default:
      charge=0;
    }
    

    fprintf(fid,"%d %d %f\n",i+2,types[i+2],particleList[i+2].m());

    fprintf(fid,"   %d %f %f %f %f\n",charge,
	    particleList[i+2].px(),particleList[i+2].py(),
	    particleList[i+2].pz(),particleList[i+2].t());

  }
  
  m_eventCounter++;
}


