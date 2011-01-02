#if !defined(ROOTDATAREADER)
#define ROOTDATAREADER

#include "IUAmpTools/DataReader.h"
#include "IUAmpTools/Kinematics.h"

#include "TString.h"
#include "TFile.h"
#include "TTree.h"

#include <string>

using namespace std;

class ROOTDataReader : public DataReader
{
	
public:
  
	ROOTDataReader( const string& inFileName, 
                 const string& inTreeName );
	
	virtual Kinematics* getEvent();
	virtual void resetSource();
	
	virtual unsigned int numEvents() const;
	
private:
	
	TFile* m_inFile;
	TTree* m_inTree;
  unsigned int m_eventCounter;
	
	int m_nPart;
	float m_e[Kinematics::kMaxParticles];
	float m_px[Kinematics::kMaxParticles];
	float m_py[Kinematics::kMaxParticles];
	float m_pz[Kinematics::kMaxParticles];
	float m_eBeam;
	float m_pxBeam;
	float m_pyBeam;
	float m_pzBeam;
	float m_eRecoil;
	float m_pxRecoil;
	float m_pyRecoil;
	float m_pzRecoil;
};

#endif
