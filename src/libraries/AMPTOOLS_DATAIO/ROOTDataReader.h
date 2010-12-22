#if !defined(ROOTDATAREADER)
#define ROOTDATAREADER

#include "IUAmpTools/DataReader.h"

#include "TString.h"
#include "TFile.h"
#include "TTree.h"

#include <string>

using namespace std;

class ROOTDataReader : public DataReader
{
	
public:
	
	enum { kMaxFSParticles = 10 };
	
	ROOTDataReader( const string& inFileName, 
                  const string& inTreeName = "kin" );
	
	Kinematics* getEvent();
	void resetSource();
	
	unsigned int numEvents() const;
	int eventCounter() const { return m_eventCounter; }
	
private:
	
	TFile* m_inFile;
	TTree* m_inTree;
    int m_eventCounter;
	
	int m_nPart;
	float m_e[kMaxFSParticles];
	float m_px[kMaxFSParticles];
	float m_py[kMaxFSParticles];
	float m_pz[kMaxFSParticles];
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
