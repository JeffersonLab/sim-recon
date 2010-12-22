#if !defined(ROOTDATAWRITER)
#define ROOTDATAWRITER

#include "IUAmpTools/Kinematics.h"

#include "TTree.h"
#include "TFile.h"

class ROOTDataWriter
{

public:
	
	enum { kMaxFSParticles = 10 };
	
	ROOTDataWriter( const string& outFile );
	~ROOTDataWriter();
	
	void writeEvent( const Kinematics& kin );
	
	int eventCounter() const { return m_eventCounter; }
	
private:
		
	TFile* m_outFile;
	TTree* m_outTree;
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
