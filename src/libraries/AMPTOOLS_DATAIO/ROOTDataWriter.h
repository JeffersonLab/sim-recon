#if !defined(ROOTDATAWRITER)
#define ROOTDATAWRITER

#include "IUAmpTools/Kinematics.h"

#include "TTree.h"
#include "TFile.h"

class ROOTDataWriter
{

public:
	
  /**
   * Constructor for ROOTDataWriter.
   *
   * \param[in] outFile name of output file
   * \param[in] outTreeName (optional) name of tree name in which to store events. Default: "kin"
   * \param[in] overwrite (optional) boolean parameter specifying whether to
   *               overwrite (default) or update the ROOT file.
   * \param[in] writeWeight (optional) enables writing of the event weight in the ROOT file
   */
  ROOTDataWriter( const string& outFile,
                  const string& outTreeName="kin",
                  bool overwrite=true, bool writeWeight=false )
  {
    IOinit(outFile, outTreeName, overwrite, writeWeight);
  };
 
  ~ROOTDataWriter();
  
  void writeEvent( const Kinematics& kin );
  
  int eventCounter() const { return m_eventCounter; }
  
private:
	
  
  void IOinit( const string& outFile,
	       const string& outTreeName,
	       bool overwrite, bool writeWeight);  
  
  TFile* m_outFile;
  TTree* m_outTree;
  int m_eventCounter;
  float m_weight;
  
  int m_nPart;
  float m_e[Kinematics::kMaxParticles];
  float m_px[Kinematics::kMaxParticles];
  float m_py[Kinematics::kMaxParticles];
  float m_pz[Kinematics::kMaxParticles];
  
  float m_eBeam;
  float m_pxBeam;
  float m_pyBeam;
  float m_pzBeam;
};

#endif
