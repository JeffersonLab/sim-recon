#if !defined(ROOTDATAREADERBOOTSTRAP)
#define ROOTDATAREADERBOOTSTRAP

#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/UserDataReader.h"

#include "TString.h"
#include "TRandom2.h"
#include "TFile.h"
#include "TTree.h"

#include <string>

using namespace std;

class ROOTDataReaderBootstrap : public UserDataReader< ROOTDataReaderBootstrap >
{
  
public:
  
  /**
   * Default constructor for ROOTDataReaderBootstrap
   */
  ROOTDataReaderBootstrap() : UserDataReader< ROOTDataReaderBootstrap >(), m_inFile( NULL ) { }
  
  ~ROOTDataReaderBootstrap();
  
  /**
   * Constructor for ROOTDataReaderBootstrap
   * \param[in] args vector of string arguments
   * arguments:
   *   0:  file name
   *   1:  random seeD
   *   2:  tree name (optional; deafult: "kin")
   */
  ROOTDataReaderBootstrap( const vector< string >& args );
  
  string name() const { return "ROOTDataReaderBootstrap"; }
  
  virtual Kinematics* getEvent();
  virtual void resetSource();
  
  /**
   * This function returns a true if the file was open
   * with weight-reading enabled and had this tree branch,
   * false, if these criteria are not met.
   */
  virtual bool hasWeight(){ return m_useWeight; };
  virtual unsigned int numEvents() const;
  
private:
  
  TFile* m_inFile;
  TTree* m_inTree;
  unsigned int m_eventCounter;
  bool m_useWeight;
  
  TRandom2* m_randGenerator;
  
  int m_nPart;
  float m_e[Kinematics::kMaxParticles];
  float m_px[Kinematics::kMaxParticles];
  float m_py[Kinematics::kMaxParticles];
  float m_pz[Kinematics::kMaxParticles];
  float m_eBeam;
  float m_pxBeam;
  float m_pyBeam;
  float m_pzBeam;
  float m_weight;
};

#endif
