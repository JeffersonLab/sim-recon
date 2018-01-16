#if !defined(ROOTDATAREADERWITHTCUT)
#define ROOTDATAREADERWITHTCUT

#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/UserDataReader.h"

#include "TString.h"
#include "TFile.h"
#include "TTree.h"

#include <string>

using namespace std;

class ROOTDataReaderWithTCut : public UserDataReader< ROOTDataReaderWithTCut >
{
	
public:
  
  /**
   * Default constructor for ROOTDataReaderWithTCut
   */
  ROOTDataReaderWithTCut() : UserDataReader< ROOTDataReaderWithTCut >(), m_inFile( NULL ) { }
  
  ~ROOTDataReaderWithTCut();
  
  /**
   * Constructor for ROOTDataReaderWithTCut
   * \param[in] args vector of string arguments
   */
  ROOTDataReaderWithTCut( const vector< string >& args );
  
  string name() const { return "ROOTDataReaderWithTCut"; }
  
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
  unsigned int m_eventCounter,m_numEvents;
  bool m_useWeight, m_RangeSpecified;
  double m_tMin,m_tMax;
  
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
