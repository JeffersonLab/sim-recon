#if !defined(HDDMDATAWRITER)
#define HDDMDATAWRITER

#include <cstdlib>
#include <vector>
#include <cassert>

using namespace std;

#include "AMPTOOLS_DATAIO/HDDMDataWriter.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "IUAmpTools/Kinematics.h"
#include "HDDM/hddm_s.h"

using namespace CLHEP;

class HDDMDataWriter
{

public:
	
  HDDMDataWriter( const string& outFile, int runNumber=9000);
  ~HDDMDataWriter();
  
  void writeEvent( const Kinematics& kin, vector<int> ptype, 
		   bool centeredVertex=false);
  void writeEvent( const Kinematics& kin, vector<int> ptype,
		   float vx, float vy, float vz_min, float vz_max);
  void writeEvent( const Kinematics& kin, vector<int> ptype,
		   float vx, float vy, float vz);
    
  int eventCounter() const { return m_eventCounter; }
  bool FileOpen(){return m_OutputStream;};
  
private:

  s_iostream_t *m_OutputStream;
  int m_eventCounter, m_runNumber;

};

#endif
