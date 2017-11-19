#if !defined(HDDMDATAWRITER)
#define HDDMDATAWRITER

#include <cstdlib>
#include <vector>
#include <cassert>
#include <fstream>

using namespace std;

#include "IUAmpTools/Kinematics.h"
#include "HDDM/hddm_s.hpp"

class HDDMDataWriter
{

public:
	
  HDDMDataWriter( const string& outFile, int runNumber=9000, int seed=0);
  ~HDDMDataWriter();
  
  void writeEvent( const Kinematics& kin, const vector<int>& ptype,
		   bool centeredVertex=false);
  void writeEvent( const Kinematics& kin, const vector<int>& ptype,
		   float vx, float vy, float vz_min, float vz_max);
  void writeEvent( const Kinematics& kin, const vector<int>& ptype,
		   float vx, float vy, float vz);
    
  int eventCounter() const { return m_eventCounter; }
  bool FileOpen() { return m_OutputStream; }
  
private:

  std::ofstream *m_OutputFile;        // output hddm file ofstream
  hddm_s::ostream *m_OutputStream;    // provides hddm layer on top of ofstream
  int m_eventCounter, m_runNumber;

};

#endif
