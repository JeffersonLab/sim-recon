#ifndef _DLorentzMapCalibDB_
#define _DLorentzMapCalibDB_

#include "DLorentzDeflections.h"

#include <JANA/JApplication.h>
#include <JANA/JCalibration.h>
using namespace jana;

class DLorentzMapCalibDB:public DLorentzDeflections{
 public:
  DLorentzMapCalibDB(JApplication *japp);
  DLorentzMapCalibDB(JCalibration *jcalib);
  ~DLorentzMapCalibDB(){};
  
  unsigned int GetLorentzDeflections(unsigned int runnumber=1);

 protected:
  JCalibration *jcalib;
  
};

#endif // _DLorentzMapCalibDB_
