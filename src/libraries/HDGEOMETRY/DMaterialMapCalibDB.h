#ifndef _DMaterialMapCalibDB_
#define _DMaterialMapCalibDB_

#include "DMaterialMap.h"

#include <JANA/JApplication.h>
#include <JANA/JCalibration.h>
using namespace jana;

class DMaterialMapCalibDB:public DMaterialMap{
 public:
  DMaterialMapCalibDB(JApplication *japp);
  DMaterialMapCalibDB(JCalibration *jcalib);
  ~DMaterialMapCalibDB(){};
  
  unsigned int GetMaterialMap(unsigned int runnumber=1);

 protected:
  JCalibration *jcalib;
  
};

#endif // _DMaterialMapCalibDB_
