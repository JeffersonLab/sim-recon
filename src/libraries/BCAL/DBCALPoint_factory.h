#ifndef _DBCALPoint_factory_
#define _DBCALPoint_factory_

#include <vector>
#include <map>
using namespace std;

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>

using namespace jana;

#include "BCAL/DBCALPoint.h"
#include "BCAL/DBCALUnifiedHit.h"

#include <TTree.h>

typedef vector< vector<double> >  attenuation_parms_t;
typedef vector< double >          effective_vel_t;

class DBCALHit;

class DBCALPoint_factory : public JFactory<DBCALPoint> {

 public:
  DBCALPoint_factory() {
    PRINTCALIBRATION = false;
    if(gPARMS){
      gPARMS->SetDefaultParameter("BCALPOINT:PRINTCALIBRATION", PRINTCALIBRATION, "Print the calibration parameters.");
    }
  }
  ~DBCALPoint_factory() {}

 private:
  class cellHits{
   public:
    vector<const DBCALUnifiedHit*> uphits;
    vector<const DBCALUnifiedHit*> dnhits;
  };

  double m_z_target_center;
  attenuation_parms_t attenuation_parameters;
  effective_vel_t effective_velocities;
 
  static const int BCAL_NUM_MODULES  = 48;
  static const int BCAL_NUM_LAYERS   =  4;
  static const int BCAL_NUM_SECTORS  =  4;

  bool PRINTCALIBRATION;

  jerror_t brun(JEventLoop *loop, int runnumber);
  jerror_t evnt(JEventLoop *loop, int eventnumber);

  const int GetCalibIndex( int module, int layer, int sector ) const {
	  return BCAL_NUM_LAYERS*BCAL_NUM_SECTORS*(module-1) + BCAL_NUM_SECTORS*(layer-1) + (sector-1);
  }

  bool GetAttenuationParameters(int id, double &attenuation_length,
				double &attenuation_L1, double &attenuation_L2);
  double GetEffectiveVelocity(int id);
};

#endif //_DBCALPoint_factory_
