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

typedef map<int, vector<double> >  attenuation_parms_t;

class DBCALHit;

class DBCALPoint_factory : public JFactory<DBCALPoint> {

 public:
  DBCALPoint_factory() {}
  ~DBCALPoint_factory() {}

 private:
  class cellHits{
   public:
    vector<const DBCALUnifiedHit*> uphits;
    vector<const DBCALUnifiedHit*> dnhits;
  };

  double m_z_target_center;

  attenuation_parms_t attenuation_parameters;
 
  jerror_t brun(JEventLoop *loop, int runnumber);
  jerror_t evnt(JEventLoop *loop, int eventnumber);

  bool GetAttenuationParameters(int id, double &attenuation_length,
				double &attenuation_L1, double &attenuation_L2);
};

#endif //_DBCALPoint_factory_
