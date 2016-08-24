
#ifndef _DPSGeometry_
#define _DPSGeometry_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

#include <string>

class DPSGeometry : public JObject {
  
 public:
  
  JOBJECT_PUBLIC(DPSGeometry);
  
  DPSGeometry(JEventLoop *loop);
  ~DPSGeometry();

  enum Arm { kNorth, kSouth };
  
  static const int NUM_ARMS = 2;
  
  // number of channels in coarse and fine detectors in each arm
  static const int NUM_COARSE_COLUMNS = 8;
  static const int NUM_FINE_COLUMNS = 145;
  
  // columns are numbered 1..NUM_FINE_COLUMNS; arm is 0:North or 1:South
  double getElow(int arm,int column) const;
  double getEhigh(int arm,int column) const;

 private:
  double m_energy_low[NUM_ARMS][NUM_FINE_COLUMNS];
  double m_energy_high[NUM_ARMS][NUM_FINE_COLUMNS];
};

#endif // _DPSGeometry_
