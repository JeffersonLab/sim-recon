
#ifndef _DLumi_
#define _DLumi_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

#include "TAGGER/DTAGHGeometry.h"
#include "TAGGER/DTAGMGeometry.h"

using namespace jana;
using namespace std;

#include <string>

class DLumi : public JObject {
  
 public:
  
  JOBJECT_PUBLIC(DLumi);
  
  DLumi(JEventLoop *loop);
  ~DLumi();
  
  static const int DETECTORS =  2;
  static const int TAGM_CH   =  102;
  static const int TAGH_CH   =  274;

  double m_psc_accept[3];
  double m_ps_accept[3];
  
  double tagm_tagged[TAGM_CH];
  double tagh_tagged[TAGH_CH];

  double tagm_lumi[TAGM_CH];
  double tagh_lumi[TAGH_CH];

  double Ebeam;

  void CalcLumi();  
  void PrintLumi();
  void SaveLumi();
  
  void CalcTAGHEff();
  void CalcTAGMEff();

 private: 

  vector<const DTAGHGeometry*> taghGeomVect;
  vector<const DTAGMGeometry*> tagmGeomVect;

  int compute_lumi;
  

};

#endif // _DLumi_
