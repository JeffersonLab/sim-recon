//********************************************************************* 
// DFDCCathodeCluster.h: defintion for a cathode strip distribution
// Author: Craig Bookwalter (craigb at jlab.org)
// Date: Apr 2006
//*********************************************************************

#ifndef DFDCCATHODECLUSTER_H
#define DFDCCATHODECLUSTER_H

#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

#include "DFDCHit.h"
#include <JANA/JObject.h>
using namespace jana;

#define HIT_TIME_DIFF_MIN 10.0

class DFDCCathodeCluster : public JObject {
 public:
  JOBJECT_PUBLIC(DFDCCathodeCluster);		/// DANA identifier
  
  vector<const DFDCHit*> members; ///< DFDCHits that make up this cluster
  int plane;			  ///< V=1, U=3
  int gLayer;   ///< #1-24, which detection layer of all FDC  
  int gPlane;	///< #1-74, which plane out of all FDC modules
  float q_tot;		   ///< total energy/charge deposited in the cluster
		
  /// Return a sensible string representation of this object
  void toStrings(vector<pair<string,string> > &items)const{
    AddString(items, "Nmembers", "%d", members.size());
    AddString(items, "plane", "%d", plane);
    AddString(items, "gLayer", "%d", gLayer); 
    AddString(items, "gPlane", "%d", gPlane);
    AddString(items, "q_tot", "%f", q_tot);
  }
};

#endif //DFDCCATHODECLUSTER_H
