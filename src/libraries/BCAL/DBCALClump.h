#ifndef _DBCALClump_
#define _DBCALClump_

/*
 *  DBCALClump.h
 *
 *  Created by Beni Zihlmann Tue Mar 12 2013
 *
 */


#include "BCAL/DBCALGeometry.h"
#include "BCAL/DBCALHit.h"

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

#include <vector>

using namespace jana;
using namespace std;

class DBCALClump : public JObject {
  
 public:
  
  JOBJECT_PUBLIC( DBCALClump );
  DBCALClump(vector <const DBCALHit*>, vector <const DBCALHit*>);
  
  vector <const DBCALHit*> HitsU;  // up stream hits of this Clump
  vector <const DBCALHit*> HitsD;  // down stream hits of this Clump
  vector <float> MeanTime;  // list of mean times 
  vector <float> DeltaTime; // list of time differences in [cm]
  vector <int> Sector;
  vector <int> Layer;
  double ProfileU[60];  // up stream profile of the Clump
  double ProfileD[60];  // down stream profile of the Clump 
  double ProfileMT[60]; // average mean time in sector
  double ProfileTD[60]; // average time difference in sector
  
  // the following values are vectors for future upgrade of the code do 
  // take into account possible overlapping showers
  vector <float> ClumpE;     // comined energy of the Clump from up and downstream
  vector <float> ClumpMT;    // mean time of the Clump
  vector <float> ClumpPos;   // Clump position along the BCAL
  vector <float> ClumpPhi;   // azimutal angle of Clump

  void resetProfiles(void);
  void fillArrays(float*, float*);
  void AnalyzeClump();

 private:

  
};

#endif // _DBCALClump_
