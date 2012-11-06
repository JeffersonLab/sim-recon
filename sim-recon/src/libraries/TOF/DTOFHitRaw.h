// $Id: DTOFHitRaw.h Fri Dec  3 08:37:04 EST 2010
//
/// File: DTOFHitRaw.h
/// Created: Fri Dec  3 08:37:04 EST 2010
/// Creator: B. Zihlmann
/// Purpose: Container class to hold Monte Carlo data, unsmeared and 
///          smeared with the MC tag.
//

#ifndef _DTOFHitRaw_
#define _DTOFHitRaw_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DTOFHitRaw:public JObject{
  
 public:
  JOBJECT_PUBLIC(DTOFHitRaw);
  
  int plane;		// plane (0: vertical, 1: horizontal)
  int bar;		// bar number
  int ptype;		// GEANT particle type
  int itrack;           // Track number of primary particle causing the hit
  float t_north;	// time of light at end of bar
  float dE_north;	// attenuated energy deposition
  float t_south;	// time of light at end of bar
  float dE_south;	// attenuated energy deposition
  float dist;           // hit position in paddel (distance from center)
  float x;              // hit location in global coordiantes
  float y;
  float z;
  float px;		// particle momentum
  float py;
  float pz;
  float E;		// particle Energy
  
  void toStrings(vector<pair<string,string> > &items)const{
    AddString(items, "bar", "%d", bar);
    AddString(items, "plane", "%d", plane);
    AddString(items, "t_north", "%12.4e", t_north);
    AddString(items, "dE_north", "%12.4e", dE_north);
    AddString(items, "t_south", "%12.4e", t_south);
    AddString(items, "dE_south", "%12.4e", dE_south);
    AddString(items, "dist", "%12.4e", dist);
    AddString(items, "x", "%12.4e", x);
    AddString(items, "y", "%12.4e", y);
    AddString(items, "z", "%12.4e", z);
    AddString(items, "px", "%12.4e", px);
    AddString(items, "py", "%12.4e", py);
    AddString(items, "pz", "%12.4e", pz);
    AddString(items, "E", "%12.4e", E);
    AddString(items, "ptype", "%d", ptype);
    AddString(items, "itrack", "%d", itrack);
  }
};

#endif // _DTOFHitRaw_

