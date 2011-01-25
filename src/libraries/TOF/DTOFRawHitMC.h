// $Id: DTOFHitRawMC.h Wed Jan 19 14:22:41 EST 2011
//
/// File:    DTOFRawHitMC.h
/// Created: Wed Jan 19 14:22:41 EST 2011
/// Creator: B. Zihlmann
/// Purpose: Container class to hold Monte Carlo track data, 
///          like track id number, particle type ect.
//

#ifndef _DTOFRawHitMC_
#define _DTOFRawHitMC_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DTOFRawHitMC:public JObject{
  
 public:
  JOBJECT_PUBLIC(DTOFRawHitMC);
  
  int plane;		// plane (0: vertical, 1: horizontal)
  int bar;		// bar number
  int lr;               // north/south or left/right
  int ptype;		// GEANT particle type
  int itrack;           // Track number of primary particle causing the hit
  float dist;           // Hit distance from center of paddle (or x=0)
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
    AddString(items, "lr", "%d", lr);
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

#endif // _DTOFRawHitMC_

