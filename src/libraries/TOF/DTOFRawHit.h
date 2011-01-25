// $Id: DTOFRawHit.h Tue Jan 18 16:15:26 EST 2011
//
/// File:    DTOFRawHit.h
/// Created: Tue Jan 18 16:15:26 EST 2011
/// Creator: B. Zihlmann
/// Purpose: Container class to hold Monte Carlo data, unsmeared and 
///          smeared with the MC tag.
//

#ifndef _DTOFRawHit_
#define _DTOFRawHit_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DTOFRawHit:public JObject{
  
 public:
  JOBJECT_PUBLIC(DTOFRawHit);
  
  int plane;		// plane (0: vertical, 1: horizontal)
  int bar;		// bar number
  int lr;               // left/right 0/1 or North/South 0/1
  float dE;      	// attenuated energy deposition
  float t;	        // time of light at end of bar
  
  void toStrings(vector<pair<string,string> > &items)const{
    AddString(items, "bar", "%d", bar);
    AddString(items, "plane", "%d", plane);
    AddString(items, "lr", "%d", lr);
    AddString(items, "dE", "%12.4e", dE);
    AddString(items, "t", "%12.4e", t);
  }
};

#endif // _DTOFRawHit_

