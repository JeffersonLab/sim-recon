// $Id$
//
//    File: DTOFGeometry.h
// Created: Mon Jul 18 11:43:31 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8 i686)
//

#ifndef _DTOFGeometry_
#define _DTOFGeometry_

#include <math.h>

#include "JANA/JObject.h"
#include "JANA/JFactory.h"
using namespace jana;

class DTOFGeometry : public JObject {

 public:
  JOBJECT_PUBLIC(DTOFGeometry);

  int NLONGBARS;        ///> number of long scintillator bars
  int NWIDEBARS;        ///> number of long scintillator bars
  int NBARS;            ///> number of long scintillator bars
  int NSHORTBARS;       ///> number of short scintillator bars
  float LONGBARLENGTH;  ///> length of the long scintillators
  float SHORTBARLENGTH; ///> length of the short scintillators
  float BARWIDTH;       ///> width of the scintillator bars
  float YPOS[50];       ///> y position for bar number
  int FirstShortBar;    ///> bar number of first short bar
  int LastShortBar;     ///> bar number of last short bar of same type north

  float CenterVPlane;  /// center z position of Vertical Plane
  float CenterHPlane;  /// center z position of Horizontal Plane
  float CenterMPlane;  /// center z position between the two Plane

  int NLAYERS;         /// number of scintillator layers
  int NENDS;           /// maximum number of ends that are read out (should be 2!)

  float bar2y(int bar, int end=0)  const ///> convert bar number to the
  ///> position of the center of the
  ///> bar in local coordinations
  {
    float y;
    y = YPOS[bar];

    if (bar>=FirstShortBar && bar<=LastShortBar && end != 0) y *= -1.0;

    return y;
  }
  
  
  int y2bar(double y) const   ///> convert local position y to bar number
  {
    int jm=1;
    if (y>YPOS[44]) jm=44;
    else if (y>YPOS[1]){
      int jl=-1;
      int ju=44;
      while(ju-jl>1){
	jm=(ju+jl)>>1;
        if (y>=YPOS[jm])
	  jl=jm;
        else
	  ju=jm;
      }     
      if (fabs(y-YPOS[jm-1])<fabs(y-YPOS[jm])) return jm-1;
      if (fabs(y-YPOS[jm+1])<fabs(y-YPOS[jm])) return jm+1;
    }

    return jm;
  }

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "NBARS", "%d", NBARS);
			AddString(items, "NLONGBARS", "%d", NLONGBARS);
			AddString(items, "NWIDEBARS", "%d", NWIDEBARS);
			AddString(items, "NSHORTBARS", "%d", NSHORTBARS);
			AddString(items, "LONGBARLENGTH", "%6.3f", LONGBARLENGTH);
			AddString(items, "SHORTBARLENGTH", "%6.3f", SHORTBARLENGTH);
			AddString(items, "BARWIDTH", "%6.3f", BARWIDTH);
		}
};

#endif // _DTOFGeometry_

