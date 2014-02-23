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

class DTOFGeometry:public JObject{

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


  float bar2y(int bar, int orientation)  const ///> convert bar number to the
  ///> position of the center of the
  ///> bar in local coordinations
  {
    float y;
    y = YPOS[bar];

    //if (orientation == 0) y *= -1.0;

    return y;
  }
  
  
  int y2bar(float x, float y, int orientation) const   ///> convert local x,y to 
  ///> bar number
  {
    int N1 = NWIDEBARS/2;
    int N2 = (NLONGBARS - NWIDEBARS )/2;
    int N3 = NSHORTBARS/4;
    int N4 = (NLONGBARS - NWIDEBARS )/2;

    int bar;
    if (orientation == 0) {float temp = y; y = -1.0*x; x = temp;}

    if (fabs(y)<N3*BARWIDTH){ // short bars on each side of the beam hole
      
      bar = N1 + N2 + (int)((y+N3*BARWIDTH)/BARWIDTH) ;
      
    } else if (fabs(y)<(N3*BARWIDTH + N2*BARWIDTH/2.)){ // 2 long narrow bars on each side of the beam 
      if (y<0){
	bar = N1+N2 + ((int)((y + N3*BARWIDTH)/(BARWIDTH/2.))-1);
      } else {
	bar = N1+N2+N3+N3 + (int)((y - N3*BARWIDTH)/(BARWIDTH/2.));
      }
      
    } else { // 19 long wide bars on the both sides of the beam
      if (y<0){
	bar = N1 + ((int)((y + N2*BARWIDTH/2 + N3*BARWIDTH)/BARWIDTH)-1);
      } else {
	bar = N1+N2+N3+N3+N4 + ((int)((y - N2*BARWIDTH/2 - N3*BARWIDTH)/BARWIDTH));
      }
    }
    

    return bar;
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

