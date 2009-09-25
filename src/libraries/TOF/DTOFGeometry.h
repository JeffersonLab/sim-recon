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
  int NSHORTBARS;       ///> number of short scintillator bars
  float LONGBARLENGTH;  ///> length of the long scintillators
  float SHORTBARLENGTH; ///> length of the short scintillators
  float BARWIDTH;       ///> width of the scintillator bars

    
  float bar2y(int bar, int orientation)  const ///> convert bar number to the
  ///> position of the center of the
  ///> bar in local coordinations
  {
    float y;
    if (bar <= NLONGBARS){
      y = ((NLONGBARS/2.0)+(NSHORTBARS/4.0)-(bar-0.5))*BARWIDTH;
      if (bar > (NLONGBARS/2)) y -= (NSHORTBARS/2.0)*BARWIDTH;
    }
    else{
      y = ((NSHORTBARS/4.0)-(((bar-NLONGBARS+1)/2)-0.5))*BARWIDTH;
    }
    if (orientation == 0) y *= -1.0;
    return y;
  }
  
  
  int y2bar(float x, float y, int orientation) const   ///> convert local x,y to 
  ///> bar number
  {
    int bar;
    if (orientation == 0) {float temp = y; y = -1.0*x; x = temp;}
    if (fabs(y) < (NSHORTBARS/4.0)*BARWIDTH){
      bar = NLONGBARS + (int)((NSHORTBARS/4.0*BARWIDTH-y)/BARWIDTH)*2 + 1;
      if (x > 0) bar++;
    }
    else{
      bar = (int)(((NLONGBARS/2.0+NSHORTBARS/4.0)*BARWIDTH-y)/BARWIDTH) + 1;
      if (bar > NLONGBARS/2) bar -= NSHORTBARS/2;
    }
    return bar;
  }

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "", "%d", NLONGBARS);
			AddString(items, "", "%d", NSHORTBARS);
			AddString(items, "", "%6.3f", LONGBARLENGTH);
			AddString(items, "", "%6.3f", SHORTBARLENGTH);
			AddString(items, "", "%6.3f", BARWIDTH);
		}
};

#endif // _DTOFGeometry_

