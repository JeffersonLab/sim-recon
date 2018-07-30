// $Id$
//
//    File: DTOFGeometry.h
// Created: Mon Jul 18 11:43:31 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8 i686)
//

#ifndef _DTOFGeometry_
#define _DTOFGeometry_

#include <HDGEOMETRY/DGeometry.h>
#include <DANA/DApplication.h>

#include <math.h>

#include "JANA/JObject.h"
#include "JANA/JFactory.h"
using namespace jana;

class DTOFGeometry : public JObject {

 public:
  JOBJECT_PUBLIC(DTOFGeometry);
  
  DTOFGeometry(const DGeometry* locGeometry);
  
  // Get functions
  int Get_NLayers() const { return NLAYERS; };
  int Get_NPlanes() const { return NLAYERS; };
  int Get_NEnds() const { return NENDS; };
  // NOTE THAT Get_NBars() is the number of bars along one end!
  // Therefore Get_NBars() != Get_NLongBars()+Get_NShortBars()
  int Get_NBars() const { return NINSTALLBARS; }
  int Get_NLongBars() const { return NLONGBARS; }
  int Get_NShortBars() const { return NSHORTBARS; }

  int Get_FirstShortBar() const { return FirstShortBar; }
  int Get_LastShortBar() const { return LastShortBar; }

  float Get_LongBarLength() const { return LONGBARLENGTH; }
  float Get_HalfLongBarLength() const { return HALFLONGBARLENGTH; }
  float Get_ShortBarLength() const { return SHORTBARLENGTH; }
  float Get_HalfShortBarLength() const { return HALFSHORTBARLENGTH; }
  float Get_BarWidth() const { return BARWIDTH; }

  float Get_CenterVertPlane() const { return CenterVPlane; };  
  float Get_CenterHorizPlane() const { return CenterHPlane; };
  float Get_CenterMidPlane() const { return CenterMPlane; };

  float bar2y(int bar, int end=0) const;  
  int y2bar(double y) const;

  void toStrings(vector<pair<string,string> > &items) const {
		AddString(items, "NBARS", "%d", Get_NBars() );
		AddString(items, "NLONGBARS", "%d",  Get_NLongBars() );
		AddString(items, "NSHORTBARS", "%d", Get_NShortBars() );
		AddString(items, "LONGBARLENGTH", "%6.3f", Get_LongBarLength() );
		AddString(items, "SHORTBARLENGTH", "%6.3f", Get_ShortBarLength() );
		AddString(items, "BARWIDTH", "%6.3f", Get_BarWidth() );
  }
  
 private:
  int NLAYERS;         /// number of scintillator layers
  int NENDS;           /// maximum number of ends that are read out (should be 2!)

  int NLONGBARS;        ///> number of long scintillator bars
  //int NWIDEBARS;        ///> number of wide long+short scintillator bars (deprecated)
  int NSHORTBARS;       ///> number of short scintillator bars
  int NBARS;            ///> number of long scintillator bars
  int NINSTALLBARS;     ///> number of bars vertically = NLONGBARS + NSHORTBARS/2

  int FirstShortBar;    ///> bar number of first short bar
  int LastShortBar;     ///> bar number of last short bar of same type north

  float LONGBARLENGTH;  ///> length of the long scintillators
  float HALFLONGBARLENGTH;  ///> middle of the long scintillators
  float SHORTBARLENGTH; ///> length of the short scintillators
  float HALFSHORTBARLENGTH;  ///> middle of the short scintillators
  float BARWIDTH;       ///> width of the scintillator bars

  float CenterVPlane;  /// center z position of Vertical Plane
  float CenterHPlane;  /// center z position of Horizontal Plane
  float CenterMPlane;  /// center z position between the two Plane

  vector<double> YPOS;  ///> y (perpendicular) position for bar number
 
};

#endif // _DTOFGeometry_

