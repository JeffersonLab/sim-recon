// $Id$
//
//    File: DBCALGeometry.h
// Created: Thu Nov 17 15:10:51 CST 2005
// Creator: gluexuser (on Linux hydra.phys.uregina.ca 2.4.20-8smp i686)
//

#ifndef _DBCALGeometry_
#define _DBCALGeometry_

#include "DObject.h"
#include "DFactory.h"

class DBCALGeometry:public DObject{
	public:
		HDCLASSDEF(DBCALGeometry);

		int NBCALMODS; ///> number of modules
		int NBCALLAYS1; ///> number of layers in first 10 ccm 
		int NBCALLAYS2; ///> number of layers in last  15 cm 
		int NBCALSECS1; ///> number of sectors in first 10cm of Mod 
		int NBCALSECS2; ///> number of sectors in last 15cm of Mod 
                float BCALINNERRAD;   ///> innner radius of BCAL in cm
                float BCALMIDRAD;   ///> mid radius of BCAL in cm
                float BCALOUTERRAD;   ///> outer radius of BCAL in cm
                float BCALFIBERLENTH;   ///> BCAL Scintilator fiber lenth in cm

};

#endif // _DBCALGeometry_

