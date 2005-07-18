// $Id$
//
//    File: DTOFGeometry.h
// Created: Mon Jul 18 11:43:31 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8 i686)
//

#ifndef _DTOFGeometry_
#define _DTOFGeometry_

#include "DFactory.h"

class DTOFGeometry{
	public:
		HDCLASSDEF(DTOFGeometry);

                int NLONGBARS;        ///> number of long scintillator bars
                int NSHORTBARS;       ///> number of short scintillator bars
                float LONGBARLENGTH;  ///> length of the long scintillators
                float SHORTBARLENGTH; ///> length of the short scintillators
                float BARWIDTH;       ///> width of the scintillator bars

		
};

#endif // _DTOFGeometry_

