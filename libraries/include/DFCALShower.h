// $Id$
//
//    File: DFCALShower.h
// Created: Tue May 17 11:57:50 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#ifndef _DFCALShower_
#define _DFCALShower_

#include "DObject.h"
#include "DFactory.h"

class DFCALShower:public DObject{
	public:
		HDCLASSDEF(DFCALShower);

                float x;  ///< x position of the shower center
                float y;  ///< y position of the shower center
                float E;  ///< Energy of the shower
                float t;  ///< Time of the shower
		
};

#endif // _DFCALShower_

