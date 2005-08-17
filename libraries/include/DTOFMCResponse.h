// $Id$
//
//    File: DTOFMCResponse.h
// Created: Mon Aug 15 11:33:45 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8 i686)
//

#ifndef _DTOFMCResponse_
#define _DTOFMCResponse_

#include "DFactory.h"

class DTOFMCResponse:public DObject{
	public:
		HDCLASSDEF(DTOFMCResponse);
		
		float x;
		float y;
		float dE;
		float t;
		int orientation;	///< 0=vertical  1=horizontal
		int end;				///< 0=left/top 1=right/bottom
};

#endif // _DTOFMCResponse_

