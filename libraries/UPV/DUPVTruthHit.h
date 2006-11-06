// $Id$
//
//    File: DUPVTruthHit.h
// Created: Mon Nov  6 09:58:35 EST 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.8.0 powerpc)
//

#ifndef _DUPVTruthHit_
#define _DUPVTruthHit_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DUPVTruthHit:public JObject{
	public:
		HDCLASSDEF(DUPVTruthHit);
		
		float	E;
		bool	primary;
		float	t;
		int	track;
		float	x;
		float	y;
		float	z;
};

#endif // _DUPVTruthHit_

