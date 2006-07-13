// $Id$
//
//    File: DMCTrajectoryPoint.h
// Created: Mon Jun 12 09:29:58 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.6.0 powerpc)
//

#ifndef _DMCTrajectoryPoint_
#define _DMCTrajectoryPoint_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DMCTrajectoryPoint:public JObject{
	public:
		HDCLASSDEF(DMCTrajectoryPoint);
		
		float x,y,z,t;
		float px,py,pz;
		float E, dE;
		int track;
		int part;
};

#endif // _DMCTrajectoryPoint_

