// $Id$
//
//    File: DDIRCLut_factory.h
//

#ifndef _DDIRCLut_factory_
#define _DDIRCLut_factory_

#include "JANA/JFactory.h"
#include "DDIRCLut.h"

class DDIRCLut_factory : public JFactory<DDIRCLut> {

public:
	
	DDIRCLut_factory() {}
	~DDIRCLut_factory(){}

private:
	
	jerror_t brun(JEventLoop *loop, int32_t runnumber);	
	jerror_t erun(void);	
};

#endif // _DDIRCLut_factory_

