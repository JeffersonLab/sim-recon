// $Id$
//
//    File: DDIRCLut_factory.h
//

#ifndef _DDIRCLut_factory_
#define _DDIRCLut_factory_

#include <JANA/JFactory.h>
#include "DDIRCLut.h"

class DDIRCLut_factory:public JFactory<DDIRCLut> {

public:
	
	DDIRCLut_factory(){};
	~DDIRCLut_factory(){};

private:
	jerror_t evnt(jana::JEventLoop *loop, uint64_t eventnumber){
		
		DDIRCLut *dDIRCLut = new DDIRCLut(loop);
		SetFactoryFlag(PERSISTANT);
		ClearFactoryFlag(WRITE_TO_OUTPUT);
		_data.push_back(dDIRCLut);
		
		return NOERROR;
	}
};

#endif // _DDIRCLut_factory_

