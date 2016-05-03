// $Id$
//
//    File: DFactoryGenerator_DESSkimData.h

#ifndef _DFactoryGenerator_DESSkimData_
#define _DFactoryGenerator_DESSkimData_

#include <JANA/jerror.h>
#include <JANA/JFactory.h>
#include <JANA/JFactoryGenerator.h>

#include "DESSkimData.h"

class DFactoryGenerator_DESSkimData : public jana::JFactoryGenerator
{
        public:
                virtual const char* className(void){return static_className();}
                static const char* static_className(void){return "DFactoryGenerator_DESSkimData";}
                
                jerror_t GenerateFactories(jana::JEventLoop* locEventLoop)
                {
                        locEventLoop->AddFactory(new JFactory<DESSkimData>());
                        return NOERROR;
                }
};

#endif // _DFactoryGenerator_DESSkimData_
