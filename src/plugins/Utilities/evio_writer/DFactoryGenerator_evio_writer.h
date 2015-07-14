#ifndef _DFactoryGenerator_evio_writer_
#define _DFactoryGenerator_evio_writer_

#include <JANA/jerror.h>
#include <JANA/JFactoryGenerator.h>

#include "DEventWriterEVIO_factory.h"

class DFactoryGenerator_evio_writer : public jana::JFactoryGenerator
{
	public:
		DFactoryGenerator_evio_writer(){}
		virtual ~DFactoryGenerator_evio_writer(){}
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DFactoryGenerator_evio_writer";}
		
		jerror_t GenerateFactories(jana::JEventLoop *loop)
		{
			loop->AddFactory(new DEventWriterEVIO_factory());
			return NOERROR;
		}
};

#endif // _DFactoryGenerator_evio_writer_

