// $Id$
//
//    File: DEventWriterROOT_factory_CLcomp.h
// Created: Thu Feb  1 13:48:42 EST 2018
// Creator: aebarnes (on Linux egbert 2.6.32-696.13.2.el6.x86_64 x86_64)
//

#ifndef _DEventWriterROOT_factory_CLcomp_
#define _DEventWriterROOT_factory_CLcomp_

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>

#include "DEventWriterROOT_CLcomp.h"

class DEventWriterROOT_factory_CLcomp : public jana::JFactory<DEventWriterROOT>
{
	public:
		DEventWriterROOT_factory_CLcomp(){use_factory = 1;}; //prevents JANA from searching the input file for these objects
		~DEventWriterROOT_factory_CLcomp(){};
		const char* Tag(void){return "CLcomp";}

	private:
		jerror_t brun(jana::JEventLoop *locEventLoop, int locRunNumber)
		{
			// Create single DEventWriterROOT_CLcomp object and marks the factory as persistent so it doesn't get deleted every event.
			SetFactoryFlag(PERSISTANT);
			ClearFactoryFlag(WRITE_TO_OUTPUT);

			_data.push_back(new DEventWriterROOT_CLcomp());
			_data.back()->Initialize(locEventLoop);
			return NOERROR;
		}
};

#endif // _DEventWriterROOT_factory_CLcomp_
