// $Id$
//
//    File: DEventProcessor_mcthrown_tree.h
// Created: Mon Apr  3 11:38:03 EDT 2006
// Creator: pmatt (on Darwin swire-b241.jlab.org 8.4.0 powerpc)
//

#ifndef _DEventProcessor_mcthrown_tree_
#define _DEventProcessor_mcthrown_tree_

#include "JANA/JEventProcessor.h"
#include "DANA/DApplication.h"
#include <ANALYSIS/DEventWriterROOT.h>

using namespace jana;

class DEventProcessor_mcthrown_tree : public JEventProcessor
{
	public:
		DEventProcessor_mcthrown_tree(){};
		~DEventProcessor_mcthrown_tree(){};
		const char* className(void){return "DEventProcessor_mcthrown_tree";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t evnt(JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.

		bool dTagCheck;
};

#endif // _DEventProcessor_mcthrown_tree_

