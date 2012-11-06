#ifndef _DEventProcessor_DCdEdxStudy_tree_
#define _DEventProcessor_DCdEdxStudy_tree_

#include <JANA/JEventProcessor.h>
using namespace jana;

#include <TFile.h>
#include <TTree.h>
#include <DVector3.h>
#include <particleType.h>

#include <DANA/DApplication.h>
#include <TRACKING/DMCThrown.h>
#include <TRACKING/DTrackTimeBased.h>
#include <PID/DParticleID.h>
#include <DCdEdxInformation.h>

class DEventProcessor_DCdEdxStudy_tree:public JEventProcessor{
	public:
		DEventProcessor_DCdEdxStudy_tree(){};
		~DEventProcessor_DCdEdxStudy_tree(){};
		const char* className(void){return "DEventProcessor_DCdEdxStudy_tree";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		DCdEdxInformation *dDCdEdxInformation;
		TTree* dPluginTree_DCdEdxInformation;

		DParticleID *dPIDAlgorithm;
};

#endif // _DEventProcessor_DCdEdxStudy_tree_

