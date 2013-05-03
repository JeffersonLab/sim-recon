#ifndef _DEventWriterREST_
#define _DEventWriterREST_

#include <math.h>
#include <vector>

#include <HDDM/hddm_r.hpp>

#include <JANA/JObject.h>
#include <JANA/JEventLoop.h>
#include <JANA/JApplication.h>

#include <DVector3.h>
#include <DMatrixDSym.h>
#include <DMatrix.h>

#include "PID/DMCReaction.h"
#include "TRACKING/DMCThrown.h"
#include "TAGGER/DTagger.h"
#include "FCAL/DFCALShower.h"
#include "PID/DNeutralShower.h"
#include "BCAL/DBCALShower.h"
#include "TOF/DTOFPoint.h"
#include "START_COUNTER/DSCHit.h"
#include "TRACKING/DTrackTimeBased.h"
#include "TRIGGER/DMCTrigger.h"

using namespace std;
using namespace jana;

class DEventWriterREST : public JObject
{
	public:
		JOBJECT_PUBLIC(DEventWriterREST);

		DEventWriterREST(JEventLoop* locEventLoop, string locOutputFileBaseName);
		~DEventWriterREST(void);

		bool Write_RESTEvent(JEventLoop* locEventLoop, string locOutputFileNameSubString) const;
		string Get_OutputFileName(string locOutputFileNameSubString) const;

	private:
		bool Write_RESTEvent(string locOutputFileName, hddm_r::HDDM& locRecord) const;

		string dOutputFileBaseName;
};

#endif //_DEventWriterREST_


