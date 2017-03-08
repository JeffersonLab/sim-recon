// $Id$
//
//    File: DNeutralParticleHypothesis.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DNeutralParticleHypothesis_
#define _DNeutralParticleHypothesis_

#include <vector>
#include <string>
#include <memory>

#include <JANA/JObject.h>
#include <PID/DKinematicData.h>
#include <PID/DNeutralShower.h>

using namespace std;

class DNeutralParticleHypothesis : public DKinematicData
{
	public:
		JOBJECT_PUBLIC(DNeutralParticleHypothesis);

		//CONSTRUCTORS & OPERATORS
		DNeutralParticleHypothesis(void);
		DNeutralParticleHypothesis(const DNeutralParticleHypothesis& locSourceData, bool locShareTimingFlag = false,
				bool locShareNonKinematicsFlag = false, bool locShareKinematicsFlag = false);
		DNeutralParticleHypothesis& operator= (const DNeutralParticleHypothesis& locSourceData);

		//GETTERS
		const DNeutralShower* Get_NeutralShower(void) const{return dNeutralShower;}

		//Timing
		double t0(void) const{return dTimingInfo->dt0;}
		double t0_err(void) const{return dTimingInfo->dt0_err;}
		DetectorSystem_t t0_detector(void) const{return dTimingInfo->dt0_detector;}
		DetectorSystem_t t1_detector(void) const{return dNeutralShower->dDetectorSystem;}

		//totals for overall PID determination
		unsigned int Get_NDF(void) const{return dTimingInfo->dNDF;}
		double Get_ChiSq(void) const{return dTimingInfo->dChiSq;}
		double Get_FOM(void) const{return dTimingInfo->dFOM;}

		//SETTERS
		void Set_NeutralShower(const DNeutralShower* locNeutralShower){dNeutralShower = locNeutralShower;}

		//Timing
		void Set_T0(double locT0, double locT0Error, DetectorSystem_t locT0Detector);
		void Set_ChiSq_Overall(double locChiSq, unsigned int locNDF, double locFOM);

		void toStrings(vector<pair<string,string> > &items) const
		{
			DKinematicData::toStrings(items);
			AddString(items, "PID_ChiSq", "%f", dChiSq);
			AddString(items, "PID_FOM", "%f", dFOM);
		}

	private:

		struct DTimingInfo
		{
			DTimingInfo(void);

			double dt0;
			double dt0_err;
			DetectorSystem_t dt0_detector;

			float dChiSq;
			unsigned int dNDF;
			float dFOM;
		};

		const DNeutralShower* dNeutralShower;
		shared_ptr<DTimingInfo> dTimingInfo;
};

/************************************************************** CONSTRUCTORS & OPERATORS ***************************************************************/

inline DNeutralParticleHypothesis::DNeutralParticleHypothesis(void) :
dTimingInfo(make_shared<DTimingInfo>()), dTrackingInfo(make_shared<DTrackingInfo>())
{}

inline DNeutralParticleHypothesis::DNeutralParticleHypothesis(const DNeutralParticleHypothesis& locSourceData,
		bool locShareTimingFlag, bool locShareNonKinematicsFlag, bool locShareKinematicsFlag) :
		DKinematicData(locSourceData, locShareNonKinematicsFlag, locShareKinematicsFlag)
{
	//Default is NOT to share: create a new, independent copy of the input data (tracked separately from input so it can be modified)
	dTimingInfo = locShareTimingFlag ? locSourceData->dTimingInfo : make_shared<DTimingInfo>(*(locSourceData->dTimingInfo));
}

inline DNeutralParticleHypothesis& DNeutralParticleHypothesis::operator=(const DNeutralParticleHypothesis& locSourceData)
{
	//Replace current data with a new, independent copy of the input data: tracked separately from input so it can be modified
	dTimingInfo = make_shared<DTimingInfo>(*(locSourceData->dTimingInfo));
	dNeutralShower = locSourceData->dNeutralShower;
}

inline DNeutralParticleHypothesis::DTimingInfo::DTimingInfo(void) :
dt0(0.0), dt0_err(0.0), dt0_detector(SYS_NULL), dNDF(0), dChiSq(0.0), dFOM(0.0)
{}

/********************************************************************** GETTERS ************************************************************************/


/********************************************************************** SETTERS ************************************************************************/

inline void DNeutralParticleHypothesis::Set_T0(double locT0, double locT0Error, DetectorSystem_t locT0Detector)
{
	dTimingInfo->dt0 = locT0;
	dTimingInfo->dt0_err = locT0Error;
	dTimingInfo->dt0_detector = locT0Detector;
}

inline void DNeutralParticleHypothesis::Set_ChiSq_Overall(double locChiSq, unsigned int locNDF, double locFOM)
{
	dTimingInfo->dChiSq = locChiSq;
	dTimingInfo->dNDF = locNDF;
	dTimingInfo->dFOM = locFOM;
}

#endif // _DNeutralParticleHypothesis_
