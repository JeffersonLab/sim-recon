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

#include "DResettable.h"
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
		DNeutralParticleHypothesis(const DNeutralParticleHypothesis& locSourceData, bool locShareTimingFlag = false, bool locShareKinematicsFlag = false);
		DNeutralParticleHypothesis& operator= (const DNeutralParticleHypothesis& locSourceData);

		//SHARE RESOURCES
		void Share_FromInput(const DNeutralParticleHypothesis* locSourceData, bool locShareTimingFlag, bool locShareKinematicsFlag);

		void Reset(void);
		void Release(void);

		//GETTERS
		const DNeutralShower* Get_NeutralShower(void) const{return dNeutralShower;}

		//Timing
		double t0(void) const{return dTimingInfo->dt0;}
		double t0_err(void) const{return dTimingInfo->dt0_err;}
		double t1(void) const{return dNeutralShower->dSpacetimeVertex.T();}
		double t1_err(void) const{return (*(dNeutralShower->dCovarianceMatrix))(4, 4);}
		DetectorSystem_t t0_detector(void) const{return dTimingInfo->dt0_detector;}
		DetectorSystem_t t1_detector(void) const{return dNeutralShower->dDetectorSystem;}
		double Get_PathLength(void) const{return (dNeutralShower->dSpacetimeVertex.Vect() - position()).Mag();}
		double measuredBeta(void) const{return ((Get_PathLength()/(t1() - t0())))/29.9792458;}

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
			AddString(items, "PID_ChiSq", "%f", Get_ChiSq());
			AddString(items, "PID_FOM", "%f", Get_FOM());
		}

		class DTimingInfo : public DResettable
		{
			public:
				void Reset(void);
				void Release(void){};

				double dt0 = 0.0;
				double dt0_err = 0.0;
				DetectorSystem_t dt0_detector = SYS_NULL;

				unsigned int dNDF = 0;
				float dChiSq = 0.0;
				float dFOM = 0.0;
		};

	private:
		shared_ptr<DTimingInfo> dTimingInfo = nullptr;
		const DNeutralShower* dNeutralShower = nullptr;

		//RESOURCE POOL
		static thread_local shared_ptr<DResourcePool<DTimingInfo>> dResourcePool_TimingInfo;
};


/************************************************************** CONSTRUCTORS & OPERATORS ***************************************************************/

inline DNeutralParticleHypothesis::DNeutralParticleHypothesis(void) :
dTimingInfo(dResourcePool_TimingInfo->Get_SharedResource()), dNeutralShower(nullptr) {}

inline DNeutralParticleHypothesis::DNeutralParticleHypothesis(const DNeutralParticleHypothesis& locSourceData,
		bool locShareTimingFlag, bool locShareKinematicsFlag) :
		DKinematicData(locSourceData, locShareKinematicsFlag), dNeutralShower(nullptr)
{
	//Default is NOT to share: create a new, independent copy of the input data (tracked separately from input so it can be modified)
	if(locShareTimingFlag)
		dTimingInfo = locSourceData.dTimingInfo;
	else
	{
		dTimingInfo = dResourcePool_TimingInfo->Get_SharedResource();
		*dTimingInfo = *(locSourceData.dTimingInfo);
	}
}

inline DNeutralParticleHypothesis& DNeutralParticleHypothesis::operator=(const DNeutralParticleHypothesis& locSourceData)
{
	//Replace current data with a new, independent copy of the input data: tracked separately from input so it can be modified
	DKinematicData::operator=(locSourceData);
	if(dTimingInfo == locSourceData.dTimingInfo)
		return *this; //guard against self-assignment
	dTimingInfo = dResourcePool_TimingInfo->Get_SharedResource();
	*dTimingInfo = *(locSourceData.dTimingInfo);
	dNeutralShower = locSourceData.dNeutralShower;
	return *this;
}

/********************************************************************** GETTERS ************************************************************************/

inline void DNeutralParticleHypothesis::Share_FromInput(const DNeutralParticleHypothesis* locSourceData, bool locShareTimingFlag, bool locShareKinematicsFlag)
{
	if(locShareTimingFlag)
		dTimingInfo = const_cast<DNeutralParticleHypothesis*>(locSourceData)->dTimingInfo;
	if(locShareKinematicsFlag)
		Share_FromInput_Kinematics(static_cast<const DKinematicData*>(locSourceData));
}

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

inline void DNeutralParticleHypothesis::Reset(void)
{
	DKinematicData::Reset();
	dTimingInfo = dResourcePool_TimingInfo->Get_SharedResource(); //not safe to reset individually, since you don't know what it's shared with
	dNeutralShower = nullptr;
}

inline void DNeutralParticleHypothesis::Release(void)
{
	DKinematicData::Release();
	dTimingInfo = nullptr;
	dNeutralShower = nullptr;
}

inline void DNeutralParticleHypothesis::DTimingInfo::Reset(void)
{
	dt0 = 0.0;
	dt0_err = 0.0;
	dt0_detector = SYS_NULL;
	dNDF = 0;
	dChiSq = 0.0;
	dFOM = 0.0;
}

#endif // _DNeutralParticleHypothesis_
