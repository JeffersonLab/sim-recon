// $Id$
//
//    File: DChargedTrackHypothesis_factory.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DChargedTrackHypothesis_
#define _DChargedTrackHypothesis_

#include <vector>
#include <string>
#include <memory>

#include <PID/DKinematicData.h>
#include <PID/DDetectorMatches.h>

using namespace std;

class DChargedTrackHypothesis : public DKinematicData
{
	public:
		JOBJECT_PUBLIC(DChargedTrackHypothesis);

		//CONSTRUCTORS & OPERATORS
		DChargedTrackHypothesis(void);
		DChargedTrackHypothesis(const DChargedTrackHypothesis& locSourceData, bool locShareTrackingFlag = false,
				bool locShareTimingFlag = false, bool locShareNonKinematicsFlag = false, bool locShareKinematicsFlag = false);
		DChargedTrackHypothesis(const DTrackTimeBased* locSourceData);
		DChargedTrackHypothesis& operator= (const DChargedTrackHypothesis& locSourceData);

		//GETTERS

		//Tracking
		unsigned int Get_NDF_DCdEdx(void) const{return dTrackingInfo->dNDF_DCdEdx;}
		double Get_ChiSq_DCdEdx(void) const{return dTrackingInfo->dChiSq_DCdEdx;}
		const DTrackTimeBased* Get_TrackTimeBased(void) const{return dTrackingInfo->dTrackTimeBased;}

		//Timing
		double t0(void) const{return dTimingInfo->dt0;}
		double t0_err(void) const{return dTimingInfo->dt0_err;}
		DetectorSystem_t t0_detector(void) const{return dTimingInfo->dt0_detector;}
		DetectorSystem_t t1_detector(void) const;
		unsigned int Get_NDF_Timing(void) const{return dTimingInfo->dNDF_Timing;}
		double Get_ChiSq_Timing(void) const{return dTimingInfo->dChiSq_Timing;}

		//totals for overall PID determination
		unsigned int Get_NDF(void) const{return dTimingInfo->dNDF;}
		double Get_ChiSq(void) const{return dTimingInfo->dChiSq;}
		double Get_FOM(void) const{return dTimingInfo->dFOM;}

		//Match params (return nullptr if no match)
		shared_ptr<const DSCHitMatchParams> Get_SCHitMatchParams(void) const{return dTrackingInfo->dSCHitMatchParams;}
		shared_ptr<const DTOFHitMatchParams> Get_TOFHitMatchParams(void) const{return dTrackingInfo->dTOFHitMatchParams;}
		shared_ptr<const DBCALShowerMatchParams> Get_BCALShowerMatchParams(void) const{return dTrackingInfo->dBCALShowerMatchParams;}
		shared_ptr<const DFCALShowerMatchParams> Get_FCALShowerMatchParams(void) const{return dTrackingInfo->dFCALShowerMatchParams;}

		//SETTERS

		//Timing
		void Set_T0(double locT0, double locT0Error, DetectorSystem_t locT0Detector);
		void Set_ChiSq_Timing(double locChiSq, unsigned int locNDF);
		void Set_ChiSq_Overall(double locChiSq, unsigned int locNDF, double locFOM);

		//Tracking
		void Set_TrackTimeBased(const DTrackTimeBased* locTrackTimeBased){dTrackingInfo->dTrackTimeBased = locTrackTimeBased;}
		void Set_ChiSq_Tracking(double locChiSq, unsigned int locNDF);

		//Match params
		void Set_SCHitMatchParams(shared_ptr<const DSCHitMatchParams> locMatchParams){dTrackingInfo->dSCHitMatchParams = locMatchParams;}
		void Set_TOFHitMatchParams(shared_ptr<const DTOFHitMatchParams> locMatchParams){dTrackingInfo->dTOFHitMatchParams = locMatchParams;}
		void Set_BCALShowerMatchParams(shared_ptr<const DBCALShowerMatchParams> locMatchParams){dTrackingInfo->dBCALShowerMatchParams = locMatchParams;}
		void Set_FCALShowerMatchParams(shared_ptr<const DFCALShowerMatchParams> locMatchParams){dTrackingInfo->dFCALShowerMatchParams = locMatchParams;}

		void toStrings(vector<pair<string,string> > &items) const
		{
			AddString(items, "candidate","%d", dTrackingInfo->dTrackTimeBased->candidateid);
			DKinematicData::toStrings(items);
			AddString(items, "Track_ChiSq", "%f", dTrackingInfo->dTrackTimeBased->chisq);
			AddString(items, "dEdx_ChiSq", "%f", dTrackingInfo->dChiSq_DCdEdx);
			AddString(items, "TOF_ChiSq", "%f", dTimingInfo->dChiSq_Timing);
			AddString(items, "PID_ChiSq", "%f", dTimingInfo->dChiSq);
			AddString(items, "PID_FOM", "%f", dTimingInfo->dFOM);
		}

	private:

		struct DTimingInfo
		{
			DTimingInfo(void);

			double dt0;
			double dt0_err;
			DetectorSystem_t dt0_detector;

			unsigned int dNDF_Timing;
			double dChiSq_Timing;

			//technically, these can depend on the tracking chisq also, but no one in their right mind would change the tracking dE/dx info
			unsigned int dNDF; //total NDF used for PID determination
			double dChiSq; //total chi-squared used for PID determination
			double dFOM; //overall FOM for PID determination
		};

		struct DTrackingInfo
		{
			DTrackingInfo(void);

			unsigned int dNDF_DCdEdx;
			double dChiSq_DCdEdx;

			const DTrackTimeBased* dTrackTimeBased; //can get candidateid from here

			shared_ptr<const DSCHitMatchParams> dSCHitMatchParams;
			shared_ptr<const DTOFHitMatchParams> dTOFHitMatchParams;
			shared_ptr<const DBCALShowerMatchParams> dBCALShowerMatchParams;
			shared_ptr<const DFCALShowerMatchParams> dFCALShowerMatchParams;
		};

		//memory of object in shared_ptr is managed automatically: deleted automatically when no references are left
		shared_ptr<DTimingInfo> dTimingInfo;
		shared_ptr<DTrackingInfo> dTrackingInfo;
};

/************************************************************** CONSTRUCTORS & OPERATORS ***************************************************************/

inline DChargedTrackHypothesis::DChargedTrackHypothesis(void) :
dTimingInfo(make_shared<DTimingInfo>()), dTrackingInfo(make_shared<DTrackingInfo>())
{}

inline DChargedTrackHypothesis::DChargedTrackHypothesis(const DChargedTrackHypothesis& locSourceData, bool locShareTrackingFlag,
		bool locShareTimingFlag, bool locShareNonKinematicsFlag, bool locShareKinematicsFlag) :
		DKinematicData(locSourceData, locShareNonKinematicsFlag, locShareKinematicsFlag)
{
	//Default is NOT to share: create a new, independent copy of the input data (tracked separately from input so it can be modified)
	dTrackingInfo = locShareTrackingFlag ? locSourceData->dTrackingInfo : make_shared<DTrackingInfo>(*(locSourceData->dTrackingInfo));
	dTimingInfo = locShareTimingFlag ? locSourceData->dTimingInfo : make_shared<DTimingInfo>(*(locSourceData->dTimingInfo));
}

inline DChargedTrackHypothesis::DChargedTrackHypothesis(const DTrackTimeBased* locSourceData) :
		DKinematicData(*static_cast<DKinematicData*>(locSourceData), true, false)
{
	//Default is NOT to share: create a new, independent copy of the input data (tracked separately from input so it can be modified)
	dTrackingInfo = make_shared<DTrackingInfo>();
	dTimingInfo = make_shared<DTimingInfo>();
	dTrackingInfo->dTrackTimeBased = locSourceData;
}

inline DChargedTrackHypothesis& DChargedTrackHypothesis::operator=(const DChargedTrackHypothesis& locSourceData)
{
	//Replace current data with a new, independent copy of the input data: tracked separately from input so it can be modified
	dTimingInfo = make_shared<DTimingInfo>(*(locSourceData->dTimingInfo));
	dTrackingInfo = make_shared<DTrackingInfo>(*(locSourceData->dTrackingInfo));
}

inline DChargedTrackHypothesis::DTimingInfo::DTimingInfo(void) :
dt0(0.0), dt0_err(0.0), dt0_detector(SYS_NULL), dNDF_Timing(0), dChiSq_Timing(0.0), dNDF(0), dChiSq(0.0), dFOM(0.0)
{}

inline DChargedTrackHypothesis::DTrackingInfo::DTrackingInfo(void) :
dNDF_DCdEdx(0), dChiSq_DCdEdx(0.0), dTrackTimeBased(nullptr),
dSCHitMatchParams(nullptr), dTOFHitMatchParams(nullptr), dBCALShowerMatchParams(nullptr), dFCALShowerMatchParams(nullptr)
{}

/********************************************************************** GETTERS ************************************************************************/

inline DetectorSystem_t DChargedTrackHypothesis::t1_detector(void) const
{
	if(Get_BCALShowerMatchParams() != nullptr)
		return SYS_BCAL;
	else if(Get_TOFHitMatchParams() != nullptr)
		return SYS_TOF;
	else if(Get_FCALShowerMatchParams() != nullptr)
		return SYS_FCAL;
	return SYS_NULL;
}

/********************************************************************** SETTERS ************************************************************************/

inline void DChargedTrackHypothesis::Set_T0(double locT0, double locT0Error, DetectorSystem_t locT0Detector)
{
	dTimingInfo->dt0 = locT0;
	dTimingInfo->dt0_err = locT0Error;
	dTimingInfo->dt0_detector = locT0Detector;
}

inline void DChargedTrackHypothesis::Set_ChiSq_Timing(double locChiSq, unsigned int locNDF)
{
	dTimingInfo->dChiSq_Timing = locChiSq;
	dTimingInfo->dNDF_Timing = locNDF;
}

inline void DChargedTrackHypothesis::Set_ChiSq_Tracking(double locChiSq, unsigned int locNDF)
{
	dTrackingInfo->dChiSq_DCdEdx = locChiSq;
	dTrackingInfo->dNDF_DCdEdx = locNDF;
}

inline void DChargedTrackHypothesis::Set_ChiSq_Overall(double locChiSq, unsigned int locNDF, double locFOM)
{
	dTimingInfo->dChiSq = locChiSq;
	dTimingInfo->dNDF = locNDF;
	dTimingInfo->dFOM = locFOM;
}

#endif // _DChargedTrackHypothesis_
