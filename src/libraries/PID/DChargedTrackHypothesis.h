// $Id$
//
//    File: DChargedTrackHypothesis_factory.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DChargedTrackHypothesis_
#define _DChargedTrackHypothesis_

#include <vector>
#include <PID/DKinematicData.h>
#include <PID/DDetectorMatches.h>

using namespace std;

class DChargedTrackHypothesis : public DKinematicData
{
	public:
		JOBJECT_PUBLIC(DChargedTrackHypothesis);

		const DTrackTimeBased* dTrackTimeBased; //can get candidateid from here

		unsigned int dNDF_DCdEdx;
		double dChiSq_DCdEdx;
		
		unsigned int dNDF_Timing;
		double dChiSq_Timing;

		unsigned int dNDF; //total NDF used for PID determination
		double dChiSq; //total chi-squared used for PID determination
		double dFOM; //overall FOM for PID determination

		//THESE RETURN NULL IF NO MATCH TO THAT SYSTEM
		shared_ptr<const DSCHitMatchParams> Get_SCHitMatchParams(void) const{return dSCHitMatchParams;}
		shared_ptr<const DTOFHitMatchParams> Get_TOFHitMatchParams(void) const{return dTOFHitMatchParams;}
		shared_ptr<const DBCALShowerMatchParams> Get_BCALShowerMatchParams(void) const{return dBCALShowerMatchParams;}
		shared_ptr<const DFCALShowerMatchParams> Get_FCALShowerMatchParams(void) const{return dFCALShowerMatchParams;}

		void Set_SCHitMatchParams(shared_ptr<const DSCHitMatchParams> locMatchParams){dSCHitMatchParams = locMatchParams;}
		void Set_TOFHitMatchParams(shared_ptr<const DTOFHitMatchParams> locMatchParams){dTOFHitMatchParams = locMatchParams;}
		void Set_BCALShowerMatchParams(shared_ptr<const DBCALShowerMatchParams> locMatchParams){dBCALShowerMatchParams = locMatchParams;}
		void Set_FCALShowerMatchParams(shared_ptr<const DFCALShowerMatchParams> locMatchParams){dFCALShowerMatchParams = locMatchParams;}

		void toStrings(vector<pair<string,string> > &items) const
		{
			AddString(items, "candidate","%d",dTrackTimeBased->candidateid);
			DKinematicData::toStrings(items);
			AddString(items, "Track_ChiSq", "%f", dTrackTimeBased->chisq);
			AddString(items, "dEdx_ChiSq", "%f", dChiSq_DCdEdx);
			AddString(items, "TOF_ChiSq", "%f", dChiSq_Timing);
			AddString(items, "PID_ChiSq", "%f", dChiSq);
			AddString(items, "PID_FOM", "%f", dFOM);
		}

	private:

		shared_ptr<const DSCHitMatchParams> dSCHitMatchParams;
		shared_ptr<const DTOFHitMatchParams> dTOFHitMatchParams;
		shared_ptr<const DBCALShowerMatchParams> dBCALShowerMatchParams;
		shared_ptr<const DFCALShowerMatchParams> dFCALShowerMatchParams;
};

#endif // _DChargedTrackHypothesis_

