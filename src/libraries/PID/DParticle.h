// $Id$
//
//    File: DParticle.h
// Created: Thu Sep  4 14:02:44 EDT 2008
// Creator: davidl (on Darwin harriet.jlab.org 8.11.1 i386)
//

#ifndef _DParticle_
#define _DParticle_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include <PID/DKinematicData.h>

class DReferenceTrajectory;


class DParticle:public DKinematicData{
	public:
		JOBJECT_PUBLIC(DParticle);
		
		oid_t trackid;			///< id of DTrack this came from
		float chisq;			///< Chi-squared for the track (not chisq/dof!)
		int Ndof;				///< Number of degrees of freedom in the fit
		//float dE;				///< Total energy deposited in straws
		//float ds;				///< Total pathlength through straws contributing to dE
		//float err_dE;			///< Error on value of dE
		//float err_ds;			///< Error on value of ds
		const DReferenceTrajectory *rt; ///< pointer to reference trjectory representing this track

		void toStrings(vector<pair<string,string> > &items)const{
			DKinematicData::toStrings(items);
			AddString(items, "trackid", "0x%x", trackid);
			AddString(items, "chisq", "%f", chisq);
			AddString(items, "Ndof", "%d", Ndof);
		}
};

#endif // _DParticle_

