// $Id$
//
//    File: DPiMinus.h
// Created: Sat Apr 14 12:13:05 EDT 2012
// Creator: davidl (on Darwin genmacbook.local 11.3.0 i386)
//

#ifndef _DPiMinus_
#define _DPiMinus_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include <PID/DKinematicData.h>
#include <TMath.h>

class DPiMinus:public DKinematicData{
	public:
		JOBJECT_PUBLIC(DPiMinus);
	
		DPiMinus(const DKinematicData *kd):DKinematicData(*kd),hypoth(NULL){}
		DPiMinus(const DChargedTrackHypothesis *hypoth):DKinematicData(*hypoth),hypoth(hypoth){}
		
		/// This class is used to hold a copy of the reconstructed
		/// objects for which the pi- hypothesis is the most probable. 
		/// It should not be used for an in-depth physics analysis
		/// because there may be other hypotheses for this track that
		/// are nearly as probable.

		// The following are convienience method for accessing values in the 
		// DChargedTrackHypothesis object on which this is based
		double GetFOM(void) const {
			if(!hypoth)return -1.0;
			return (double)hypoth->dFOM;
		}	
	
		double GetChiSq(void) const {
			if(!hypoth)return -1.0;
			return (double)hypoth->dChiSq;
		}
	
		double GetNDF(void) const {
			if(!hypoth)return -1.0;
			return (double)hypoth->dNDF;
		}	

		double GetProb(void) const {
			if(!hypoth)return -1.0;
			return TMath::Prob((double)hypoth->dChiSq, (int)hypoth->dNDF);
		}
	
		double GetProbTrack(void) const {
			if(!hypoth)return -1.0;
			return TMath::Prob((double)hypoth->dChiSq_Track, (int)hypoth->dNDF_Track);
		}
	
		const DChargedTrackHypothesis* GetChargedTrackHypothesis(void) const {
			return hypoth;
		}
	
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			this->DKinematicData::toStrings(items);
		}
		
	protected:
		const DChargedTrackHypothesis *hypoth;
	
};

#endif // _DPiMinus_

