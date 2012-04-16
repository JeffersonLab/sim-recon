// $Id$
//
//    File: DPiPlus.h
// Created: Sat Apr 14 12:13:05 EDT 2012
// Creator: davidl (on Darwin genmacbook.local 11.3.0 i386)
//

#ifndef _DPiPlus_
#define _DPiPlus_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include <PID/DKinematicData.h>

class DPiPlus:public DKinematicData{
	public:
		JOBJECT_PUBLIC(DPiPlus);
	
		DPiPlus(const DKinematicData *kd):DKinematicData(*kd){}
		
		/// This class is used to hold a copy of the reconstructed
		/// objects for which the pi+ hypothesis is the most probable. 
		/// It should not be used for an in-depth physics analysis
		/// because there may be other hypotheses for this track that
		/// are nearly as probable.
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			this->DKinematicData::toStrings(items);
		}
		
};

#endif // _DPiPlus_

