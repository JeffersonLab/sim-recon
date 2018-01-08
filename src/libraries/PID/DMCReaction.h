// $Id$
//
//    File: DMCReaction.h
// Created: Sun Aug 28 18:41:08 EDT 2011
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//

#ifndef _DMCReaction_
#define _DMCReaction_

#include <JANA/jerror.h>
#include <JANA/JObject.h>

#include <PID/DKinematicData.h>

class DMCReaction:public JObject{
	public:
		DMCReaction(){}
		virtual ~DMCReaction(){}
		JOBJECT_PUBLIC(DMCReaction);
		
		int type;
		double weight;
		DKinematicData target;
		DKinematicData beam;
		
		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "type", "%2d", type);
			AddString(items, "weight", "%3.1f", weight);
			AddString(items, "mass target(GeV)", "%3.1f", target.mass());
			AddString(items, "energy beam(GeV/c^2)", "%f", beam.energy());
		}

	protected:
	
	
	private:

};

#endif // _DMCReaction_

