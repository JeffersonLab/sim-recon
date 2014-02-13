// $Id$
//
//    File: DMCThrown.h
// Created: Sun Apr  3 12:22:09 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DMCThrown_
#define _DMCThrown_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

#include "PID/DKinematicData.h"

class DMCThrown:public DKinematicData{
	public:
		JOBJECT_PUBLIC(DMCThrown);
		
		int type;			///< GEANT particle ID
		int pdgtype;		///< PDG particle type (not used by GEANT)
		int myid;			///< id of this particle from original generator
		int parentid;		///< id of parent of this particle from original generator
		int mech;			///< production mechanism of this partcle (generator specific)

		void toStrings(vector<pair<string,string> > &items)const{
			DKinematicData::toStrings(items);
			AddString(items, "pdgtype", "%d", pdgtype);
			AddString(items, "myid", "%d", myid);
			AddString(items, "parentid", "%d", parentid);
			AddString(items, "mech", "%d", mech);
		}

};

#endif // _DMCThrown_

