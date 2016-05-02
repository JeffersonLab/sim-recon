// $Id$
//
//    File: DESSkimData.cc
// Creator: sdobbs 
//

#include "DESSkimData.h"


//---------------------------------
// DESSkimData
//---------------------------------
DESSkimData::DESSkimData(JEvent &event, vector<string> &in_skim_list, int in_base_skim_index) 
{
	// initialize members
	skim_list = in_skim_list;
	BASE_SKIM_INDEX = in_base_skim_index;
	MAX_SKIM_INDEX = 64 - BASE_SKIM_INDEX;
	if(skim_list.size() < MAX_SKIM_INDEX)
		MAX_SKIM_INDEX = BASE_SKIM_INDEX + skim_list.size();

	//cout << "TEST" << endl;
	//for(int i=0; i<skim_list.size(); i++)
	//	cout << i << " " << skim_list[i] << endl;

	// scan through the available status bits, to see if this event is
	// a member of any skims
	for(int i = BASE_SKIM_INDEX; i < MAX_SKIM_INDEX; i++) {
		//cout << i << " " << skim_list[i - BASE_SKIM_INDEX] << endl;
		if(event.GetStatusBit(i))
			event_skims.push_back(skim_list[i - BASE_SKIM_INDEX]);
	}
}
