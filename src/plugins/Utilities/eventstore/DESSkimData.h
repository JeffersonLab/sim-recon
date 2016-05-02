// $Id$
//
//    File: DESSkimData.h
// Creator: sdobbs 
//


#ifndef _DESSkimData_
#define _DESSkimData_

#include <JANA/JFactory.h>
#include <JANA/JObject.h>
#include <JANA/JEvent.h>
using namespace jana;

#include <vector>
#include <string>
#include <iostream>
using namespace std;

class DESSkimData : public JObject {

  public:
	
	JOBJECT_PUBLIC(DESSkimData);

	DESSkimData(JEvent &event, vector<string> &in_skim_list, int in_base_skim_index);
	~DESSkimData() {}

	inline const vector<string>& GetEventSkims() const { return event_skims; } 
	inline const vector<string>& GetAllSkims() const { return skim_list; } 

	void Print(string mode="") const { 
		if(mode=="all") {
			cout << endl << "These skims are available:" << endl;
			//cout << " N = " << skim_list.size();
			for(vector<string>::const_iterator it = skim_list.begin();
				it != skim_list.end(); it++) 
			 	cout << "  " << *it << endl;
		}
		
		cout << endl << "Event satisfies these skims:" << endl;
		//cout << " N = " << event_skims.size();
		for(vector<string>::const_iterator it = event_skims.begin();
			it != event_skims.end(); it++) 
			 cout << "  " << *it << endl;
	} 

  protected:
	// We tag which skims JEvents belong to using JEvent::SetStatusBit()
	// We can get away with this now, since no one else is using fields above 16 yet
	// Probably the scheme needs to change or we need our own fields
	int BASE_SKIM_INDEX;           // the first status bit that we use for EventStore
	int MAX_SKIM_INDEX;            // we can store 64 bits in the JEvent, so 64 - MAX_SKIM_INDEX

	vector<string> skim_list;
	vector<string> event_skims;
};

#endif  // _DESSkimData_