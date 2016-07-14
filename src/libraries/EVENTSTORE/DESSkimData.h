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
#include <set>
#include <string>
#include <iostream>
using namespace std;

class DESSkimData : public JObject {

  public:
	
	JOBJECT_PUBLIC(DESSkimData);

	DESSkimData(const set<string> &in_event_skims, const vector<string> &in_skim_list) :
			skim_list(in_skim_list), event_skims(in_event_skims) {}
	~DESSkimData() {}

	inline const vector<string>& GetAllSkims() const { return skim_list; } 
	inline const set<string>& GetEventSkims() const { return event_skims; } 
	bool Get_IsEventSkim(string locSkim) const{return (event_skims.find(locSkim) != event_skims.end());}


	void Print(string mode="") const { 
		if(mode=="all") {
			cout << endl << "These skims are available:" << endl;
			//cout << " N = " << skim_list.size() << endl;
			for(vector<string>::const_iterator it = skim_list.begin();
				it != skim_list.end(); it++) 
			 	cout << "  " << *it << endl;
		}
		
		cout << endl << "Event satisfies these skims:" << endl;
		//cout << " N = " << event_skims.size() << endl;
		for(set<string>::const_iterator it = event_skims.begin();
			it != event_skims.end(); it++) 
			 cout << "  " << *it << endl;
	} 

  protected:
	vector<string> skim_list;
	set<string> event_skims;
};

#endif  // _DESSkimData_
