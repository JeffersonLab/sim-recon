// $Id$
//
//    File: DEPICSStore.h
// Created: Tue Nov 18 15:44:17 EST 2014
// Creator: sdobbs (on Linux ifarm1102 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#ifndef _DEPICSSTORE_H_
#define _DEPICSSTORE_H_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>
#include <JANA/JCalibration.h>
#include <JANA/JException.h>
//using namespace jana;

#include <stdint.h>
#include <string>
#include <vector>
#include <map>
using namespace std;

#include <DAQ/DEPICSvalue.h>


struct DEPICSvalue_data {
DEPICSvalue_data():value(NULL),first_time(0),cumulative_average(0.) {}

	DEPICSvalue *value;
	
	time_t first_time;
	double cumulative_average;
};
typedef struct DEPICSvalue_data DEPICSvalue_data_t;


class DEPICSstore : public jana::JObject {
public:
	JOBJECT_PUBLIC(DEPICSstore);
	//DEPICSstore() {}

	const DEPICSvalue *GetValue(string key);
	const double GetAverage(string name);
	void AddValue(const DEPICSvalue *new_value);
	
	vector<string> GetNames();
	map<string, DEPICSvalue_data_t> &GetStore() { return stored_values; }

	void ClearAverages();
	void ResetStartTimes(time_t new_start_time);
	

private:
	map<string, DEPICSvalue_data_t> stored_values;
	map<string, string> pretty_names;                   // name mapping for pretty printing
};

#endif   // _DEPICSSTORE_H_
