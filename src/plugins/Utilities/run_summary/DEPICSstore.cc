// $Id$
//
//    File: DEPICSStore.cc
// Created: Tue Nov 18 15:44:17 EST 2014
// Creator: sdobbs (on Linux ifarm1102 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#include <DEPICSstore.h>

//----------------
// GetValue
//   get last read EPICS value, if we have come across it
//----------------
const DEPICSvalue *DEPICSstore::GetValue(string name)
{
	map<string, DEPICSvalue_data_t>::iterator result_itr = stored_values.find(name);
	if(result_itr == stored_values.end())
		return NULL;
	else
		return result_itr->second.value;
}

const double DEPICSstore::GetAverage(string name)
{
	map<string, DEPICSvalue_data_t>::iterator result_itr = stored_values.find(name);
	if(result_itr == stored_values.end())
		return -1.;
	else
		return result_itr->second.cumulative_average;
}

void DEPICSstore::AddValue(const DEPICSvalue *new_value)
{
	map<string, DEPICSvalue_data_t>::iterator result_itr = stored_values.find(new_value->name);
	if(result_itr == stored_values.end()) {
		stored_values[new_value->name] = DEPICSvalue_data_t();
		stored_values[new_value->name].value = new DEPICSvalue(*new_value);
		stored_values[new_value->name].first_time = new_value->timestamp;
	} else {
		stored_values[new_value->name].cumulative_average += stored_values[new_value->name].value->fval / static_cast<double>( new_value->timestamp - stored_values[new_value->name].value->timestamp );
		if(stored_values[new_value->name].value != NULL)
			delete stored_values[new_value->name].value;
		stored_values[new_value->name].value = new DEPICSvalue(*new_value);
	}
}

void DEPICSstore::ClearAverages()
{
	for(map<string, DEPICSvalue_data_t>::iterator val_itr = stored_values.begin();
	    val_itr != stored_values.end(); val_itr++) {
		val_itr->second.cumulative_average = 0.;
	}
}


void DEPICSstore::ResetStartTimes(time_t new_start_time)
{
	for(map<string, DEPICSvalue_data_t>::iterator val_itr = stored_values.begin();
	    val_itr != stored_values.end(); val_itr++) {
		val_itr->second.first_time = new_start_time;
	}
}


vector<string> DEPICSstore::GetNames()
{
	vector<string> names;
	for(map<string, DEPICSvalue_data_t>::iterator val_itr = stored_values.begin();
	    val_itr != stored_values.end(); val_itr++) {
		names.push_back(val_itr->first);
	}
	return names;
}
