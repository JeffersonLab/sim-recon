// $Id$
// $HeadURL$
//
//    File: DEPICSvalue.h
// Created: Tue Nov 11 21:40:17 EST 2014
// Creator: davidl (on Linux gluon104.jlab.org 2.6.32-358.23.2.el6.x86_64 x86_64)
//

#ifndef _DEPICSvalue_
#define _DEPICSvalue_

#include <string>
#include <time.h>
#include <stdint.h>
using std::string;

#include <JANA/jerror.h>
#include <JANA/JObject.h>

/// A DEPICSvalue object holds information for a single
/// EPICS value read from the data stream. Values are
/// stored in the data stream as strings of the form:
/// key=value. This string is available in the "nameval"
/// field. It is parsed however so that one may access the
/// name via "name" and the value as a sting via "sval".
/// The "ival", "uval", and "fval" hold the value converted
/// into a an "int", "uint32_t", and "double" respectively.
/// This is done using the stringstream class an is done
/// for convience. The EPICS values are inserted into the
/// EVIO file during data taking by the epics2et program
/// which should be started automatically by the DAQ system.
/// Source for that can be found here:
///
///  https://halldsvn.jlab.org/repos/trunk/online/packages/etUtils/src/epics2et


class DEPICSvalue:public jana::JObject{
	public:
		JOBJECT_PUBLIC(DEPICSvalue);
		DEPICSvalue(time_t timestamp, string &nameval){
			this->timestamp = timestamp;
			this->nameval = nameval;
			size_t pos = nameval.find("=");
			if(pos != nameval.npos){
				name = nameval.substr(0, pos);
				sval = nameval.substr(pos+1);
				
				stringstream ss(sval);
				ss >> ival;
				ss.str(sval);
				ss >> uval;
				ss.str(sval);
				ss >> fval;
				
				// stringstream will set fval to 0.0 if no decimal
				// point is in the string. Use ival in these cases
				if( (fval==0.0) && (ival!=0.0) ) fval = (double)ival;
			}
		}
		virtual ~DEPICSvalue(){}
		
		time_t   timestamp;
		string   nameval;
		string   name;
		string   sval;
		int      ival;
		uint32_t uval;
		double   fval;

		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			string timestr = ctime(&timestamp);
			timestr[timestr.length()-1] = 0;
			AddString(items, "timestamp", "%d", timestamp);
			AddString(items, "name", "%s", name.c_str());
			AddString(items, "sval", "%s", sval.substr(0, 255).c_str());
			AddString(items, "ival", "%d", ival);
			AddString(items, "fval", "%f", (float)fval);
			AddString(items, "t", "%s", timestr.c_str());
		}
};

#endif // _DEPICSvalue_

