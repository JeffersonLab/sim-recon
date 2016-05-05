// $Id$
//
//    File: DESDBProvider.h
// Creator: sdobbs
//
// Interface for EventStore database connection
//

#ifndef _DESDBProvider_
#define _DESDBProvider_

#include <iostream>
#include <string>
#include <vector>

#include <JANA/jerror.h>
#include <JANA/JException.h>

#include <DEventStoreDefs.h>

using namespace jana;
using namespace std;


class DESDBProvider {
	public:
		DESDBProvider(string connection_str) {} 
		virtual ~DESDBProvider() {}
		
		virtual bool Open() = 0;
		
		virtual vector<string> GetGrades() = 0;
		virtual vector<string> GetSkims(string timestamp, string grade) = 0;
		virtual vector<string> GetTimestamps(string grade) = 0;

		virtual vector< pair<EventStore::RunRange,int> > GetRunVersions(string timestamp, string grade);	
		virtual vector<int32_t> GetRunList(EventStore::RunRange run_range,
											int graphid, string & view) = 0;
		virtual vector< pair<int32_t,int32_t> > GetRunUidList(EventStore::RunRange run_range,
											  					int graphid, string &view) = 0;
		virtual string GetKeyFileName(int graphid, string &view, 
									  int32_t run, int32_t uid=0) = 0;
		virtual vector< pair<string,string> > GetDataFileNameTypePairs(int graphid, string &view, 
									  						 int32_t run, int32_t uid=0) = 0;
											
	
		// utility functions
		virtual string GetFileName(int32_t fid) = 0;
		virtual int32_t GetFID(string &filename) = 0;
		virtual pair<string,string> GetFileNameAndType(int fid) = 0;

	protected:
};

#endif   // _DESDBProvider_