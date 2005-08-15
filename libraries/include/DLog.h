//****************************************************************
// DLog.h - Error-logging class for use with the Hall D analysis 
// framework.
// 
// At some point, there will be a global log that all subsystems will
// log errors to. For now, in your own subsystem, create yourself a 
// DLog object in one of the following ways:
//
//	DLog* log = new DLog();    // Only display "ERROR"-level messages (on cerr)
//  
//
//	ofstream info ("info.log");
//  ofstream warn ("warn.log");
//	DLog* log = new DLog(info, warn, cerr);  // Put "INFO"-, "WARN"-, and 
//											 // "ERR"-level messages into these
//											 // respective streams.
//
// After that, to log a message to your DLog, use:
// 
//	log("A blankity-blank error occurred.", ERR);
//  
//
// 
// Author: Craig Bookwalter
// July 2005
//
// TODO: Add an ID string to identify who the log belongs to?
// TODO: use DOxygen to do the docs.
// TODO: operator()?
//****************************************************************
#ifndef _DLOG_H_
#define _DLOG_H_

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <string>

enum {
	INFO,
	WARN,
	ERR
};

class DLog
{
	public:
		DLog();
		DLog(std::ostream& info, 
			 std::ostream& warn, 
			 std::ostream& err,
			 int level = 0);
		virtual ~DLog();
		void setLoggingLevel(int level);		
		inline void log(std::string msg, int level) {
			if (level < __level)
				return;
			std::time_t ts = std::time(0);
			char* t = ctime(&ts); 
			t[strlen(t)-1] = 0;
			switch (level) {
				case INFO : 
					*__info << "<INFO " 
							<< t 
							<< "> "
							<< msg 
							<< std::endl;
				case WARN :
					*__warn << "<WARNING " 
							<< t
							<< "> "
							<< msg
							<< std::endl;
				case ERR  :
					*__err  << "<ERROR "
							<< t
							<< "> "
							<< msg
							<< std::endl;
			}
		}
		inline void log(char* msg, int level){
			log(std::string(msg), level);
		}
		inline void operator() (std::string msg, int level){
			log(msg, level);
		}
		inline void operator() (char* msg, int level) {
			log(msg, level);
		}
			
	private:
		std::ostream* __info;
		std::ostream* __warn;
		std::ostream* __err;
		int __level;
};

// Global variable used for logging
extern DLog dlog;

#endif /*_DLOG_H_*/
