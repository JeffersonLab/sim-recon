// TODO: Fix the way the time output looks--make the damned thing all on one 
// line
// TODO: Add an ID string to identify who the log belongs to?

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
			
	private:
		std::ostream* __info;
		std::ostream* __warn;
		std::ostream* __err;
		int __level;
};

#endif /*_DLOG_H_*/
