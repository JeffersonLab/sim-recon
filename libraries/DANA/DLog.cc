#include "DLog.h"
#include <stdexcept>
#include <time.h>

DLog::DLog() {
	__info  = &std::cout;
	__warn  = &std::cout;
	__err   = &std::cerr;
	__level = ERR;
}

DLog::DLog(std::ostream& info, 
		   std::ostream& warn, 
		   std::ostream& err, 
		   int level) 
{
	__info  = &info;
	__warn  = &warn;
	__err   = &err;
	__level = level;
}

DLog::~DLog() {}

void DLog::setLoggingLevel(int level) {
	__level = level;
}
