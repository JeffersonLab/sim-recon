
#include "DStreamLog.h"

DStreamLog::DStreamLog(std::streambuf* buf, const char* tag) :
std::ostream(new DStreamLogBuffer(buf, tag))
{}

DStreamLog::DStreamLog(const std::ostream& os, const char* tag) :
std::ostream(new DStreamLogBuffer(os.rdbuf(), tag))
{}

DStreamLog::DStreamLog(const std::ostream& os, const std::string& tag) :
std::ostream(new DStreamLogBuffer(os.rdbuf(), tag.c_str()))
{}

DStreamLog::~DStreamLog() {
	delete rdbuf();
}

std::ostream& endMsg(std::ostream& dSL) {
	dSL << static_cast<char>(6); 
	return dSL;
}
