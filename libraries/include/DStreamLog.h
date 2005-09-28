//*************************************************************
// DStreamLog.h - Header file for stream-based logging.
// Author:	Craig Bookwalter
// Date:	Aug 2005
// Notes: 	Much of this was derived from examples given by
// 	D. Kuehl at http://www.inf.uni-konstanz.de/~kuehl/iostream/
// 	Also, many thanks to J. Hardie for his assistance with 
//  obscure protected-member rules. 
//*************************************************************


#ifndef _DSTREAMLOG_H_
#define _DSTREAMLOG_H_

#include <iostream>
#include <string>
#include "DStreamLogBuffer.h"

/// DStreamLog provides an intuitive interface for logging errors
/// and messages. To use:
///	
///
/// #include <iostream>
/// #include <fstream>	
///	#include "DStreamLog.h" 
///
/// using namespace std;
///
/// int main() {
///		DStreamLog info(cout, "INFO");
///			
///		info << "Some information. " << endl;
///		
///		ofstream f("errors.log");		
///		DStreamLog err(f, "ERR");
///		
///		err << "There was an error. " << endl;
/// 	return 0;
/// }
///
///	Some details:
/// The instantiation of a DStreamLog requires that an output
/// stream and a tag be provided. To write a single message to
/// a DLogStream object, be sure to pass in std::endl at the end
/// of your text so that the buffer is flushed. 
/// 

class DStreamLog : public std::ostream
{
	public:
		DStreamLog(const std::ostream& os, const char* tag);
		DStreamLog(const std::ostream& os, const std::string& tag);
		DStreamLog(std::streambuf* buf, const char* tag);
		virtual ~DStreamLog();
};

#endif //_DSTREAMLOG_H_
