//************************************************************
// DStreamLogBuffer.h - streambuf-derived buffer for use with
// DStreamLog.
// Author:	Craig Bookwalter
// Date: 	Aug 2005
// Notes:	Much of this was derived from examples given by
// 	D. Kuehl at http://www.inf.uni-konstanz.de/~kuehl/iostream/
// 	Also, many thanks to J. Hardie for his assistance with 
//  obscure protected-member rules. 
//*************************************************************

#ifndef _DSTREAMLOGBUFFER_H_
#define _DSTREAMLOGBUFFER_H_

#include <iostream>


/// DStreamLogBuffer is a streambuf-derived buffer for use
/// with DLogStream. Unless you're writing your own stream,
/// you shouldn't ever need this class. It's basically a 
/// wrapper for a passed-in streambuf with the tweak that
/// it appends a prefix to each output statement, namely
/// a status label and a timestamp. 
///

class DStreamLogBuffer : public std::streambuf
{
	private:
		std::streambuf* __sbuf;
		char*			__tag;
		bool			__newline;
		
	protected:
		int overflow(int c);
		int sync();
		const char* getTimeStamp();
	
	public:
		DStreamLogBuffer(std::streambuf* buf, const char* tag);
		virtual ~DStreamLogBuffer();
};

#endif //_DSTREAMLOGBUFFER_H_
