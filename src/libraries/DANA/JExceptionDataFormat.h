// $Id$
//
//    File: JExceptionDataFormat.h
// Created: Tue Mar  6 15:57:29 EST 2018
// Creator: davidl (on Linux gluon119.jlab.org 2.6.32-642.3.1.el6.x86_64 x86_64)
//

#ifndef _JExceptionDataFormat_
#define _JExceptionDataFormat_

#include <JANA/jerror.h>
#include <JANA/JException.h>

/// This is a subclass of JException that is used to indicate a
/// parsing error. This was motivated by hdmon needing to catch
/// this specific type of error and set an alarm. hdmon is used
/// for online monitoring and the source is kept in subversion.

class JExceptionDataFormat: public JException{
	public:
		JExceptionDataFormat(const std::string &txt):JException(txt){}
		JExceptionDataFormat(const std::string &txt, const char *file, int line):JException(txt, file, line){}
		virtual ~JExceptionDataFormat(){}
		
	protected:

};

#endif // _JExceptionDataFormat_

