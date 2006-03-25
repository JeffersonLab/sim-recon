/// 
/// DException - exception definition for DANA
/// Author:  Craig Bookwalter (craigb at jlab.org) 
/// Date:    December 2005
/// Usage: 
/// 	
///     #include <DException.h>
///
///     void someFunction() {
///	        ...
///         throw DException("Error details.");
///		}
///
///	 Notes:	
///   - you must compile with the -g option (g++) to get readable output
///	  - you must be running an executable in the current directory or from
///	  some directory on your path in order for exceptions to function fully
///
///	 To do:
///	  - protect against executables that cannot be located
///   - respond intelligently to executables not compiled with debug info
///   - test and adapt to Solaris platform
///	  - make it thread-safe 	
///

#ifndef DEXCEPTION_H
#define DEXCEPTION_H

#include <iostream>
#include <exception>
#include <cstdlib>
#include <string>
#include <sstream>
#include <new>

#include <execinfo.h>
#include <unistd.h>
#include <limits.h>

class DException : public std::exception 
{
	public :
		DException(std::string msg="");
		virtual	~DException() throw();
		const char* what() const throw();
		friend std::ostream& 
			operator<<(std::ostream& os, const DException& d);
	private:
		std::string _msg;
		std::string _trace;
};	
	
#endif //DEXCEPTION_H
