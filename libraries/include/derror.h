//
/// derror.h
/// 
/// This file contains error codes for errors specific to the
/// Hall-D analysis code. Many Hall-D functions return values
/// of type derror_t. This header should be included in all
/// files which must deal with this type.

#ifndef _DERROR_H_
#define _DERROR_H_


enum derror_t{
	NOERROR = 0,
	UNKNOWN_ERROR = -1000,
	
	MAX_EVENT_PROCESSORS_EXCEEDED
};

#endif //_DERROR_H_

