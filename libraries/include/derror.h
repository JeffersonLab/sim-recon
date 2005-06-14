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
	
	MAX_EVENT_PROCESSORS_EXCEEDED,
	
	ERROR_OPENING_EVENT_SOURCE,
	ERROR_CLOSING_EVENT_SOURCE,
	NO_MORE_EVENTS_IN_SOURCE,
	NO_MORE_EVENT_SOURCES,
	EVENT_NOT_IN_MEMORY,
	EVENT_SOURCE_NOT_OPEN,
	OBJECT_NOT_AVAILABLE,
	DEVENT_OBJECT_DOES_NOT_EXIST,
	
	MEMORY_ALLOCATION_ERROR,
	
	RESOURCE_UNAVAILABLE,
	VALUE_OUT_OF_RANGE,
	
	INFINITE_RECURSION,
	UNRECOVERABLE_ERROR,
	
	FILTER_EVENT_OUT
};

#endif //_DERROR_H_

