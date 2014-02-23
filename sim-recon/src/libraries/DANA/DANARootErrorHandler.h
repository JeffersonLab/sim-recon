// $Id$
//
//    File: DANARootErrorHandler.h
// Created: Mon Nov  9 16:09:14 EST 2009
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//

#ifndef _DANARootErrorHandler_
#define _DANARootErrorHandler_

#include <iostream>

#include <JANA/jerror.h>

void InitDANARootErrorHandler(int my_ROOT_ERROR_LEVEL_SUPRESS=10000);
void DANARootErrorHandler(int lvl, bool abt, const char* loc, const char* msg);


#endif // _DANARootErrorHandler_

