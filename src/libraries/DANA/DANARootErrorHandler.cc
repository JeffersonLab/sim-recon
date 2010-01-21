// $Id$
//
//    File: DANARootErrorHandler.cc
// Created: Mon Nov  9 16:09:14 EST 2009
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//

// This is based on info found here:
//
// http://root.cern.ch/root/roottalk/roottalk02/0514.html
//

#include <iostream>
using namespace std;

#include <TError.h>
#include "DANARootErrorHandler.h"


static int ROOT_ERROR_LEVEL_SUPRESS=10000;

//---------------------------------
// InitDANARootErrorHandler
//---------------------------------
void InitDANARootErrorHandler(int my_ROOT_ERROR_LEVEL_SUPRESS)
{
	ROOT_ERROR_LEVEL_SUPRESS = my_ROOT_ERROR_LEVEL_SUPRESS;

	// Register handler with ROOT
	SetErrorHandler(DANARootErrorHandler);
}
 
//---------------------------------
// DANARootErrorHandler
//---------------------------------
void DANARootErrorHandler(int lvl, bool abt, const char* loc, const char* msg)
{
	// Here is the entry point for the error handler

	if(lvl>ROOT_ERROR_LEVEL_SUPRESS || abt){
		// Use normal ROOT error handler 
		DefaultErrorHandler(lvl, abt, loc, msg);
	}
}

