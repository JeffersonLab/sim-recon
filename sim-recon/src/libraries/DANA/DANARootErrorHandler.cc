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


static int ROOT_ERROR_LEVEL_SUPRESS=-1;

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
	
		// Check for certain messages we wish to ignore no matter what
		string sloc(loc);
		string smsg(msg);
		
		// The following is an error that occurs on occasion in the
		// tracking code while inverting a matrix. The condition is
		// checked for and handled there so printing of an error is
		// unnecessary (and even a bit misleading since it is not 
		// really an error condition)
		if(sloc == "TDecompLU::DecomposeLUCrout" && smsg == "matrix is singular") return;

		// The DRootGeom class generates several "errors" in when
		// loading the ROOT geometry such as from the following:
		//    TGeoVolumeMulti::CheckShapes
		//    TGeoEltu::GetMakeRuntimeShape
		//     ... etc...
		// These may real issues with the ROOT geometry, but it is
		// not used for core reconstruction at this time so the
		// errors are distracting.
		if(sloc.find("TGeo") == 0) return;

		// Use normal ROOT error handler 
		DefaultErrorHandler(lvl, abt, loc, msg);
	}
}

