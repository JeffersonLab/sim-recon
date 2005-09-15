// $Id$
//
//    File: DException.cc
// Created: Wed Sep 14 07:42:51 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include "DException.h"

//---------------------------------
// DException    (Constructor)
//---------------------------------
DException::DException(const char* filename, int line, string message)
{
	this->filename = filename;
	this->line = line;
	this->message = message;
}

//---------------------------------
// ~DException    (Destructor)
//---------------------------------
DException::~DException()
{

}
