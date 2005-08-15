// $Id$
//
//    File: DParameter.cc
// Created: Fri Aug 12 14:58:18 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include "DParameter.h"

//---------------------------------
// DParameter    (Constructor)
//---------------------------------
DParameter::DParameter(string my_key, string my_value)
{
	key = my_key;
	value = my_value;
	isdefault = false;
}

//---------------------------------
// ~DParameter    (Destructor)
//---------------------------------
DParameter::~DParameter()
{

}
