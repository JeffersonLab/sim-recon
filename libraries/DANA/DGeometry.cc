// $Id$
//
//    File: DGeometry.cc
// Created: Thu Aug 11 22:17:29 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include "DGeometry.h"

//---------------------------------
// DGeometry    (Constructor)
//---------------------------------
DGeometry::DGeometry(unsigned int run_number)
{
	/// The value of run_number is the run for which the geometry is
	/// requested. The values of min_run_number and max_run_number
	/// should be set here to reflect the range of runs for which
	/// this geometry is valid. These are checked using the inline
	/// method IsInRange() from DApplication when a geometry object
	/// is requested.

	/// (Placeholder for now. This should come from the database ...)
	min_run_number = max_run_number = run_number;
}

//---------------------------------
// ~DGeometry    (Destructor)
//---------------------------------
DGeometry::~DGeometry()
{

}
