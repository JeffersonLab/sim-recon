// $Id$
//
//    File: DFactory_DTOFGeometry.cc
// Created: Mon Jul 18 11:43:31 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8 i686)
//

#include "DFactory_DTOFGeometry.h"

derror_t DFactory_DTOFGeometry::init(void)
{

//   This is the geometry as it currently stands inside HDGeant.
//   This does not represent the current design, but will be
//   used just to get things going.  Eventually these will be changed
//   to:
//
//     nlongbars=40
//     nshortbars=4
//     longbarlength=252.0
//     shortbarlength=120.0
//     barwidth=6.0


    DTOFGeometry *myDTOFGeometry = new DTOFGeometry;

    myDTOFGeometry->NLONGBARS        = 42;
    myDTOFGeometry->NSHORTBARS       = 2;
    myDTOFGeometry->LONGBARLENGTH    = 258.0;
    myDTOFGeometry->SHORTBARLENGTH   = 126.0;
    myDTOFGeometry->BARWIDTH         = 6.0;

    _data.push_back(myDTOFGeometry);

    return NOERROR;
}  


//------------------
// evnt
//------------------
derror_t DFactory_DTOFGeometry::evnt(DEventLoop *loop, int eventnumber)
{
	// Code to generate factory data goes here. Add it like:
	//
	// DTOFGeometry *myDTOFGeometry = new DTOFGeometry;
	// myDTOFGeometry->x = x;
	// myDTOFGeometry->y = y;
	// ...
	// _data.push_back(myDTOFGeometry);
	//
	// Note that the objects you create here will be deleted later
	// by the system and the _data vector will be cleared automatically.

	return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DTOFGeometry::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
    Get();
    if(_data.size()<=0)return string(); // don't print anything if we have no data!

	// Put the class specific code to produce nicely formatted ASCII here.
	// The DFactory_base class has several methods defined to help. They
	// rely on positions of colons (:) in the header. Here's an example:
	//

    printheader("NLONGBARS:  NSHORTBARS:  LONGBARLENGTH:  SHORTBARLENGTH:  BARWIDTH:");
	
    for(unsigned int i=0; i<_data.size(); i++){
        DTOFGeometry *myDTOFGeometry = _data[i];
        printnewrow();
        printcol("%d",myDTOFGeometry->NLONGBARS);
        printcol("%d",myDTOFGeometry->NSHORTBARS);
        printcol("%6.3f",myDTOFGeometry->LONGBARLENGTH);
        printcol("%6.3f",myDTOFGeometry->SHORTBARLENGTH);
        printcol("%6.3f",myDTOFGeometry->BARWIDTH);
        printrow();
    }
	
    return _table;

}
