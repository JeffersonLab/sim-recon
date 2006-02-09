// $Id$
//
//    File: DFactory_DBCALGeometry.cc
// Created: Thu Nov 17 15:10:51 CST 2005
// Creator: gluexuser (on Linux hydra.phys.uregina.ca 2.4.20-8smp i686)
//

#include <cassert>	

#include "DFactory_DBCALGeometry.h"

//------------------
// evnt
//------------------
derror_t DFactory_DBCALGeometry::evnt(DEventLoop *loop, int eventnumber)
{

  DBCALGeometry *bcalGeom = new DBCALGeometry;

   bcalGeom->NBCALMODS  =48; 
   bcalGeom->NBCALLAYS1 =8;  
   bcalGeom->NBCALLAYS2 =4;  
   bcalGeom->NBCALSECS1 =4;  
   bcalGeom->NBCALSECS2 =3; 
   bcalGeom->BCALINNERRAD =65.0;   
   bcalGeom->BCALMIDRAD   =75.0;   
   bcalGeom->BCALOUTERRAD =90.0;   
   bcalGeom->BCALFIBERLENTH =395.0;   

  _data.push_back(bcalGeom);


	return NOERROR;
}

//------------------
// toString
//------------------
const string DFactory_DBCALGeometry::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)
        return string(); // don't print anything if we have no data!

       printheader("mod: layn1: layn2: secn1: secn2: inr: midr: outr: length:");

      for(unsigned int i = 0; i < _data.size(); i++) {
       DBCALGeometry *s = _data[i];
       printnewrow();
       printcol("%d",s->NBCALMODS);
       printcol("%d",s->NBCALLAYS1);  
       printcol("%d",s->NBCALLAYS2);  
       printcol("%d",s->NBCALSECS1);  
       printcol("%d",s->NBCALSECS2);
       printcol("%6.3f",s->BCALINNERRAD);   
       printcol("%6.3f",s->BCALMIDRAD);   
       printcol("%6.3f",s->BCALOUTERRAD);   
       printcol("%6.3f",s->BCALFIBERLENTH);   
	cout<<"NBCALSECS2= "<<s->NBCALSECS2<<"\n";
       printrow();
     }
     printnewrow();
     printrow();

     return _table;

}
