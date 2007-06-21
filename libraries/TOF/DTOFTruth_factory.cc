// $Id$
//
//    File: DTOFTruth_factory.cc
// Created: Mon Oct 17 13:58:02 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#include <cassert>	

#include "DTOFTruth_factory.h"


//------------------
// toString
//------------------
const string DTOFTruth_factory::toString(void)
{
  // Ensure our Get method has been called so _data is up to date
  Get();
  if(_data.size()<=0)return string(); // don't print anything if we have no data!
  
  printheader("id:  primary: track:      x:         y:         z:        t:");
	
  for(unsigned int i=0; i<_data.size(); i++){

    DTOFTruth *truth = _data[i];

    printnewrow();
    printcol("%d",	truth->id);
    printcol("%d",	truth->primary);
    printcol("%d",	truth->track);
    printcol("%d",	truth->ptype);
    printcol("%1.3f",	truth->x);
    printcol("%1.3f",	truth->y);
    printcol("%1.3f",	truth->z);
    printcol("%1.3f",	truth->t);
    printcol("%1.3f",	truth->px);
    printcol("%1.3f",	truth->py);
    printcol("%1.3f",	truth->pz);
    printcol("%1.3f",	truth->E);
    printrow();
  }

  return _table;
}
