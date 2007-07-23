/*
 *  DMCBCALHit_factory.cc
 *  Hall D
 *
 *  Created by Matthew Shepherd on 7/23/07.
 *
 */

#include "BCAL/DMCBCALHit_factory.h"

const string DMCBCALHit_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!
    
	printheader("id:   module:  layer:    sector:      E:      t:    zLocal:");
    
	for(unsigned int i = 0; i < _data.size(); i++) {
		DHDDMBCALHit *s = _data[i];
        
		printnewrow();
		printcol("%i",s->id);
		printcol("%5i",s->module);
		printcol("%5i",s->layer);
		printcol("%5i",s->sector);
		printcol("%4.3f",s->E);
		printcol("%4.3f",s->t);
        printcol("%4.3f",s->zLocal);
		printrow();
	}
    
    printnewrow();
    printrow();
    
    return _table;
    
}
