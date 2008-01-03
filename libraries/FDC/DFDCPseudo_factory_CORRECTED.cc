#include "DFDCPseudo_factory_CORRECTED.h"

// Factory for creating Lorentz-corrected pseudopoints.  Currently the actual
// work takes place in the DFDCSegment factory.


//------------------
// toString
//------------------
const string DFDCPseudo_factory_CORRECTED::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("layer: wire: time(ns):      w(cm):     s(cm):   status:");

	for(unsigned int i=0; i<_data.size(); i++){
		DFDCPseudo *hit = _data[i];
		//const DFDCWire *w = hit->wire;

		printnewrow();
		printcol("%d",		hit->wire->layer);
		printcol("%d",		hit->wire->wire);
		printcol("%3.1f",	hit->time);
		printcol("%3.1f",	hit->w);
		printcol("%1.4f",	hit->s);
		printcol("%d",		hit->status);
		printrow();
	}

	return _table;
}
