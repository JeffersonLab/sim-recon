//************************************************
// DFDCHit.h: A class defining a general FDC hit
// Author: Craig Bookwalter, David Lawrence
// Date:	March 2006
//************************************************

#ifndef DFDCHIT_H
#define DFDCHIT_H

#include <sstream>
using namespace std;

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

///
/// class DFDCHit: definition for a basic FDC hit data type.
///
class DFDCHit : public JObject{
	public:
		JOBJECT_PUBLIC(DFDCHit);		
		int layer;			// 1(V), 2(X), or 3(U)
		int module;			// 1 through 8, 1 module = 3 detection layers
		int element;			// wire or strip number
	    int plane;				// for cathodes only: u(3) or v(1) plane, u@+15,v@-15  
	    int gPlane;				// 1 through 72
	    int gLayer;				// 1 through 24
	    float q;				// charge deposited
	    float pulse_height;                 // amplitude of signal
       float pulse_height_raw; //amplitude of signal without gain correction
	    float t;				// drift time
	    float r;				// perpendicular distance from 
	    					// center of chamber to wire/strip center
	    float d;                            // DOCA distance of closest approach (only for MC data on wires)
	    // Enum to take into account split cathode strips near center in 
	    // addition to wires-versus-cathodes
	    enum fdc_hit_type{
	      AnodeWire,
	      FullCathodeStrip,
	      HalfCathodeStripA,
	      HalfCathodeStripB
	    };
	    int type;		// value according to above enum
	    

	    int itrack;                         // track number causing the hit
	    int ptype;                          // particle type causing the hit

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "layer",  "%d",      layer);
			AddString(items, "module", "%d",      module);
			AddString(items, "w/s #",  "%d",      element);
			AddString(items, "plane",  "%s",      plane==1 ? "V":(plane==2 ? "X":"U"));
			AddString(items, "gPlane", "%d",      gPlane);
			AddString(items, "gLayer", "%d",      gLayer);
			AddString(items, "q",      "%10.2f",  q);
			AddString(items, "pulse height","%10.2f", pulse_height);
			AddString(items, "t",      "%10.2f",      t);
			AddString(items, "r",      "%10.2f",      r);
			AddString(items, "d",      "%f",      d);
			AddString(items, "type",   "%d",      type);
			AddString(items, "itrack", "%d",      itrack);
			AddString(items, "ptype",  "%d",      ptype);
		}
	    	
};


#endif // DFDCHIT_H

