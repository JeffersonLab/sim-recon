// $Id$

#ifndef _DBCALUnifiedHit_
#define _DBCALUnifiedHit_

#include "BCAL/DBCALGeometry.h"

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

class DBCALUnifiedHit : public JObject{

//One DBCALUnifiedHit is created for each DBCALHit. When available, TDC hits
//are also incorporated. The class provides energy in GeV rather than ADC
//units (for now the conversion is just a constant factor) and timewalk
//corrections are also applied. With this class, these corrections only need
//to be applied once.

	public:
		JOBJECT_PUBLIC(DBCALUnifiedHit);

		int module;
		int layer;
		int sector;
		DBCALGeometry::End end;
		float E;

		//If there is a associated TDC hit, t is the timewalk-corrected TDC
		//time, otherwise t is the same as t_ADC.
		float t;           ///<  Unified time, obtained from ADC and/or TDC and used for further analysis
		float t_ADC;       ///<  Time from fADC 
		float t_TDC;       ///<  Time of TDC hit that is closes t to the ADC time
		bool has_TDC_hit;  ///<  Flag if the Unified Time is the TDC time
		
		int cellId;

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "module", "%d", module);
			AddString(items, "layer", "%d", layer);
			AddString(items, "sector", "%d", sector);
			AddString(items, "end", "%s", end==0 ? "upstream":"downstream" );
			AddString(items, "E(GeV)", "%2.3f", E);
			AddString(items, "t(ns)", "%4.2f", t);
			AddString(items, "t_ADC(ns)", "%4.2f", t_ADC);
			AddString(items, "t_TDC(ns)", "%4.2f", t_TDC);
		}
};

#endif // _DBCALUnifiedHit_

