#ifndef _DBCALSiPMSpectrum_
#define _DBCALSiPMSpectrum_

#include <JANA/JObject.h>
using namespace jana;

#include <DHistogram.h>
#include <BCAL/DBCALGeometry.h>

class DBCALSiPMSpectrum:public JObject{

	/// This class holds the signal at the BCAL SiPM as provided by GEANT
	/// (before summing and electronic response).

	public:
		JOBJECT_PUBLIC(DBCALSiPMSpectrum);

		DBCALSiPMSpectrum() : spectrum(4000, -100.0, 300.0) {};

		int module;
		int layer; //This is SiPM layer (ranges from 1-10), not ADC layer (1-4)
		int sector;
		DBCALGeometry::End end;

		DHistogram spectrum;

		int incident_id;

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "module", "%d", module);
			AddString(items, "layer", "%d", layer);
			AddString(items, "sector", "%d", sector);
			AddString(items, "end", "%s", end==0 ? "upstream":"downstream" );
			AddString(items, "incident_id", "%d", incident_id);
		}
};

#endif // _DBCALSiPMSpectrum_
