// $Id$
//
//    File: DApplication.h
// Created: Mon Jul  3 21:46:01 EDT 2006
// Creator: davidl (on Darwin Harriet.local 8.6.0 powerpc)
//

#ifndef _DApplication_
#define _DApplication_

#include <deque>
#include <map>

#include <JANA/jerror.h>
#include <JANA/JApplication.h>
using namespace jana;
using namespace std;

#include "HDGEOMETRY/DGeometry.h"
#include <HDDM/hddm_r.hpp>

class DMagneticFieldMap;
class DLorentzDeflections;
class DGeometry;
//class DMaterialMap;
class DRootGeom;

class DApplication:public JApplication{

/// The DApplication class extends the JApplication class
/// by adding the default event source generators and
/// factory generators that are HAll-D specific.

	public:
		DApplication(int narg, char* argv[]);
		virtual ~DApplication();
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DApplication";}
		
		jerror_t Init(void);
		
		DMagneticFieldMap* GetBfield(void);
		DLorentzDeflections *GetLorentzDeflections(void);
		//DMaterialMap *GetMaterialMap(void){return material;}
		DGeometry* GetDGeometry(unsigned int run_number);
		DRootGeom *GetRootGeom();

		//ONLY CALL THESE FUNCTIONS FROM WITHIN A "RESTWriter" WRITE LOCK!
		inline bool Find_RESTOutputFilePointers(string locOutputFileName) const{return (dRESTOutputFilePointers.find(locOutputFileName) != dRESTOutputFilePointers.end());}
		inline void Set_RESTOutputFilePointers(string locOutputFileName, pair<ofstream*, hddm_r::ostream*> locOutputFilePointers){dRESTOutputFilePointers[locOutputFileName] = locOutputFilePointers;}
		inline pair<ofstream*, hddm_r::ostream*> Get_RESTOutputFilePointers(string locOutputFileName){return dRESTOutputFilePointers[locOutputFileName];}
		pair<ofstream*, hddm_r::ostream*> Get_RESTOutputFilePointers(string locOutputFileName) const;
		void Erase_RESTOutputFilePointers(string locOutputFileName);
		void Get_OpenRESTOutputFileNames(deque<string>& locFileNames) const;

	protected:
	
		DMagneticFieldMap *bfield;
		DLorentzDeflections *lorentz_def;
		JEventSourceGenerator *event_source_generator;
		JFactoryGenerator *factory_generator;
  		//DMaterialMap *material;
	 	DRootGeom *RootGeom;	
		vector<DGeometry*> geometries;

		//REST Output Globals
		map<string, pair<ofstream*, hddm_r::ostream*> > dRESTOutputFilePointers;

		pthread_mutex_t mutex;
		
	private:

};

#endif // _DApplication_

