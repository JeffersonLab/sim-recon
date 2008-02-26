// $Id$
//
//    File: DApplication.h
// Created: Mon Jul  3 21:46:01 EDT 2006
// Creator: davidl (on Darwin Harriet.local 8.6.0 powerpc)
//

#ifndef _DApplication_
#define _DApplication_

#include "JANA/jerror.h"
#include "JANA/JApplication.h"

class DMagneticFieldMap;
class DLorentzDeflections;

class DApplication:public JApplication{

/// The DApplication class extends the JApplication class
/// by adding the default event source generators and
/// factory generators that are HAll-D specific.

	public:
		DApplication(int narg, char* argv[]);
		virtual ~DApplication();
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DApplication";}
		
		DMagneticFieldMap* GetBfield(void){return bfield;}
		DLorentzDeflections *GetLorentzDeflections(void){
		  return lorentz_def;
		}
								

	protected:
	
		DMagneticFieldMap *bfield;
		DLorentzDeflections *lorentz_def;
	private:

};

#endif // _DApplication_

