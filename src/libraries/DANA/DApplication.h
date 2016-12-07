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

#include "TMatrixFSym.h"

#include <JANA/jerror.h>
#include <JANA/JApplication.h>
using namespace jana;
using namespace std;

#include "HDGEOMETRY/DGeometry.h"

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
		
		DMagneticFieldMap* GetBfield(unsigned int run_number=1);
		DLorentzDeflections *GetLorentzDeflections(unsigned int run_number=1);
		DGeometry* GetDGeometry(unsigned int run_number);
		DRootGeom *GetRootGeom(unsigned int run_number);

		size_t Get_NumCovarianceMatrices(void);
		TMatrixFSym* Get_CovarianceMatrixResource(unsigned int locNumMatrixRows);
		TMatrixFSym* Get_CovarianceMatrixResource(unsigned int locNumMatrixRows, uint64_t locEventNumber);
		void Recycle_CovarianceMatrices(const deque<const TMatrixFSym*>& locMatrices);

	protected:
	
		DMagneticFieldMap *bfield;
		DLorentzDeflections *lorentz_def;
		JEventSourceGenerator *event_source_generator;
		JFactoryGenerator *factory_generator;
	 	DRootGeom *RootGeom;	
		vector<DGeometry*> geometries;

		pthread_mutex_t mutex;
		pthread_mutex_t matrix_mutex;
		
		size_t dTargetMaxNumAvailableMatrices;
		map<pthread_t, uint64_t> dEventNumberMap;
		map<pthread_t, set<TMatrixFSym*> > dUsedMatrixMap;
		deque<TMatrixFSym*> dAvailableMatrices;
};

#endif // _DApplication_
