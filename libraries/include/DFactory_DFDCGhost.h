#ifndef DFACTORY_DFDCGHOST_H
#define DFACTORY_DFDCGHOST_H

#include "DFactory.h"
#include "DEventLoop.h"
#include "DFDCGhost.h"
#include "DFDCHit.h"
#include "DException.h"
#include "DFDCGeometry.h"

#include <algorithm>
#include <map>
#include <queue>
#include <cmath>

class DFactory_DFDCGhost : public DFactory<DFDCGhost> {
	public:
		DFactory_DFDCGhost();
		~DFactory_DFDCGhost();
		derror_t Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v);
		void conjure(	const vector<const DFDCHit*>& u, 
						const vector<const DFDCHit*>& v,
						const map<const int, const DFDCHit*>& x,
						float angle);
		const string toString(void);
	
	protected:
		//derror_t init(void);						///< Called once at program start.
		//derror_t brun(DEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		derror_t evnt(DEventLoop *eventLoop, int eventNo);	///< Called every event.
		//derror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		//derror_t fini(void);						///< Called after last event of last event source has been processed.
		
	private:
		DFDCGeometry _geo;
			
};

#endif // DFACTORY_DFDCGHOST_H

