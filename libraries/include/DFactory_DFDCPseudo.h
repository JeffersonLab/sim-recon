#ifndef DFACTORY_DFDCPSEUDO_H
#define DFACTORY_DFDCPSEUDO_H

#include "DFactory.h"
#include "DFDCPseudo.h"
#include "DFDCHit.h"
#include "DException.h"
#include "DFDCGeometry.h"
#include "DStreamLog.h"

#include <algorithm>
#include <map>
#include <cmath>

class DFactory_DFDCPseudo : public DFactory<DFDCPseudo> {
	public:
		DFactory_DFDCPseudo();
		~DFactory_DFDCPseudo();
		void conjure(	vector<const DFDCHit*>& u, 
						vector<const DFDCHit*>& v,
						map<int, const DFDCHit*>& x,
						float angle);
		void dummy(	vector<const DFDCHit*>& u, 
		  			vector<const DFDCHit*>& v,
					map<int, const DFDCHit*>& x,
					float angle,
					int evNo,
					int layerNo);
		const string toString(void) {return "";}
	
	protected:
		//derror_t init(void);						///< Called once at program start.
		//derror_t brun(DEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		derror_t evnt(DEventLoop *eventLoop, int eventNo);	///< Called every event.
		//derror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		//derror_t fini(void);						///< Called after last event of last event source has been processed.
		
	private:
		DFDCGeometry _geo;
		DStreamLog* _log;
		ofstream* logFile;
};

#endif // DFACTORY_DFDCPSEUDO_H

