/*
 * DTACHit_factory.h
 *
 *  Created on: Mar 24, 2017
 *      Author: hovanes
 */

#ifndef LIBRARIES_TAC_DTACHIT_FACTORY_H_
#define LIBRARIES_TAC_DTACHIT_FACTORY_H_

#include <vector>
using namespace std;

#include <JANA/JFactory.h>
#include "DTACHit.h"

class DTACHit_factory: public jana::JFactory<DTACHit> {
public:
	DTACHit_factory() {
		// TODO Auto-generated constructor stub

	}
	virtual ~DTACHit_factory() {
		// TODO Auto-generated destructor stub
	}

protected:
   virtual jerror_t init(void) override ;                                          ///< Called once at program start
   virtual jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber) override;    ///< Called everytime a new run number is detected
   virtual jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber) override;  ///< Called every event
   virtual jerror_t erun(void) override;                                          ///< Called everytime run number changes, if brun has been called
   virtual jerror_t fini(void) override;                                          ///< Called after last event of last event source has been processed

   virtual void Reset_Data(void);

};

#endif /* LIBRARIES_TAC_DTACHIT_FACTORY_H_ */
