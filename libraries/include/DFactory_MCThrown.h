// $Id$



#ifndef _DMCTHROWN_H_
#define _DMCTHROWN_H_

#include "DFactory.h"

typedef struct{
	int type;			///< GEANT particle ID
	float q;				///< electric charge
	float p;				///< Total momentum in GeV/c
	float E;				///< Total energy in GeV
	float theta,phi;	///< Inital theta and phi angles in radians
	float x,y,z;		///< Vertex position in cm
	float mass;			///< Mass in GeV/c^2
}MCThrown_t;

class DFactory_MCThrown:public DFactory{
	public:
		DFactory_MCThrown(DEvent *event):DFactory(event, "MCThrown", sizeof(MCThrown_t)){};
		~DFactory_MCThrown(){};
				
		derror_t Print(void);
	
	private:
		derror_t evnt(int eventnumber);		///< Called every event.
};

#endif //_DMCTHROWN_H_

