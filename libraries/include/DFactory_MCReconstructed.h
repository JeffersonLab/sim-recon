// $Id$



#ifndef _DMCRECONSTRUCTED_H_
#define _DMCRECONSTRUCTED_H_

#include "DFactory.h"

typedef struct{
	int type;			///< GEANT particle ID
	float q;				///< electric charge
	float p;				///< Total momentum in GeV/c
	float E;				///< Total energy in GeV
	float theta,phi;	///< Inital theta and phi angles in radians
	float x,y,z;		///< Vertex position in cm
	float mass;			///< Mass in GeV/c^2
}MCReconstructed_t;

class DFactory_MCReconstructed:public DFactory{
	public:
		DFactory_MCReconstructed(DEvent *event):DFactory(event, "MCReconstructed", sizeof(MCReconstructed_t)){};
		~DFactory_MCReconstructed(){};
				
		derror_t Print(void);
	
	private:
		derror_t evnt(int eventnumber);		///< Called every event.
};

#endif //_DMCRECONSTRUCTED_H_

