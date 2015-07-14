//
//    File: MCReconstructionStatuses.h
// Created: Thu Oct 29 09:49:51 EDT 2009
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//

#ifndef _MCReconstructionStatuses_
#define _MCReconstructionStatuses_

#include <TObject.h>
#include <vector>
#include <MCReconstructionStatus.h>

class MCReconstructionStatuses : public TObject{

	public:

		MCReconstructionStatuses(){};
		~MCReconstructionStatuses(){};

		// Data members
		vector<MCReconstructionStatus*> dMCReconstructionStatusVector;

	private:
		ClassDef(MCReconstructionStatuses, 1);

};

#endif // _MCReconstructionStatuses_
