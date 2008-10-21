#ifndef _DFDCPseudo_factory_CORRECTED_ 
#define _DFDCPseudo_factory_CORRECTED_ 

#include <JANA/JFactory.h>
using namespace jana;

#include "DFDCPseudo.h"

class DFDCPseudo_factory_CORRECTED:public JFactory<DFDCPseudo>{
	public:
		DFDCPseudo_factory_CORRECTED(){};
		~DFDCPseudo_factory_CORRECTED(){}; 
		const char* Tag(void){return "CORRECTED";}

	private:

};

#endif // _DFDCPseudo_factory_CORRECTED_

