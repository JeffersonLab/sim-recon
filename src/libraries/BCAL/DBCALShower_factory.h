// $Id$

#ifndef _DBCALShower_factory_
#define _DBCALShower_factory_

#include <JANA/JFactory.h>
#include <BCAL/DBCALShower.h>

class DBCALShower_factory:public jana::JFactory<DBCALShower>{
	public:
		DBCALShower_factory(){};
		~DBCALShower_factory(){};

	private:
		jerror_t evnt(jana::JEventLoop *loop, uint64_t eventnumber){

			// This is a trivial factory that simply implements the
			// IU tagged factory as the default. It is here so 
			// that the default can be changed easily by simply
			// changing the tag here or on the command line.
			vector<const DBCALShower*> showers;
			loop->Get(showers, "IU");
			for(unsigned int i=0; i<showers.size(); i++){
				_data.push_back(const_cast<DBCALShower*>(showers[i]));
			}
			SetFactoryFlag(NOT_OBJECT_OWNER);

			return NOERROR;
		}
};

#endif // _DBCALShower_factory_
