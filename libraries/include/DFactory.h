// $Id$

/// Hall-D Data Factory
///
/// All data (except that read from the data source) must
/// be derived from the data that was read from the source.
/// One DFactory object should exist for each type of data
/// that is "generated". For example: Clusters in the FCAL
/// are generated from Hits in the FCAL. Asking the factory
/// for a list of FCAL clusters will make the factory ask
/// for a list of FCAL Hits from the FCAL Hits factory which
/// may, in turn, ask for a list of raw data values from yet
/// another factory).
/// A list of DFactories is kept in the DEvent object so
/// all data may be accessed through it.

#ifndef _DFACTORY_H_
#define _DFACTORY_H_

#include "DEventProcessor.h"
#include "DContainer.h"

class DEvent;

class DFactory:public DEventProcessor{
	public:
		DFactory(DEvent *my_devent, char *my_name, int rowsize);
		~DFactory();
		DContainer* Get();  ///< Get a list of the data produced by this factory
		inline void ResetNrows(void){_data->ResetNrows();}
		virtual derror_t Print(void){cout<<" <Print method undefined for "<<name<<">"<<endl;}

		char* name;

	protected:
		DEvent *event;
};

#endif // _DFACTORY_H_
