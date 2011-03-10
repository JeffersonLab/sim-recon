// $Id$
//
//    File: DEventSourceEVIO.h
// Created: Sat May  8 13:54:35 EDT 2010
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//

#ifndef _DEventSourceEVIO_
#define _DEventSourceEVIO_

#include <JANA/jerror.h>
#include <JANA/JEventSource.h>
#include <JANA/JEvent.h>

#include <HDGEOMETRY/DMagneticFieldMap.h>
#include <HDGEOMETRY/DGeometry.h>
#include <TRACKING/DTrackTimeBased.h>

#include <evioUtil.hxx>
#include <evioFileChannel.hxx>
#include <vector>

using namespace jana;
using namespace evio;
using namespace std;


class DEventSourceEVIO:public JEventSource{
	public:
		DEventSourceEVIO(const char* source_name);
		virtual ~DEventSourceEVIO();

		jerror_t GetEvent(JEvent &event);
		void FreeEvent(JEvent &event);
		jerror_t GetObjects(JEvent &event, JFactory_base *factory);
		
		template<class T> const vector<T>& GetVector(evioDOMNodeList* nodeList, string) const;
		
	protected:
	
		map<string, pair<int, int> > tagMap; // first=tag, second=num
		
		jerror_t Extract_DTrackTimeBased(evioDOMTree *evt,  JFactory<DTrackTimeBased> *factory);
	
	private:
		evioFileChannel *chan;
		
		DMagneticFieldMap *bfield;
		DGeometry *geom;

		pthread_mutex_t rt_mutex;
		map<evioDOMTree*, vector<DReferenceTrajectory*> > rt_by_event;
		list<DReferenceTrajectory*> rt_pool;
};

//------------------
// GetVector
//------------------
template<class T>
const vector<T>& DEventSourceEVIO::GetVector(evioDOMNodeList* nodeList, string name) const
{
	/// Attempt to get the pointer to a const vector of the specified type for the node
	/// matching the given name (as appearing in tagMap). If the pointer can't be obtained
	/// for any reason (e.g. name doesn't appear in tagMap or nodelist doesn't contain
	/// the tag/num correpsonding to name) then an exception is thrown.
	
	// Check that tagMap has name in it
	map<string, pair<int, int> >::const_iterator iter = tagMap.find(name);
	if(iter==tagMap.end())throw JException(string("can't find \"")+name+"\" in nodeList");
	
	// Loop over nodes and find one with the correct tag/num
	evioDOMNodeList::const_iterator iter2;
	for(iter2=nodeList->begin(); iter2!=nodeList->end(); iter2++) {
		evioDOMNode *node = *iter2;
		if(*node == iter->second){ // iter->second is pair<int,int> which holds tag/num
			return *(node->getVector<T>());
		}
	}
	
	// Didn't find bank with tag/num for this name. Throw exception
	throw JException(string("can't find with tag/num corresponding to \"")+name+"\" in nodeList");
	
	// Unused return value to avoid compiler warnings
	const vector<T> *tmp = NULL;
	return *tmp;
}

#endif // _DEventSourceEVIO_

