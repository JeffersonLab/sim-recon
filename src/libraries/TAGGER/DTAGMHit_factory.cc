// $Id$
//
//    File: DTAGMHit_factory.cc
// Created: Sat Aug  2 12:23:43 EDT 2014
// Creator: jonesrt (on Linux gluex.phys.uconn.edu)
//

// aebarnes moved original factory to DTAGMHit_factory_Calib.cc on August 11, 2016.

#include <iostream>
#include <iomanip>
#include <TAGGER/DTAGMHit.h>
using namespace std;

// Sort by column
bool SortByCol(const DTAGMHit* a, const DTAGMHit* b)
{
	return a->column < b->column;
}

#include "DTAGMHit_factory.h"

using namespace jana;

//------------------
// init
//------------------
jerror_t DTAGMHit_factory::init(void)
{
    // Set default configuration parameters
    DELTA_T_CLUSTER_MAX = 5;
    gPARMS->SetDefaultParameter("TAGMHit:DELTA_T_CLUSTER_MAX",DELTA_T_CLUSTER_MAX,
                                "Maximum time difference in ns between hits in adjacent"
                                " columns to be merged");

    MERGE_HITS = true;
    gPARMS->SetDefaultParameter("TAGMHit:MERGE_HITS",MERGE_HITS,
                                "Merge neighboring hits when true");

    // Setting this flag makes it so that JANA does not delete the objects in _data.
    // This factory will manage this memory.
    SetFactoryFlag(NOT_OBJECT_OWNER);

    return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTAGMHit_factory::brun(jana::JEventLoop *eventLoop, int32_t runnumber)
{
    return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTAGMHit_factory::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	// Clear _data vector
	Reset_Data();

	// Get (calibrated) TAGM hits
	vector<const DTAGMHit*> hits;
	loop->Get(hits, "Calib");

	sort(hits.begin(), hits.end(), SortByCol);

	set<uint32_t> locHitIndexUsedSoFar;

	for(uint32_t i = 0; i < hits.size(); ++i)
	{
		const DTAGMHit *hit_i = hits[i];

		if (hit_i->row > 0 || !MERGE_HITS)
		{
			_data.push_back(const_cast<DTAGMHit*>(hit_i));
			continue;
		}
		if (!hit_i->has_fADC) continue;

		// check if column has been paired
		if (locHitIndexUsedSoFar.find(i) != locHitIndexUsedSoFar.end()) continue;

		for (uint32_t j = i+1; j < hits.size(); ++j)
		{
			const DTAGMHit *hit_j = hits[j];

			if (!hit_j->has_fADC) continue;
			if (hit_j->row > 0) continue;

			int colDiff = hit_i->column - hit_j->column;
			double deltaT = hit_i->t - hit_j->t;

			if (abs(colDiff) == 1 && fabs(deltaT) <= DELTA_T_CLUSTER_MAX)
			{
				locHitIndexUsedSoFar.insert(j);
				break;
			}
		}
		_data.push_back(const_cast<DTAGMHit*>(hit_i));
	}


    return NOERROR;
}

//------------------
// Reset_Data()
//------------------
void DTAGMHit_factory::Reset_Data(void)
{
	// Clear _data vector
	_data.clear();
}

//------------------
// erun
//------------------
jerror_t DTAGMHit_factory::erun(void)
{
    return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTAGMHit_factory::fini(void)
{
    return NOERROR;
}
