/*
 * DCereHit.h
 *
 *  Created on: Oct 10, 2012
 *      Author: yqiang
 */

#ifndef DCEREHIT_H_
#define DCEREHIT_H_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

class DCereHit: public JObject {

public:
	JOBJECT_PUBLIC (DCereHit);

	int sector; // sector hit by photons
	float pe;	// number of photon electrons
	float t;	// time

	void toStrings(vector<pair<string, string> >&items) const {
		AddString(items, "sector", "%d", sector);
		AddString(items, "pe", "%1.3f", pe);
		AddString(items, "t", "%1.3f", t);
	}
};

#endif /* DCEREHIT_H_ */
