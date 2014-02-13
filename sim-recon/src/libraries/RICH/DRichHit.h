/*
 * DRichHit.h
 *
 *  Created on: Oct 11, 2012
 *      Author: yqiang
 */

#ifndef DRICHHIT_H_
#define DRICHHIT_H_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

class DRichHit: public JObject {

public:
	JOBJECT_PUBLIC (DRichHit);

	float x, y, z;	// true point of intersection
	float t;	// time

	void toStrings(vector<pair<string, string> >&items) const {
		AddString(items, "x", "%1.3f", x);
		AddString(items, "y", "%1.3f", y);
		AddString(items, "z", "%1.3f", z);
		AddString(items, "t", "%1.3f", t);
	}
};

#endif /* DRICHHIT_H_ */
