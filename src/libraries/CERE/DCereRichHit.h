/*
 * DCereRichHit.h
 *
 *  Created on: Oct 3, 2012
 *      Author: yqiang
 */

#ifndef DCERERICHHIT_H_
#define DCERERICHHIT_H_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

class DCereRichHit: public JObject {

public:
	JOBJECT_PUBLIC (DCereRichHit);

	float x, y, z;	// true point of intersection
	float t;	// time

	void toStrings(vector<pair<string, string> >&items) const {
		AddString(items, "x", "%1.3f", x);
		AddString(items, "y", "%1.3f", y);
		AddString(items, "z", "%1.3f", z);
		AddString(items, "t", "%1.3f", t);
	}
};


#endif /* DCERERICHHIT_H_ */
