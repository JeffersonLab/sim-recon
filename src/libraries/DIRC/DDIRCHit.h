/*
 * DDIRCHit.h
 *
 *  Created on: Oct 11, 2012
 *      Author: yqiang
 */

#ifndef DDIRCHIT_H_
#define DDIRCHIT_H_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

class DDIRCHit: public JObject {

public:
	JOBJECT_PUBLIC (DDIRCHit);

	float x, y, z;	// true point of intersection
	float t;	// time
	float E;	// energy

	void toStrings(vector<pair<string, string> >&items) const {
		AddString(items, "x", "%1.3f", x);
		AddString(items, "y", "%1.3f", y);
		AddString(items, "z", "%1.3f", z);
		AddString(items, "t", "%1.3f", t);
	}
};

#endif /* DDIRCHIT_H_ */
