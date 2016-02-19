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

	int pixel;	// pixel index
	float t;	// time

	void toStrings(vector<pair<string, string> >&items) const {
		AddString(items, "pixel", "%d", pixel);
		AddString(items, "t", "%1.3f", t);
	}
};

#endif /* DDIRCHIT_H_ */
