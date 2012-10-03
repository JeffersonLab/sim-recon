/*
 * DCereTruth.h
 *
 *  Created on: Oct 3, 2012
 *      Author: yqiang
 */

#ifndef DCERETRUTH_H_
#define DCERETRUTH_H_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

class DCereTruth: public JObject {

public:
	JOBJECT_PUBLIC (DCereTruth);

	int track;		// track index
	int primary;	// 0: secondary, 1: primary
	int ptype;		// particle type
	float x, y, z;	// true point of intersection
	float px, py, pz;	// momentum
	float t;	// time

	void toStrings(vector<pair<string, string> >&items) const {
		AddString(items, "track", "%d", track);
		AddString(items, "primary", "%d", primary);
		AddString(items, "ptype", "%d", ptype);
		AddString(items, "x", "%1.3f", x);
		AddString(items, "y", "%1.3f", y);
		AddString(items, "z", "%1.3f", z);
		AddString(items, "t", "%1.3f", t);
		AddString(items, "px", "%1.3f", px);
		AddString(items, "py", "%1.3f", py);
		AddString(items, "pz", "%1.3f", pz);
	}
};


#endif /* DCERETRUTH_H_ */
