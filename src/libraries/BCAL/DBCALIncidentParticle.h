// $Id: DBCALIncidentParticle.h 9441 2012-08-06 21:56:11Z davidl $
//
//    File: DBCALIncidentParticle.h
// Created: Thu Jun  9 10:14:35 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DBCALIncidentParticle_
#define _DBCALIncidentParticle_

#include "BCAL/DBCALGeometry.h"

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

class DBCALIncidentParticle:public JObject{

	/// This class holds data originating from the fADC250
	/// modules connected to the BCAL

	public:
		JOBJECT_PUBLIC(DBCALIncidentParticle);
		
		int ptype;
		float px;
		float py;
		float pz;
		float x;
		float y;
		float z;

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "ptype", "%d", ptype);
			AddString(items, "px[GeV]", "%7.5f", px);
			AddString(items, "py[GeV]", "%7.5f", py);
			AddString(items, "pz[GeV]", "%7.5f", pz);
			AddString(items, "x[GeV]", "%7.5f", x);
			AddString(items, "y[GeV]", "%7.5f", y);
			AddString(items, "z[GeV]", "%7.5f", z);
		}
};

#endif // _DBCALIncidentParticle_

