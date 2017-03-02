// $Id$
//
//    File: DBeamCurrent.h
// Created: Tue Feb 21 04:25:03 EST 2017
// Creator: davidl (on Linux gluon48.jlab.org 2.6.32-431.20.3.el6.x86_64 x86_64)
//

#ifndef _DBeamCurrent_
#define _DBeamCurrent_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DBeamCurrent:public jana::JObject{
	public:
		JOBJECT_PUBLIC(DBeamCurrent);
		
		// This gives the approximate electron beam current
		// for the current event. The values are obtained 
		// from the CCDB which was populated from values
		// in the EPICS archive. The values were obtained
		// from EPICS using the mySampler program so they
		// are interpolations between values to estimate the
		// current at 1 sec. intervals throughout the run.
		// These values were further reduced by looking
		// for boundaries where the current changes by more
		// than 3nA. Thus, the beam current cannot be considered
		// accurate to better than 3nA.
		//
		// The t_prev and t_next values give the time in
		// seconds of the previous and next beam trips.
		// The method IsFiducial() can be used to easily check 
		// if this event is in a region far from any beam
		// trips and that the beam current was greater then
		// some minumum (default is 5nA).
	
		double Ibeam;  // electron beam current in nA
		double t;      // time relative to start of run in seconds
		double t_prev; // time since last trip in sec
		double t_next; // time to next trip in sec
		bool   is_fiducial; // true if this event is in fiducial region of time
		
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "Ibeam", "%5.1f", Ibeam);
			AddString(items, "t",      "%f", t);
			AddString(items, "t_prev", "%f", t_prev);
			AddString(items, "t_next", "%f", t_next);
			AddString(items, "is_fiducial", "%d", is_fiducial ? 1:0);
		}
		
};

#endif // _DBeamCurrent_

