// $Id$
//
//    File: DMCTrigger.h
// Created: Tue Jun  7 10:15:05 EDT 2011
// Creator: davidl (on Darwin eleanor.jlab.org 10.7.0 i386)
//

#ifndef _DMCTrigger_
#define _DMCTrigger_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DMCTrigger:public jana::JObject{
	public:
		JOBJECT_PUBLIC(DMCTrigger);
		
		bool L1a_fired; // BCAL + 4FCAL >2 GeV && BCAL > 200 MeV && FCAL > 30 MeV
		bool L1b_fired; // BCAL + 4FCAL >2 GeV && BCAL > 30 MeV && FCAL > 30 MeV && NSC>0
		bool L1c_fired; // FCAL > 250MeV

        // energy sums with per-channel thresholds applied
		double Ebcal;
		double Efcal;		
        // energy sums without per-channel thresholds
		double Ebcal_all;   
		double Efcal_all;		
        // number of hits in fast sub-detectors
		unsigned int Nschits;
		unsigned int Ntofhits;

		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "L1a_fired", "%d", L1a_fired ? 1:0);
			AddString(items, "L1b_fired", "%d", L1b_fired ? 1:0);
			AddString(items, "Ebcal", "%5.3f", Ebcal);
			AddString(items, "Efcal", "%5.3f", Efcal);
			AddString(items, "Ebcal_all", "%5.3f", Ebcal_all);
			AddString(items, "Efcal_all", "%5.3f", Efcal_all);
			AddString(items, "Nschits", "%2d", Nschits);
			AddString(items, "Ntofhits", "%2d", Ntofhits);
		}
		
};

#endif // _DMCTrigger_

