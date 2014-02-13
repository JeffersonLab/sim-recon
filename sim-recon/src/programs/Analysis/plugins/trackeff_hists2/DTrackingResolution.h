// $Id$
//
//    File: DTrackingResolution.h
// Created: Mon Feb 25 15:06:17 EST 2008
// Creator: davidl (on Darwin fwing-dhcp13.jlab.org 8.11.1 i386)
//

#ifndef _DTrackingResolution_
#define _DTrackingResolution_

#include <TRandom3.h>
#include <TVector3.h>

class DTrackingResolution{
	public:
		DTrackingResolution();
		virtual ~DTrackingResolution();
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DTrackingResolution";}

		// Virtual methods that must be supplied by subclass
		// Momenta are in units of GeV/c and angular resolutions
		// are in units of milliradians.
		virtual void GetResolution(int geanttype, const TVector3 &mom, double &pt_res, double &theta_res, double &phi_res)=0;
		virtual double GetEfficiency(int geanttype, const TVector3 &mom)=0;

		// Methods implemented in this class
		bool Smear(int geanttype, TVector3 &mom);
		bool Efficiency(int geanttype, const TVector3 &mom);
	
	private:
		TRandom3 rnd;
};

#endif // _DTrackingResolution_

