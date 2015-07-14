// $Id$
//
//    File: DTrackingResolutionGEANT.h
// Created: Mon Feb 25 15:06:17 EST 2008
// Creator: davidl (on Darwin fwing-dhcp13.jlab.org 8.11.1 i386)
//

#ifndef _DTrackingResolutionGEANT_
#define _DTrackingResolutionGEANT_

#include <TH2.h>
#include <TFile.h>

#include "DTrackingResolution.h"

class DTrackingResolutionGEANT:public DTrackingResolution{
	public:

		class TableInfo{
			public:
				TFile *file;
				TH2D* pt_res_hist;
				TH2D* theta_res_hist;
				TH2D* phi_res_hist;
				TH2D* efficiency_hist;
		};

		DTrackingResolutionGEANT();
		virtual ~DTrackingResolutionGEANT();
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DTrackingResolutionGEANT";}

		void ReadTableInfo(const char *fname, TableInfo &ti);
		
		// Accessor methods called through virtual method of DTrackingResolution
		void GetResolution(int geanttype, const TVector3 &mom, double &pt_res, double &theta_res, double &phi_res);
		double GetEfficiency(int geanttype, const TVector3 &mom);

		// Workhorse methods that actually do the work
		void GetResolution(TableInfo &ti, int geanttype, const TVector3 &mom, double &pt_res, double &theta_res, double &phi_res);
		double GetEfficiency(TableInfo &ti, int geanttype, const TVector3 &mom);

	private:
		TableInfo pion_info;
		TableInfo proton_info;
};

#endif // _DTrackingResolutionGEANT_

