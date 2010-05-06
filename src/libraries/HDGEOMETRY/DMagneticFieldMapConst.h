// $Id$
//
//    File: DMagneticFieldMapConst.h
// Created: Fri Nov  7 04:01:28 EST 2008
// Creator: davidl (on Darwin Amelia.local 8.11.1 i386)
//

#ifndef _DMagneticFieldMapConst_
#define _DMagneticFieldMapConst_

#include <JANA/jerror.h>

#include "DMagneticFieldMap.h"

#include <vector>
#include <string>
using std::vector;
using std::string;

#include <JANA/JApplication.h>
#include <JANA/JCalibration.h>
using namespace jana;

class DMagneticFieldMapConst:public DMagneticFieldMap{
	public:
		DMagneticFieldMapConst(JApplication *japp, string namepath = "Magnets/Solenoid/solenoid_const");
		DMagneticFieldMapConst(JCalibration *jcalib, string namepath = "Magnets/Solenoid/solenoid_const");
		virtual ~DMagneticFieldMapConst();
		
		int GetValues(string namepath, int runnumber=1, string context="");
		
		void GetField(double x, double y, double z, double &Bx, double &By, double &Bz, int method=0) const;
		double GetBz(double x,double y, double z) const {return Bz;};
		
		void GetFieldGradient(double x, double y, double z,
				      double &dBxdx, double &dBxdy,
				      double &dBxdz,
				      double &dBydx, double &dBydy,
				      double &dBydz,		
				      double &dBzdx, double &dBzdy,
				      double &dBzdz) const;
		
		void GetFieldBicubic(double x,double y,double z,
				     double &Bx,double &By,double &Bz) const;
		void GetFieldAndGradient(double x,double y,double z,
					 double &Bx,double &By,
					 double &Bz,
					 double &dBxdx, double &dBxdy,
					 double &dBxdz,
					 double &dBydx, double &dBydy,
					 double &dBydz,
					 double &dBzdx, double &dBzdy,
					 double &dBzdz) const;

	protected:
		
		JCalibration *jcalib;
		double Br, Bphi, Bz;

};

#endif // _DMagneticFieldMapConst_

