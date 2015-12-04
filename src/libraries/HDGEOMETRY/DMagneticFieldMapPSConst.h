// $Id$
//
//    File: DMagneticFieldMapPSConst.h
//

#ifndef _DMagneticFieldMapPSConst_
#define _DMagneticFieldMapPSConst_

#include <JANA/jerror.h>

#include "DMagneticFieldMapPS.h"

#include <vector>
#include <string>
using std::vector;
using std::string;

#include <JANA/JApplication.h>
#include <JANA/JCalibration.h>
using namespace jana;

class DMagneticFieldMapPSConst:public DMagneticFieldMapPS{
	public:
		DMagneticFieldMapPSConst(JApplication *japp, string namepath = "Magnets/PairSpectrometer/PS_const_field");
		DMagneticFieldMapPSConst(JCalibration *jcalib, string namepath = "Magnets/PairSpectrometer/PS_const_field");
		DMagneticFieldMapPSConst(double Bx, double By, double Bz);		
		virtual ~DMagneticFieldMapPSConst();
		
		int GetValues(string namepath, int32_t runnumber=1, string context="");
		void GetField(const DVector3 &pos,DVector3 &Bout) const;
		void GetField(double x, double y, double z, double &Bx, double &By, double &Bz, int method=0) const;
		
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
		double Bx, By, Bz;

};

#endif // _DMagneticFieldMapPSConst_

