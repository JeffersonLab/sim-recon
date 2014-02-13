// $Id$
//
//    File: DMagneticFieldMapNoField.h
// Created: Fri Nov  7 04:01:28 EST 2008
// Creator: davidl (on Darwin Amelia.local 8.11.1 i386)
//

#ifndef _DMagneticFieldMapNoField_
#define _DMagneticFieldMapNoField_

#include <JANA/jerror.h>

#include "DMagneticFieldMap.h"

#include <vector>
#include <string>
using std::vector;
using std::string;

#include <JANA/JApplication.h>
#include <JANA/JCalibration.h>
using namespace jana;

class DMagneticFieldMapNoField:public DMagneticFieldMap{
	public:
		DMagneticFieldMapNoField(JApplication *japp, string namepath = "");
		DMagneticFieldMapNoField(JCalibration *jcalib, string namepath = "");
		virtual ~DMagneticFieldMapNoField();
		
 void GetField(const DVector3 &pos,DVector3 &Bout) const{
    Bout.SetXYZ(0.,0.,0.);
  };
  void GetField(double x, double y, double z, double &Bx, double &By, double &Bz, int method=0) const{
    Bx=0.;
    By=0.;
    Bz=0.;
  };
  double GetBz(double x,double y, double z) const {return 0.;};
		
  void GetFieldGradient(double x, double y, double z,
			double &dBxdx, double &dBxdy,
			double &dBxdz,
			double &dBydx, double &dBydy,
			double &dBydz,		
				      double &dBzdx, double &dBzdy,
			double &dBzdz) const{
    dBxdx = 0.0;
    dBxdy = 0.0;
    dBxdz = 0.0;
    dBydx = 0.0;
    dBydy = 0.0;
    dBydz = 0.0;
    dBzdx = 0.0;
    dBzdy = 0.0;
    dBzdz = 0.0;
  };
  
  void GetFieldAndGradient(double x,double y,double z,
			   double &Bx,double &By,
			   double &Bz,
			   double &dBxdx, double &dBxdy,
			   double &dBxdz,
			   double &dBydx, double &dBydy,
			   double &dBydz,
			   double &dBzdx, double &dBzdy,
			   double &dBzdz) const{
    Bx=0.;
    By=0.;
    Bz=0.;
    dBxdx = 0.0;
    dBxdy = 0.0;
    dBxdz = 0.0;
    dBydx = 0.0;
    dBydy = 0.0;
    dBydz = 0.0;
    dBzdx = 0.0;
    dBzdy = 0.0;
    dBzdz = 0.0;
    
  };
  void GetFieldBicubic(double x,double y,double z,
		       double &Bx,double &By,double &Bz) 
  const{
    Bx=0.;
    By=0.;
    Bz=0.;
  };
};

#endif // _DMagneticFieldMapNoField_

