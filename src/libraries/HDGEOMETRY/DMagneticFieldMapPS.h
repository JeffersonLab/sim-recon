// $Id$
//
//    File: DMagneticFieldMapPS.h
//

#ifndef _DMagneticFieldMapPS_
#define _DMagneticFieldMapPS_

#include <JANA/jerror.h>
#include <DVector3.h>

class DMagneticFieldMapPS{
	public:
	
		DMagneticFieldMapPS(){}
		virtual ~DMagneticFieldMapPS(){}

		virtual void GetField(const DVector3 &pos,DVector3 &Bout) const=0;
		virtual void GetField(double x, double y, double z, double &Bx, double &By, double &Bz, int method=0) const = 0;

		virtual void GetFieldGradient(double x, double y, double z,
                                      double &dBxdx, double &dBxdy,
                                      double &dBxdz,
                                      double &dBydx, double &dBydy,
                                      double &dBydz,
                                      double &dBzdx, double &dBzdy,
				      double &dBzdz) const = 0;
		virtual	void GetFieldBicubic(double x,double y,double z,
				     double &Bx,double &By,double &Bz) const=0;

		virtual void GetFieldAndGradient(double x,double y,double z,
						 double &Bx,double &By,
						 double &Bz,
					         double &dBxdx, 
						 double &dBxdy,
						 double &dBxdz,
						 double &dBydx, 
						 double &dBydy,
						 double &dBydz,
						 double &dBzdx, 
						 double &dBzdy,
						 double &dBzdz) const = 0;


};

#endif // _DMagneticFieldMapPS_

