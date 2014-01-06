// $Id$
//
//    File: DMagneticFieldMapSpoiled.h
// Created: Wed Mar 25 04:04:16 EDT 2009
// Creator: davidl (on Darwin Amelia.local 9.6.0 i386)
//

#ifndef _DMagneticFieldMapSpoiled_
#define _DMagneticFieldMapSpoiled_

#include <JANA/jerror.h>

#include <HDGEOMETRY/DMagneticFieldMap.h>

#include <JANA/JApplication.h>
#include <JANA/JCalibration.h>
using namespace jana;

class DMagneticFieldMapCalibDB;

class DMagneticFieldMapSpoiled:public DMagneticFieldMap{
	public:
		DMagneticFieldMapSpoiled(JApplication *japp, unsigned int run_number=1, string namepath = "Magnets/Solenoid/solenoid_1500");
		DMagneticFieldMapSpoiled(JCalibration *jcalib, string namepath = "Magnets/Solenoid/solenoid_1500");
		virtual ~DMagneticFieldMapSpoiled();

		void GetField(const DVector3 &pos,DVector3 &Bout) const;
		void GetField(double x, double y, double z, double &Bx, double &By, double &Bz, int method=0) const;
		
		double GetBz(double x,double y, double z) const;
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
		
		void Init(void);
	
		bool initialized;
		DMagneticFieldMapCalibDB *bfield;
		double phi_amp;		///< Amplitude of phi spoiler (fraction of magnitude)
		double phi_omega;		///< Angular frequency of phi spoiler (radians/radian)
		double r_amp;			///< Amplitude of r spoiler (fraction of magnitude)
		double r_omega;			///< Angular frequency of r spoiler (radians/cm)
		double z_amp;			///< Amplitude of z spoiler (fraction of magnitude)
		double z_omega;			///< Angular frequency of z spoiler (radians/cm)
	
	private:

};

#endif // _DMagneticFieldMapSpoiled_

