// $Id$
//
//    File: DMagneticFieldMapParameterized.h
// Created: Tue Oct 20 14:06:19 EDT 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.8.0 i386)
//

#ifndef _DMagneticFieldMapParameterized_
#define _DMagneticFieldMapParameterized_

#include <JANA/jerror.h>

#include <string>
using std::string;

#include <JANA/JApplication.h>
#include <JANA/JCalibration.h>

#include <DMatrix.h>
#include <HDGEOMETRY/DMagneticFieldMap.h>

class DMagneticFieldMapParameterized:public DMagneticFieldMap{
	public:
		DMagneticFieldMapParameterized(jana::JApplication *japp, string namepath = "Magnets/Solenoid/solenoid_1500_poisson_20090814_01_params");
		DMagneticFieldMapParameterized(jana::JCalibration *jcalib, string namepath = "Magnets/Solenoid/solenoid_1500_poisson_20090814_01_params");
		virtual ~DMagneticFieldMapParameterized();
		void Init(jana::JCalibration *jcalib, string namepath);

		void GetField(const DVector3 &pos,DVector3 &Bout) const;
		virtual void GetField(double x, double y, double z, double &Bx, double &By, double &Bz, int method=0) const;
		double GetBz(double x,double y,double z) const;
		virtual void GetFieldGradient(double x, double y, double z,
                                      double &dBxdx, double &dBxdy,
                                      double &dBxdz,
                                      double &dBydx, double &dBydy,
                                      double &dBydz,
                                      double &dBzdx, double &dBzdy,
				                          double &dBzdz) const;
		virtual	void GetFieldBicubic(double x,double y,double z, double &Bx,double &By,double &Bz) const;
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
						 double &dBzdz) const;
	protected:
		jana::JCalibration *jcalib;
		
		class Dsection{
			public:
				string namepath;
				string Bi;
				int section;
				double zmin;
				double zmax;
				double zmid;
				double znorm;
				double rmid;
				double rnorm;
				unsigned int order1;
				unsigned int order2;
				vector<vector<double> > pp;	// parameters from calibDB
				DMatrix Q;							// parameters matrix transformed by Chebyshev matrix
				vector<vector<double> > cc;	// parameters matrix transformed by Chebyshev matrix
				
				bool IsInRange(double &z) const {return z>=zmin && z<=zmax;}
				double Eval(double &r, double &z) const;
		};
		
		vector<Dsection> sections_Bx;
		vector<Dsection> sections_Bz;
};

#endif // _DMagneticFieldMapParameterized_

