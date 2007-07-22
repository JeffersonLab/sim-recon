// $Id$
//
//    File: DMagneticFieldMapCalibDB.h
// Created: Thu Jul 19 13:58:21 EDT 2007
// Creator: davidl (on Darwin fwing-dhcp61.jlab.org 8.10.1 i386)
//

#ifndef _DMagneticFieldMapCalibDB_
#define _DMagneticFieldMapCalibDB_

#include <JANA/jerror.h>

#include "DMagneticFieldMap.h"

#include <vector>
#include <string>
using std::vector;
using std::string;

class JApplication;
class JCalibration;

class DMagneticFieldMapCalibDB:public DMagneticFieldMap{
	public:
		DMagneticFieldMapCalibDB(JApplication *japp);
		DMagneticFieldMapCalibDB(JCalibration *jcalib);
		virtual ~DMagneticFieldMapCalibDB();
		
		int ReadMap(int runnumber=1, string context="");
		
		void GetField(double x, double y, double z, double &Bx, double &By, double &Bz, int method=0) const;

		typedef struct{
			float x,y,z,Bx,By,Bz;
			double dBxdx, dBxdy, dBxdz;
			double dBydx, dBydy, dBydz;
			double dBzdx, dBzdy, dBzdz;
		}DBfieldPoint_t;

	protected:
		
		JCalibration *jcalib;
		vector< vector< vector<DBfieldPoint_t> > > Btable;

		float xmin, xmax, ymin, ymax, zmin, zmax;
		int Nx, Ny, Nz;
		double dx, dy,dz;
};

#endif // _DMagneticFieldMapCalibDB_

