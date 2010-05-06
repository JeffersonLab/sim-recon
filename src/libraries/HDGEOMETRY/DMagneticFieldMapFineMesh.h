// $Id$
//
//    File: DMagneticFieldMapFineMesh.h

#ifndef _DMagneticFieldMapFineMesh_
#define _DMagneticFieldMapFineMesh_

#include <JANA/jerror.h>

#include <HDGEOMETRY/DMagneticFieldMap.h>

#include <vector>
#include <string>
using std::vector;
using std::string;

#include <JANA/JApplication.h>
#include <JANA/JCalibration.h>
using namespace jana;

class DMagneticFieldMapFineMesh:public DMagneticFieldMap{
 public:
  DMagneticFieldMapFineMesh(JApplication *japp, string namepath = "Magnets/Solenoid/solenoid_1500_poisson_20090814_01");
  DMagneticFieldMapFineMesh(JCalibration *jcalib, string namepath = "Magnets/Solenoid/solenoid_1500_poisson_20090814_01");
  virtual ~DMagneticFieldMapFineMesh();
  
  int ReadMap(string namepath, int runnumber=1, string context="");
  
  void GetField(double x, double y, double z, double &Bx, double &By, double &Bz, int method=0) const;
	double GetBz(double x, double y, double z) const; 
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
  void GetFineMeshMap(void);
  void WriteEvioFile(string evioFileName);	
  void ReadEvioFile(string evioFileName);
  void GenerateFineMesh(void);
  
  typedef struct{
    float x,y,z,Bx,By,Bz;
    double dBxdx, dBxdy, dBxdz;
    double dBydx, dBydy, dBydz;
    double dBzdx, dBzdy, dBzdz;
    double dBxdxdy,dBxdxdz,dBxdydz;
    double dBydxdy,dBydxdz,dBydydz;
    double dBzdxdy,dBzdxdz,dBzdydz;
  }DBfieldPoint_t;
  
  typedef struct{
    double Br,Bz;
    double dBrdr,dBrdz,dBzdr,dBzdz;
  }DBfieldCylindrical_t;
  
 protected:
  
  JCalibration *jcalib;
  vector< vector< vector<DBfieldPoint_t> > > Btable;
  
  float xmin, xmax, ymin, ymax, zmin, zmax;
  int Nx, Ny, Nz;
  double dx, dy,dz;
  
  vector<vector<DBfieldCylindrical_t> >mBfine;
  float zminFine,rminFine,zmaxFine,rmaxFine,drFine,dzFine;
  unsigned int NrFine,NzFine;
  
 private:
  void InterpolateField(double r,double z,double &Br,double &Bz,double &dBrdr,
			double &dBrdz,double &dBzdr,double &dBzdz) const;
};

#endif // _DMagneticFieldMapFineMesh_

