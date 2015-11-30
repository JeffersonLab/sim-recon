// $Id$
//
//    File: DMagneticFieldMapPS2DMap.h

#ifndef _DMagneticFieldMapPS2DMap_
#define _DMagneticFieldMapPS2DMap_

#include <JANA/jerror.h>

#include <HDGEOMETRY/DMagneticFieldMapPS.h>

#include <vector>
#include <string>
using std::vector;
using std::string;

#include <JANA/JApplication.h>
#include <JANA/JCalibration.h>
using namespace jana;

class DGeometry;


class DMagneticFieldMapPS2DMap:public DMagneticFieldMapPS {
 public:
  DMagneticFieldMapPS2DMap(JApplication *japp, int32_t runnumber=1, string namepath = "Magnets/PairSpectrometer/PS_1.8T_20150513_test");
  DMagneticFieldMapPS2DMap(JCalibration *jcalib, string namepath = "Magnets/PairSpectrometer/PS_1.8T_20150513_test");
  virtual ~DMagneticFieldMapPS2DMap();
  
  int ReadMap(string namepath, int32_t runnumber=1, string context="");
  
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

  
  typedef struct{
    float x,y,z,Bx,By,Bz;
    double dBxdx, dBxdy, dBxdz;
    double dBydx, dBydy, dBydz;
    double dBzdx, dBzdy, dBzdz;
    double dBxdxdy,dBxdxdz,dBxdydz;
    double dBydxdy,dBydxdz,dBydydz;
    double dBzdxdy,dBzdxdz,dBzdydz;
  }DBfieldPoint_t;
  
 protected:
  
  JCalibration *jcalib;
  JResourceManager *jresman;
  DGeometry* geom;

  vector< vector< vector<DBfieldPoint_t> > > Btable;
  
  float xmin, xmax, ymin, ymax, zmin, zmax;
  int Nx, Ny, Nz;
  double dx, dy,dz;
  double one_over_dx,one_over_dz;

  double z_shift;
  
 private:
  void InterpolateField(double r,double z,double &Br,double &Bz,double &dBrdr,
			double &dBrdz,double &dBzdr,double &dBzdz) const;
};

#endif // _DMagneticFieldMapPS2DMap_

