// $Id$
//
//    File: DMagneticFieldMapCalibDB.cc
// Created: Thu Jul 19 13:58:21 EDT 2007
// Creator: davidl (on Darwin fwing-dhcp61.jlab.org 8.10.1 i386)
//

#include <cmath>
using namespace std;

#include "DMagneticFieldMapCalibDB.h"

//---------------------------------
// DMagneticFieldMapCalibDB    (Constructor)
//---------------------------------
DMagneticFieldMapCalibDB::DMagneticFieldMapCalibDB(JApplication *japp, string namepath)
{
	int runnumber = 1;
	jcalib = japp->GetJCalibration(runnumber);

	JParameterManager *jparms = japp->GetJParameterManager();
	jparms->SetDefaultParameter("BFIELD_MAP", namepath);
	
	int Npoints = ReadMap(namepath, runnumber); 
	if(Npoints==0){
		_DBG_<<"Error getting JCalibration object for magnetic field!"<<endl;
		japp->Quit();
	}
}

//---------------------------------
// DMagneticFieldMapCalibDB    (Constructor)
//---------------------------------
DMagneticFieldMapCalibDB::DMagneticFieldMapCalibDB(JCalibration *jcalib, string namepath)
{
	this->jcalib = jcalib;
	if(ReadMap(namepath)==0){
		_DBG_<<"Error getting JCalibration object for magnetic field!"<<endl;
		exit(1);
	} 
}

//---------------------------------
// ~DMagneticFieldMapCalibDB    (Destructor)
//---------------------------------
DMagneticFieldMapCalibDB::~DMagneticFieldMapCalibDB()
{

}

//---------------------------------
// ReadMap
//---------------------------------
int DMagneticFieldMapCalibDB::ReadMap(string namepath, int runnumber, string context)
{
  /// Read the magnetic field map in from the calibration database.
  /// This will read in the map and figure out the number of grid
  /// points in each direction (x,y, and z) and the range in each.
  /// The gradiant of the field is calculated for all but the most
  /// exterior points and saved to use in later calls to GetField(...).

  // Read in map from calibration database. This should really be
  // built into a run-dependent geometry framework, but for now
  // we do it this way. 
  if(!jcalib)return 0;
  
  jout<<"Reading Magnetic field map from "<<namepath<<" ..."<<endl;
  vector< vector<float> > Bmap;
  jcalib->Get(namepath, Bmap);
  jout<<Bmap.size()<<" entries found (";
  
  // The map should be on a grid with equal spacing in x,y, and z.
  // Here we want to determine the number of points in each of these
  // dimensions and the range. Note that all of or maps right now are
  // 2-dimensional with extent only in x and z.
  // The easiest way to do this is to use a map<float, int> to make a
  // histogram of the entries by using the key to hold the extent
  // so that the number of entries will be equal to the number of
  // different values.
  map<float, int> xvals;
  map<float, int> yvals;
  map<float, int> zvals;
  xmin = ymin = zmin = 1.0E6;
  xmax = ymax = zmax = -1.0E6;
  for(unsigned int i=0; i<Bmap.size(); i++){
    vector<float> &a = Bmap[i];
    float &x = a[0];
    float &y = a[1];
    float &z = a[2];
    
    // Convert from inches to cm and shift in z by 26 inches to align
    // the origin with the GEANT defined one.
    x *= 2.54;
    y *= 2.54;
    z = (z-0.0)*2.54; // Removed 26" offset 6/22/2009 DL
    
    xvals[x] = 1;
    yvals[y] = 1;
    zvals[z] = 1;
    if(x<xmin)xmin=x;
    if(y<ymin)ymin=y;
    if(z<zmin)zmin=z;
    if(x>xmax)xmax=x;
    if(y>ymax)ymax=y;
    if(z>zmax)zmax=z;
  }
  Nx = xvals.size();
  Ny = yvals.size();
  Nz = zvals.size();
  jout<<" Nx="<<Nx;
  jout<<" Ny="<<Ny;
  jout<<" Nz="<<Nz;
  jout<<" )  at 0x"<<hex<<(unsigned long)this<<dec<<endl;
  
  // Create 3D vector so we can index the values by [x][y][z]
  vector<DBfieldPoint_t> zvec(Nz);
  vector< vector<DBfieldPoint_t> > yvec;
  for(int i=0; i<Ny; i++)yvec.push_back(zvec);
  for(int i=0; i<Nx; i++)Btable.push_back(yvec);
	
  // Distance between map points for r and z
  dx = (xmax-xmin)/(double)(Nx-1);
  dy = (ymax-ymin)/(double)(Ny-1);
  dz = (zmax-zmin)/(double)(Nz-1);
  
  // Copy values into Btable
  for(unsigned int i=0; i<Bmap.size(); i++){
    vector<float> &a = Bmap[i];
    int xindex = (int)floor((a[0]-xmin+dx/2.0)/dx); // the +dx/2.0 guarantees against round-off errors
    int yindex = (int)(Ny<2 ? 0:floor((a[1]-ymin+dy/2.0)/dy));
    int zindex = (int)floor((a[2]-zmin+dz/2.0)/dz);
    DBfieldPoint_t *b = &Btable[xindex][yindex][zindex];
    b->x = a[0];
    b->y = a[1];
    b->z = a[2];
    b->Bx = a[3];
    b->By = a[4];
    b->Bz = a[5];
  }
  
  // Calculate gradient at every point in the map.
  // The derivatives are calculated wrt to the *index* of the
  // field map, not the physical unit. This is because
  // the fractions used in interpolation are fractions
  // between indices.
  for(int index_x=0; index_x<Nx; index_x++){
    int index_x0 = index_x - (index_x>0 ? 1:0);
    int index_x1 = index_x + (index_x<Nx-1 ? 1:0);
    for(int index_z=1; index_z<Nz-1; index_z++){
      int index_z0 = index_z - (index_z>0 ? 1:0);
      int index_z1 = index_z + (index_z<Nz-1 ? 1:0);
      
      // Right now, all of our maps are 2-D with no y extent
      // so we comment out the following for-loop and hardwire
      // the index_y values to zero
      //for(int index_y=1; index_y<Ny-1; index_y++){
      //	int index_y0 = index_y-1;
      //	int index_y1 = index_y+1;
      if(Ny>1){
	_DBG_<<"Field map appears to be 3 dimensional. Code is currently"<<endl;
	_DBG_<<"unable to handle this. Exiting..."<<endl;
	exit(-1);
      }else{
	int index_y=0;
	/*
	int index_y0, index_y1;
	index_y=index_y0=index_y1=0;
	*/
	double d_index_x=double(index_x1-index_x0);
	double d_index_z=double(index_z1-index_z0); 

	DBfieldPoint_t *Bx0 = &Btable[index_x0][index_y][index_z];
	DBfieldPoint_t *Bx1 = &Btable[index_x1][index_y][index_z];
	//DBfieldPoint_t *By0 = &Btable[index_x][index_y0][index_z];
	//DBfieldPoint_t *By1 = &Btable[index_x][index_y1][index_z];
	DBfieldPoint_t *Bz0 = &Btable[index_x][index_y][index_z0];
	DBfieldPoint_t *Bz1 = &Btable[index_x][index_y][index_z1];
	
	DBfieldPoint_t *g = &Btable[index_x][index_y][index_z];
	g->dBxdx = (Bx1->Bx - Bx0->Bx)/d_index_x;
	//g->dBxdy = (By1->Bx - By0->Bx)/(double)(index_y1-index_y0);
	g->dBxdz = (Bz1->Bx - Bz0->Bx)/d_index_z;
	
	g->dBydx = (Bx1->By - Bx0->By)/d_index_x;
	//g->dBydy = (By1->By - By0->By)/(double)(index_y1-index_y0);
	g->dBydz = (Bz1->By - Bz0->By)/d_index_z;
       	
	g->dBzdx = (Bx1->Bz - Bx0->Bz)/d_index_x;
	//g->dBzdy = (By1->Bz - By0->Bz)/(double)(index_y1-index_y0);
	g->dBzdz = (Bz1->Bz - Bz0->Bz)/d_index_z;

	DBfieldPoint_t *B11 = &Btable[index_x1][index_y][index_z1];
	DBfieldPoint_t *B01 = &Btable[index_x0][index_y][index_z1];	
	DBfieldPoint_t *B10 = &Btable[index_x1][index_y][index_z0];
	DBfieldPoint_t *B00 = &Btable[index_x0][index_y][index_z0];
	
	g->dBxdxdz=(B11->Bx - B01->Bx - B10->Bx + B00->Bx)/d_index_x/d_index_z;
	g->dBzdxdz=(B11->Bz - B01->Bz - B10->Bz + B00->Bz)/d_index_x/d_index_z;
      }
    }
  }
  
  return Bmap.size();
}

// Use bicubic interpolation to find the field at the point (x,y).  
//See Numerical Recipes in C (2nd ed.), pp.125-127.
void DMagneticFieldMapCalibDB::GetFieldBicubic(double x,double y,double z,
						   double &Bx_,double &By_,
						   double &Bz_) const{
  // table of weight factors for bicubic interpolation
  static const int wt[16][16]=
    { {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
      {-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0},
      {2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0},
      {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
      {0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1},
      {0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1},
      {-3,3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0},
      {9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2},
      {-6,6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2},
      {2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0},
      {-6, 6,-6,6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1},
      {4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1}};
  double temp[16];
  double cl[16],coeff[4][4];

  // Get closest indices for this point
  double r = sqrt(x*x + y*y);
  int index_x = (int)floor((r-xmin)/dx + 0.5);
  //if(index_x<0 || index_x>=Nx)return;
  if (index_x>=Nx) return;
  else if (index_x<0) index_x=0;
  int index_z = (int)floor((z-zmin)/dz + 0.5);	
  if(index_z<0 || index_z>=Nz)return;
  
  int index_y = 0;
  
  int i, j, k, m=0;
  int index_x1 = index_x + (index_x<Nx-1 ? 1:0);
  int index_z1 = index_z + (index_z<Nz-1 ? 1:0);

  // Pointers to magnetic field structure
  const DBfieldPoint_t *B00 = &Btable[index_x][index_y][index_z];
  const DBfieldPoint_t *B01 = &Btable[index_x][index_y][index_z1];
  const DBfieldPoint_t *B11 = &Btable[index_x1][index_y][index_z1]; 
  const DBfieldPoint_t *B10 = &Btable[index_x1][index_y][index_z];
  
  // First compute the interpolation for Br
  temp[0]=B00->Bx;
  temp[1]=B01->Bx;
  temp[2]=B11->Bx;
  temp[3]=B10->Bx;

  temp[8]=B00->dBxdx;
  temp[9]=B01->dBxdx;
  temp[10]=B11->dBxdx;
  temp[11]=B10->dBxdx;

  temp[4]=B00->dBxdz;
  temp[5]=B01->dBxdz;
  temp[6]=B11->dBxdz;
  temp[7]=B10->dBxdz;
   
  temp[12]=B00->dBxdxdz;
  temp[13]=B01->dBxdxdz;
  temp[14]=B11->dBxdxdz;
  temp[15]=B10->dBxdxdz;

  for (i=0;i<16;i++){
    double tmp2=0.0;
    for (k=0;k<16;k++) tmp2+=wt[i][k]*temp[k];
    cl[i]=tmp2;
  }
  for (i=0;i<4;i++)
    for (j=0;j<4;j++) coeff[i][j]=cl[m++];

  double t=(z - B00->z)/dz;
  double u=(r - B00->x)/dx;
  double Br_=0.;
  for (i=3;i>=0;i--){
    Br_=t*Br_+((coeff[i][3]*u+coeff[i][2])*u+coeff[i][1])*u+coeff[i][0];
  }
  
  // Next compute the interpolation for Bz
  temp[0]=B00->Bz;
  temp[1]=B01->Bz;
  temp[2]=B11->Bz;
  temp[3]=B10->Bz;

  temp[8]=B00->dBzdx;
  temp[9]=B01->dBzdx;
  temp[10]=B11->dBzdx;
  temp[11]=B10->dBzdx;

  temp[4]=B00->dBzdz;
  temp[5]=B01->dBzdz;
  temp[6]=B11->dBzdz;
  temp[7]=B10->dBzdz;
   
  temp[12]=B00->dBzdxdz;
  temp[13]=B01->dBzdxdz;
  temp[14]=B11->dBzdxdz;
  temp[15]=B10->dBzdxdz;

  for (i=0;i<16;i++){
    double tmp2=0.0;
    for (k=0;k<16;k++) tmp2+=wt[i][k]*temp[k];
    cl[i]=tmp2;
  }
  m=0;
  for (i=0;i<4;i++)
    for (j=0;j<4;j++) coeff[i][j]=cl[m++];

  Bz_=0.;
  for (i=3;i>=0;i--){
    Bz_=t*Bz_+((coeff[i][3]*u+coeff[i][2])*u+coeff[i][1])*u+coeff[i][0];
  }

  // Convert r back to x,y components
  double cos_theta = x/r;
  double sin_theta = y/r;
  if(r==0.0){
    cos_theta=1.0;
    sin_theta=0.0;
  }
  // Rotate back into phi direction
  Bx_=Br_*cos_theta;
  By_=Br_*sin_theta;

}


// Use bicubic interpolation to find the field and field gradient at the point 
// (x,y).  See Numerical Recipes in C (2nd ed.), pp.125-127.
void DMagneticFieldMapCalibDB::GetFieldAndGradient(double x,double y,double z,
						   double &Bx_,double &By_,
						   double &Bz_,
						   double &dBxdx_, 
						   double &dBxdy_,
						   double &dBxdz_,
						   double &dBydx_, 
						   double &dBydy_,
						   double &dBydz_,
						   double &dBzdx_, 
						   double &dBzdy_,
						   double &dBzdz_) const{
  // table of weight factors for bicubic interpolation
  static const int wt[16][16]=
    { {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
      {-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0},
      {2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0},
      {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
      {0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1},
      {0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1},
      {-3,3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0},
      {9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2},
      {-6,6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2},
      {2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0},
      {-6, 6,-6,6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1},
      {4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1}};
  double temp[16];
  double cl[16],coeff[4][4];

  // Get closest indices for this point
  double r = sqrt(x*x + y*y);
  int index_x = (int)floor((r-xmin)/dx + 0.5);
  //if(index_x<0 || index_x>=Nx)return;
  if (index_x>=Nx) return;
  else if (index_x<0) index_x=0;
  int index_z = (int)floor((z-zmin)/dz + 0.5);	
  if(index_z<0 || index_z>=Nz)return;
  
  int index_y = 0;
  
  int i, j, k, m=0;
  int index_x1 = index_x + (index_x<Nx-1 ? 1:0);
  int index_z1 = index_z + (index_z<Nz-1 ? 1:0);

  // Pointers to magnetic field structure
  const DBfieldPoint_t *B00 = &Btable[index_x][index_y][index_z];
  const DBfieldPoint_t *B01 = &Btable[index_x][index_y][index_z1];
  const DBfieldPoint_t *B11 = &Btable[index_x1][index_y][index_z1]; 
  const DBfieldPoint_t *B10 = &Btable[index_x1][index_y][index_z];
  
  // First compute the interpolation for Br
  temp[0]=B00->Bx;
  temp[1]=B01->Bx;
  temp[2]=B11->Bx;
  temp[3]=B10->Bx;

  temp[8]=B00->dBxdx;
  temp[9]=B01->dBxdx;
  temp[10]=B11->dBxdx;
  temp[11]=B10->dBxdx;

  temp[4]=B00->dBxdz;
  temp[5]=B01->dBxdz;
  temp[6]=B11->dBxdz;
  temp[7]=B10->dBxdz;
   
  temp[12]=B00->dBxdxdz;
  temp[13]=B01->dBxdxdz;
  temp[14]=B11->dBxdxdz;
  temp[15]=B10->dBxdxdz;

  for (i=0;i<16;i++){
    double tmp2=0.0;
    for (k=0;k<16;k++) tmp2+=wt[i][k]*temp[k];
    cl[i]=tmp2;
  }
  for (i=0;i<4;i++)
    for (j=0;j<4;j++) coeff[i][j]=cl[m++];

  double t=(z - B00->z)/dz;
  double u=(r - B00->x)/dx;
  double Br_=0.,dBrdx_=0.,dBrdz_=0.;
  for (i=3;i>=0;i--){
    Br_=t*Br_+((coeff[i][3]*u+coeff[i][2])*u+coeff[i][1])*u+coeff[i][0];
    dBrdx_=t*dBrdx_+(3.*coeff[i][3]*u+2.*coeff[i][2])*u+coeff[i][1];
    dBrdz_=u*dBrdz_+(3.*coeff[3][i]*t+2.*coeff[2][i])*t+coeff[1][i];
  }
  dBrdx_/=dx;
  dBrdz_/=dz;

  // Next compute the interpolation for Bz
  temp[0]=B00->Bz;
  temp[1]=B01->Bz;
  temp[2]=B11->Bz;
  temp[3]=B10->Bz;

  temp[8]=B00->dBzdx;
  temp[9]=B01->dBzdx;
  temp[10]=B11->dBzdx;
  temp[11]=B10->dBzdx;

  temp[4]=B00->dBzdz;
  temp[5]=B01->dBzdz;
  temp[6]=B11->dBzdz;
  temp[7]=B10->dBzdz;
   
  temp[12]=B00->dBzdxdz;
  temp[13]=B01->dBzdxdz;
  temp[14]=B11->dBzdxdz;
  temp[15]=B10->dBzdxdz;

  for (i=0;i<16;i++){
    double tmp2=0.0;
    for (k=0;k<16;k++) tmp2+=wt[i][k]*temp[k];
    cl[i]=tmp2;
  }
  m=0;
  for (i=0;i<4;i++)
    for (j=0;j<4;j++) coeff[i][j]=cl[m++];

  Bz_=dBzdx_=dBzdz_=0.;
  for (i=3;i>=0;i--){
    Bz_=t*Bz_+((coeff[i][3]*u+coeff[i][2])*u+coeff[i][1])*u+coeff[i][0];
    dBzdx_=t*dBzdx_+(3.*coeff[i][3]*u+2.*coeff[i][2])*u+coeff[i][1];
    dBzdz_=u*dBzdz_+(3.*coeff[3][i]*t+2.*coeff[2][i])*t+coeff[1][i];
  }
  dBzdx_/=dx;
  dBzdz_/=dz;

  // Convert r back to x,y components
  double cos_theta = x/r;
  double sin_theta = y/r;
  if(r==0.0){
    cos_theta=1.0;
    sin_theta=0.0;
  }
  // Rotate back into phi direction
  Bx_=Br_*cos_theta;
  By_=Br_*sin_theta;

  dBxdx_ =dBrdx_*cos_theta*cos_theta;
  dBxdy_ =dBrdx_*cos_theta*sin_theta;
  dBxdz_ =dBrdz_*cos_theta;
  dBydx_ = dBrdx_*sin_theta*cos_theta;
  dBydy_ = dBrdx_*sin_theta*sin_theta;
  dBydz_ = dBrdz_*sin_theta;
  dBzdx_ = dBzdx_*cos_theta;
  dBzdy_ = dBzdx_*sin_theta;
  /*
  printf("Grad %f %f %f %f %f %f %f %f %f\n",dBxdx_,dBxdy_,dBxdz_,
	 dBydx_,dBydy_,dBydz_,dBzdx_,dBzdy_,dBzdz_);
  */
}

//-------------
// GetFieldGradient
//-------------
void DMagneticFieldMapCalibDB::GetFieldGradient(double x, double y, double z,
						double &dBxdx, double &dBxdy,
						double &dBxdz,
						double &dBydx, double &dBydy,
						double &dBydz,		
						double &dBzdx, double &dBzdy,
						double &dBzdz) const{
  
  	// Get closest indices for this point
	double r = sqrt(x*x + y*y);
	int index_x = (int)floor((r-xmin)/dx + 0.5);
	//if(index_x<0 || index_x>=Nx)return;	
	if (index_x>=Nx) return;
	else if (index_x<0) index_x=0;
	int index_z = (int)floor((z-zmin)/dz + 0.5);	
	if(index_z<0 || index_z>=Nz)return;
	
	int index_y = 0;

	const DBfieldPoint_t *B = &Btable[index_x][index_y][index_z];

	// Convert r back to x,y components
	double cos_theta = x/r;
	double sin_theta = y/r;
	if(r==0.0){
		cos_theta=1.0;
		sin_theta=0.0;
	}

	// Rotate back into phi direction
	dBxdx = B->dBxdx*cos_theta*cos_theta/dx;
	dBxdy = B->dBxdx*cos_theta*sin_theta/dx;
	dBxdz = B->dBxdz*cos_theta/dz;
	dBydx = B->dBxdx*sin_theta*cos_theta/dx;
	dBydy = B->dBxdx*sin_theta*sin_theta/dx;
	dBydz = B->dBxdz*sin_theta/dz;
	dBzdx = B->dBzdx*cos_theta/dx;
	dBzdy = B->dBzdx*sin_theta/dx;
	dBzdz = B->dBzdz/dz;
	/*
	printf("old Grad %f %f %f %f %f %f %f %f %f\n",dBxdx,dBxdy,dBxdz,
	 dBydx,dBydy,dBydz,dBzdx,dBzdy,dBzdz);
	*/
}



//---------------------------------
// GetField
//---------------------------------
void DMagneticFieldMapCalibDB::GetField(double x, double y, double z, double &Bx, double &By, double &Bz, int method) const
{
	/// This calculates the magnetic field at an arbitrary point
	/// in space using the field map read from the calibaration
	/// database. It interpolates between grid points using the
	/// gradient values calculated in ReadMap (called from the
	/// constructor).

	Bx = By = Bz = 0.0;

	if(Ny>1){
		_DBG_<<"Field map appears to be 3 dimensional. Code is currently"<<endl;
		_DBG_<<"unable to handle this. Treating as phi symmetric using y=0."<<endl;
	}

	// Get closest indices for this point
	double r = sqrt(x*x + y*y);
	int index_x = (int)floor((r-xmin)/dx + 0.5);
	//if(index_x<0 || index_x>=Nx)return;
	if (index_x>=Nx) return;
	else if (index_x<0) index_x=0;
	int index_z = (int)floor((z-zmin)/dz + 0.5);	
	if(index_z<0 || index_z>=Nz)return;
	
	int index_y = 0;

	const DBfieldPoint_t *B = &Btable[index_x][index_y][index_z];

	// Fractional distance between map points.
	double ur = (r - B->x)/dx;
	double uz = (z - B->z)/dz;

	// Use gradient to project grid point to requested position
	double Br;
	Br = B->Bx + B->dBxdx*ur + B->dBxdz*uz;
	Bz = B->Bz + B->dBzdx*ur + B->dBzdz*uz;

	// Convert r back to x,y components
	double cos_theta = x/r;
	double sin_theta = y/r;
	if(r==0.0){
		cos_theta=1.0;
		sin_theta=0.0;
	}

	// Rotate back into phi direction
	Bx = Br*cos_theta;
	By = Br*sin_theta;
}


// Get the z-component of the magnetic field
double DMagneticFieldMapCalibDB::GetBz(double x, double y, double z) const{
  // radial position 
  double r = sqrt(x*x + y*y);

  // Get closest indices for this point
  int index_x = (int)floor((r-xmin)/dx + 0.5);
  if (index_x>=Nx) return 0.;
  else if (index_x<0) index_x=0;
  int index_z = (int)floor((z-zmin)/dz + 0.5);	
  if(index_z<0 || index_z>=Nz)return 0.;
  
  int index_y = 0;
  
  const DBfieldPoint_t *B = &Btable[index_x][index_y][index_z];
  
  // Fractional distance between map points.
  double ur = (r - B->x)/dx;
  double uz = (z - B->z)/dz;
  
  // Use gradient to project grid point to requested position
  return (B->Bz + B->dBzdx*ur + B->dBzdz*uz);
}
