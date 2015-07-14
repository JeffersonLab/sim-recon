// $Id$
//
//    File: DMagneticFieldMapPS2DMap.cc

#include <sys/stat.h>
#include <cmath>
using namespace std;
#ifdef HAVE_EVIO
#include <evioFileChannel.hxx>
#include <evioUtil.hxx>
using namespace evio;
#endif

#include "DMagneticFieldMapPS2DMap.h"

//class DApplication;
#include <HDGEOMETRY/DGeometry.h>
#include <DANA/DApplication.h>

//---------------------------------
// DMagneticFieldMapPS2DMap    (Constructor)
//---------------------------------
DMagneticFieldMapPS2DMap::DMagneticFieldMapPS2DMap(JApplication *japp, unsigned int runnumber, string namepath)
{
	jcalib = japp->GetJCalibration(runnumber);
	jresman = japp->GetJResourceManager(runnumber);
	DApplication *dapp = dynamic_cast<DApplication *>(japp);
	geom = dapp->GetDGeometry(runnumber);
	
	JParameterManager *jparms = japp->GetJParameterManager();
	jparms->SetDefaultParameter("PSBFIELD_MAP", namepath);
	
	int Npoints = ReadMap(namepath, runnumber); 
	if(Npoints==0){
		_DBG_<<"Error getting JCalibration object for magnetic field!"<<endl;
		japp->Quit();
	}
}


//---------------------------------
// DMagneticFieldMapPS2DMap    (Constructor)
//---------------------------------
DMagneticFieldMapPS2DMap::DMagneticFieldMapPS2DMap(JCalibration *jcalib, string namepath)
{
	this->jcalib = jcalib;
	this->geom = NULL;
	this->jresman = NULL;
	if(ReadMap(namepath)==0){
		_DBG_<<"Error getting JCalibration object for magnetic field!"<<endl;
		exit(1);
	} 
}

//---------------------------------
// ~DMagneticFieldMapPS2DMap    (Destructor)
//---------------------------------
DMagneticFieldMapPS2DMap::~DMagneticFieldMapPS2DMap()
{

}

//---------------------------------
// ReadMap
//---------------------------------
int DMagneticFieldMapPS2DMap::ReadMap(string namepath, int runnumber, string context)
{
  /// Read the magnetic field map in from the calibration database.
  /// This will read in the map and figure out the number of grid
  /// points in each direction (x,y, and z) and the range in each.
  /// The gradiant of the field is calculated for all but the most
  /// exterior points and saved to use in later calls to GetField(...).

  // Read in map from calibration database. This should really be
  // built into a run-dependent geometry framework, but for now
  // we do it this way. 
  if(!jcalib){
    jerr << "ERROR: jcalib pointer is NULL in DMagneticFieldMapPS2DMap::ReadMap() !" << endl;
    return 0;
  }
  if(!jresman){
    jerr << "ERROR: jresman pointer is NULL in DMagneticFieldMapPS2DMap::ReadMap() !" << endl;
    return 0;
  }
  if(!geom){
    jerr << "ERROR: geom pointer is NULL in DMagneticFieldMapPS2DMap::ReadMap() !" << endl;
    return 0;
  }
  
  jout<<endl;
  jout<<"Reading Pair Spectrometer Magnetic field map from "<<namepath<<" ..."<<endl;
  vector< vector<float> > Bmap;

  // Try getting it as a JANA resource.
  jresman->Get(namepath, Bmap);
  
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
  cout<<" Nx="<<Nx;
  cout<<" Ny="<<Ny;
  cout<<" Nz="<<Nz;
  cout<<" )  at 0x"<<hex<<(unsigned long)this<<dec<<endl;
  
  // Create 3D vector so we can index the values by [x][y][z]
  vector<DBfieldPoint_t> zvec(Nz);
  vector< vector<DBfieldPoint_t> > yvec;
  for(int i=0; i<Ny; i++)yvec.push_back(zvec);
  for(int i=0; i<Nx; i++)Btable.push_back(yvec);
	
  // Distance between map points for r and z
  dx = (xmax-xmin)/(double)(Nx-1);
  dy = (ymax-ymin)/(double)((Ny < 2)? 1 : Ny-1);
  dz = (zmax-zmin)/(double)(Nz-1);
  
  one_over_dx=1./dx;
  one_over_dz=1./dz;

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
    // convert to T
    b->Bx = a[3]/10000.;
    b->By = a[4]/10000.;
    b->Bz = a[5]/10000.;
    //b->Bx = a[3];
    //b->By = a[4];
    //b->Bz = a[5];
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

  // fix this
  vector<double>pair_spect_origin;
  geom->Get("//posXYZ[@volume='PairSpectrometer']/@X_Y_Z",pair_spect_origin);
  vector<double>pair_magnet_origin;
  geom->Get("//posXYZ[@volume='PairMagnet']/@X_Y_Z",pair_magnet_origin);

  z_shift = pair_spect_origin[2]-pair_magnet_origin[2];

  //cout << " PS volume origin = *" << pair_spect_origin[0] << ", " << pair_spect_origin[1] << ", " << pair_spect_origin[2] << ")" << endl;
  //cout << " PS magnet volume origin = *" << pair_magnet_origin[0] << ", " << pair_magnet_origin[1] << ", " << pair_magnet_origin[2] << ")" << endl;

  //jout << "Map Z shift = " << z_shift << endl;

  return Bmap.size();
}

// Use bicubic interpolation to find the field at the point (x,y).  
//See Numerical Recipes in C (2nd ed.), pp.125-127.
void DMagneticFieldMapPS2DMap::GetFieldBicubic(double x,double y,double z,
						   double &Bx_,double &By_,
						   double &Bz_) const{
  // transform from hall coordinates to map coordinates
  z -= z_shift;

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

  // Interpolate using coarse grid (fine grid not implemented yet)
  // Get closest indices for this point
    int index_x = (int)floor((x-xmin)*one_over_dx + 0.5);
    if(index_x<0 || index_x>=Nx)return;
    //if (index_x>=Nx) return;
    else if (index_x<0) index_x=0;
    int index_z = (int)floor((z-zmin)*one_over_dz + 0.5);	
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
    
    // First compute the interpolation for Bx
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
    
    double t=(z - B00->z)*one_over_dz;
    double u=(x - B00->x)*one_over_dx;   
    for (i=3;i>=0;i--){
      Bx_=t*Bx_+((coeff[i][3]*u+coeff[i][2])*u+coeff[i][1])*u+coeff[i][0];
    }
    
    // Next compute the interpolation for By
    temp[0]=B00->By;
    temp[1]=B01->By;
    temp[2]=B11->By;
    temp[3]=B10->By;
    
    temp[8]=B00->dBydx;
    temp[9]=B01->dBydx;
    temp[10]=B11->dBydx;
    temp[11]=B10->dBydx;
    
    temp[4]=B00->dBydy;
    temp[5]=B01->dBydy;
    temp[6]=B11->dBydy;
    temp[7]=B10->dBydy;
    
    temp[12]=B00->dBydxdy;
    temp[13]=B01->dBydxdy;
    temp[14]=B11->dBydxdy;
    temp[15]=B10->dBydxdy;
    
    for (i=0;i<16;i++){
      double tmp2=0.0;
      for (k=0;k<16;k++) tmp2+=wt[i][k]*temp[k];
      cl[i]=tmp2;
    }
    m=0;
    for (i=0;i<4;i++)
      for (j=0;j<4;j++) coeff[i][j]=cl[m++];
    
    By_=0.;
    for (i=3;i>=0;i--){
      By_=t*By_+((coeff[i][3]*u+coeff[i][2])*u+coeff[i][1])*u+coeff[i][0];
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

}


// Find the field and field gradient at the point (x,y,z).  
void DMagneticFieldMapPS2DMap::GetFieldAndGradient(double x,double y,double z,
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
  // transform from hall coordinates to map coordinates
  z -= z_shift;

  // radial distance
  double r = sqrt(x*x + y*y);
  if (r>xmax || z>zmax || z<zmin){
    Bx_=0.0,By_=0.0,Bz_=0.0;
    dBxdx_=0.0,dBxdy_=0.0,dBxdz_=0.0;
    dBydx_=0.0,dBydy_=0.0,dBydz_=0.0;
    dBzdx_=0.0,dBzdy_=0.0,dBzdz_=0.0;
    return;
  }

  // Initialize z-component
  Bz_=0.;
  
  // Get closest indices for this point
  int index_x = static_cast<int>((x-xmin)*one_over_dx);
  int index_z = static_cast<int>((z-zmin)*one_over_dz);	
  int index_y = 0;
  
  if(index_x>=0 && index_x<Nx && index_z>=0 && index_z<Nz){
    const DBfieldPoint_t *B = &Btable[index_x][index_y][index_z];
	  
    // Fractional distance between map points.
    double ux = (x - B->x)*one_over_dx;
    double uz = (z - B->z)*one_over_dz;
    
    // Use gradient to project grid point to requested position
    Bx_ = B->Bx+B->dBxdx*ux+B->dBxdz*uz;
    By_ = B->By;
    Bz_ = B->Bz+B->dBzdx*ux+B->dBzdz*uz;
    dBxdx_=B->dBxdx;
    dBxdy_=B->dBxdy;
    dBxdz_=B->dBxdz;
    dBydx_=B->dBydx;
    dBydy_=B->dBydy;
    dBydz_=B->dBydz;
    dBzdx_=B->dBzdx;
    dBzdy_=B->dBzdy;
    dBzdz_=B->dBzdz;

    dBydx_= dBxdy_;
  
  /*
  printf("Grad %f %f %f %f %f %f %f %f %f\n",dBxdx_,dBxdy_,dBxdz_,
	 dBydx_,dBydy_,dBydz_,dBzdx_,dBzdy_,dBzdz_);
  */
  }
}
  


//-------------
// GetFieldGradient
//-------------
void DMagneticFieldMapPS2DMap::GetFieldGradient(double x, double y, double z,
						double &dBxdx, double &dBxdy,
						double &dBxdz,
						double &dBydx, double &dBydy,
						double &dBydz,		
						double &dBzdx, double &dBzdy,
						double &dBzdz) const{
  
  // transform from hall coordinates to map coordinates
  z -= z_shift;

  	// Get closest indices for this point
	int index_x = (int)floor((x-xmin)*one_over_dx + 0.5);
	if(index_x<0 || index_x>=Nx)return;	
	//if (index_x>=Nx) return;
	else if (index_x<0) index_x=0;
	int index_z = (int)floor((z-zmin)*one_over_dz + 0.5);	
	if(index_z<0 || index_z>=Nz)return;
	
	int index_y = 0;

	const DBfieldPoint_t *B = &Btable[index_x][index_y][index_z];

	dBxdx = B->dBxdx;
	dBxdy = B->dBxdx;
	dBxdz = B->dBxdz;
	dBydx = B->dBydx;
	dBydy = B->dBydx;
	dBydz = B->dBydz;
	dBzdx = B->dBzdx;
	dBzdy = B->dBzdx;
	dBzdz = B->dBzdz;
}



//---------------------------------
// GetField
//---------------------------------
void DMagneticFieldMapPS2DMap::GetField(double x, double y, double z, double &Bx, double &By, double &Bz, int method) const
{
	//cout << "in GetField()..." << endl;

	/// This calculates the magnetic field at an arbitrary point
	/// in space using the field map read from the calibaration
	/// database. It interpolates between grid points using the
	/// gradient values calculated in ReadMap (called from the
	/// constructor).
  
	//cout << " hall coords:  x = " << x << "  z = " << z << endl;

	// transform from hall coordinates to map coordinates
	z -= z_shift;

	//cout << " local coords:  x = " << x << "  z = " << z << endl;

	Bx = By = Bz = 0.0;

	if(Ny>1){
		_DBG_<<"Field map appears to be 3 dimensional. Code is currently"<<endl;
		_DBG_<<"unable to handle this. Assuming y=0."<<endl;
	}

	if (x<xmin || x>xmax || z>zmax || z<zmin){
	  return;
	}

	// Interpolate using coarse grid 
	// Get closest indices for this point
	int index_x = static_cast<int>((x-xmin)*one_over_dx);
	if (index_x<0 || index_x>=Nx) return;
	
	int index_z = static_cast<int>((z-zmin)*one_over_dz);	
	if(index_z<0 || index_z>=Nz)return;
	  
	int index_y = 0;
	  
	//cout << " indexes:  x = " << index_x << "  z = " << index_z << endl;

	const DBfieldPoint_t *B = &Btable[index_x][index_y][index_z];

	//cout << " mag field:  Bx = " << B->x << "  By = " << B->y << "  Bz = " << B->z << endl;
	  
	// Fractional distance between map points.
	double ux = (x - B->x)*one_over_dx;
	double uz = (z - B->z)*one_over_dz;
	
	// Use gradient to project grid point to requested position
	Bx = B->Bx+B->dBxdx*ux+B->dBxdz*uz;
	By = B->By;
	Bz = B->Bz+B->dBzdx*ux+B->dBzdz*uz;
}

//---------------------------------
// GetField
//---------------------------------
void DMagneticFieldMapPS2DMap::GetField(const DVector3 &pos,DVector3 &Bout) const
{
	/// This calculates the magnetic field at an arbitrary point
	/// in space using the field map read from the calibaration
	/// database. It interpolates between grid points using the
	/// gradient values calculated in ReadMap (called from the
	/// constructor).
  

	double Bz=0.,Bx=0.0,By=0.0;

	if(Ny>1){
		_DBG_<<"Field map appears to be 3 dimensional. Code is currently"<<endl;
		_DBG_<<"unable to handle this. Assuming y=0."<<endl;
	}

	// radial position and angles
	double x = pos.x();
	//double y = pos.y();
	double z = pos.z();
	// transform from hall coordinates to map coordinates
	z -= z_shift;
	if (x<xmin || x>xmax || z>zmax || z<zmin){
	  return;
	}

	// Interpolate using coarse grid 
	// Get closest indices for this point
	int index_x = static_cast<int>((x-xmin)*one_over_dx);
	if (index_x<0 || index_x>=Nx) return;
	
	int index_z = static_cast<int>((z-zmin)*one_over_dz);	
	if(index_z<0 || index_z>=Nz)return;
	  
	int index_y = 0;
	  
	const DBfieldPoint_t *B = &Btable[index_x][index_y][index_z];
	  
	// Fractional distance between map points.
	double ux = (x - B->x)*one_over_dx;
	double uz = (z - B->z)*one_over_dz;
	
	// Use gradient to project grid point to requested position
	Bx = B->Bx+B->dBxdx*ux+B->dBxdz*uz;
	By = B->By;
	Bz = B->Bz+B->dBzdx*ux+B->dBzdz*uz;

	// Set result
	Bout.SetXYZ(Bx,By,Bz);
}




