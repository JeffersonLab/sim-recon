// $Id$
//
//    File: DMagneticFieldMapParameterized.cc
// Created: Tue Oct 20 14:06:19 EDT 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.8.0 i386)
//

#include <cmath>
using namespace std;

#include <JANA/JParameterManager.h>
using namespace jana;

#include "DMagneticFieldMapParameterized.h"


// Chebyshev polynomial functions
#define T0() (1)
#define T1(x_1) (x_1)
#define T2(x_2) (2*x_2-1)
#define T3(x_1,x_3) (4*x_3-3*x_1)
#define T4(x_2,x_4) (8*x_4-8*x_2+1)
#define T5(x_5,x_3,x_1) (16*x_5-20*x_3+5*x_1)
#define T6(x_6,x_4,x_2) (32*x_6-48*x_4+18*x_2-1)
#define T7(x_7,x_5,x_3,x_1) (64*x_7-112*x_5+56*x_3-7*x_1)
#define T8(x_8,x_6,x_4,x_2) (128*x_8-256*x_6+160*x_4-32*x_2+1)
#define T9(x_9,x_7,x_5,x_3,x_1) (256*x_9-576*x_7+432*x_5-120*x_3+9*x_1)

// First derivative of Chebyshev polynomial functions
#define dT0dx(x) (0)
#define dT1dx(x) (1)
#define dT2dx(x) (2*2*pow(x,2-1))
#define dT3dx(x) (3*4*pow(x,3-1)-3)
#define dT4dx(x) (4*8*pow(x,4-1)-2*8*pow(x,2-1))
#define dT5dx(x) (5*16*pow(x,5-1)-3*20*pow(x,3-1)+5)
#define dT6dx(x) (6*32*pow(x,6-1)-4*48*pow(x,4-1)+2*18*pow(x,2-1))
#define dT7dx(x) (7*64*pow(x,7-1)-5*112*pow(x,5-1)+3*56*pow(x,3-1)-7)
#define dT8dx(x) (8*128*pow(x,8-1)-6*256*pow(x,6-1)+4*160*pow(x,4-1)-2*32*pow(x,2-1))
#define dT9dx(x) (9*256*pow(x,9-1)-7*576*pow(x,7-1)+5*432*pow(x,5-1)-3*120*pow(x,3-1)+9)

//---------------------------------
// DMagneticFieldMapParameterized    (Constructor)
//---------------------------------
DMagneticFieldMapParameterized::DMagneticFieldMapParameterized(jana::JApplication *japp, string namepath)
{
	int runnumber = 1;
	jcalib = japp->GetJCalibration(runnumber);

	JParameterManager *jparms = japp->GetJParameterManager();
	jparms->SetDefaultParameter("BFIELD_MAP", namepath);
	
	Init(jcalib, namepath);
}

//---------------------------------
// DMagneticFieldMapParameterized    (Constructor)
//---------------------------------
DMagneticFieldMapParameterized::DMagneticFieldMapParameterized(jana::JCalibration *jcalib, string namepath)
{
	Init(jcalib, namepath);
}

//---------------------------------
// ~DMagneticFieldMapParameterized    (Destructor)
//---------------------------------
DMagneticFieldMapParameterized::~DMagneticFieldMapParameterized()
{

}

//---------------------------------
// Init
//---------------------------------
void DMagneticFieldMapParameterized::Init(jana::JCalibration *jcalib, string namepath)
{
	this->jcalib = jcalib;
	
	// First, get the parameters specifying the number of sections and their ranges in z, r
	vector<map<string, string> > sections_map;
	jcalib->Get(namepath, sections_map);
	
	// Loop over entries in the master list and read in those tables
	for(unsigned int k=0; k<sections_map.size(); k++){
		map<string, string> &vals = sections_map[k];

		Dsection section;
		section.namepath = vals["namepath"];
		      section.Bi = vals["Bi"];
		 section.section = atoi(vals["sec"].c_str());
		    section.zmin = atof(vals["zmin"].c_str());
		    section.zmax = atof(vals["zmax"].c_str());
		    section.zmid = atof(vals["zmid"].c_str());
		   section.znorm = atof(vals["znorm"].c_str());
		    section.rmid = atof(vals["rmid"].c_str());
		   section.rnorm = atof(vals["rnorm"].c_str());
		
		// Get parameters for this section
		jcalib->Get(section.namepath, section.pp);
		section.order1 = section.pp.size() - 1;
		section.order2 = section.pp[0].size() - 1;
		
		cout<<"---- parameterized B-field section ---------------"<<endl;
		cout<<"namepath = "<<section.namepath<<endl;
		cout<<"      Bi = "<<section.Bi<<endl;
		cout<<" section = "<<section.section<<endl;
		cout<<"    zmin = "<<section.zmin<<endl;
		cout<<"    zmax = "<<section.zmax<<endl;
		cout<<"    zmid = "<<section.zmid<<endl;
		cout<<"   znorm = "<<section.znorm<<endl;
		cout<<"    rmid = "<<section.rmid<<endl;
		cout<<"   rnorm = "<<section.rnorm<<endl;
		cout<<"  order1 = "<<section.order1<<endl;
		cout<<"  order2 = "<<section.order1<<endl;
		
		// Here we need to transform the order1 by order2 matrix that is in pp
		// using the matrix of Chebyshev coefficients into the similarly 
		// dimensioned matrix that will be used. This is to save all of the
		// CPU to recalculate this on every access.
		
		// Fill non-zero values of C1 with Chebyshev coefficients (from the definition of
		// the Chebyshev polynomials). We use the switch here to fill in only elements
		// for a matrix of order (order1+1) (note there are not breaks). The remaining
		// elements are zero.
		unsigned int &order1 = section.order1;
		DMatrix C1(order1+1, order1+1);
		C1.Zero();
		switch(order1){
			case 9:	C1(9,1)=  +9;	C1(9,3)=-120;	C1(9,5)=+432;	C1(9,7)=-576;	C1(9,9)=+256;
			case 8:	C1(8,0)=  +1;	C1(8,2)= -32;	C1(8,4)=+160;	C1(8,6)=-256;	C1(8,8)=+128;
			case 7:	C1(7,1)=  -7;	C1(7,3)= +56;	C1(7,5)=-112;	C1(7,7)= +64;
			case 6:	C1(6,0)=  -1;	C1(6,2)= +18;	C1(6,4)= -48;	C1(6,6)= +32;
			case 5:	C1(5,1)=  +5;	C1(5,3)= -20;	C1(5,5)= +16;
			case 4:	C1(4,0)=  +1;	C1(4,2)=  -8;	C1(4,4)=  +8;
			case 3:	C1(3,1)=  -3;	C1(3,3)=  +4;
			case 2:	C1(2,0)=  -1;	C1(2,2)=  +2;
			case 1:	C1(1,1)=  +1;
			case 0:	C1(0,0)=  +1;
		}
		TMatrixD C1_t(TMatrixD::kTransposed, C1);

		// Do the same for C2 which is an (order2+1) x (order2+1) matrix
		unsigned int &order2 = section.order2;
		DMatrix C2(order2+1, order2+1);
		C2.Zero();
		switch(order2){
			case 9:	C2(9,1)=  +9;	C2(9,3)=-120;	C2(9,5)=+432;	C2(9,7)=-576;	C2(9,9)=+256;
			case 8:	C2(8,0)=  +1;	C2(8,2)= -32;	C2(8,4)=+160;	C2(8,6)=-256;	C2(8,8)=+128;
			case 7:	C2(7,1)=  -7;	C2(7,3)= +56;	C2(7,5)=-112;	C2(7,7)= +64;
			case 6:	C2(6,0)=  -1;	C2(6,2)= +18;	C2(6,4)= -48;	C2(6,6)= +32;
			case 5:	C2(5,1)=  +5;	C2(5,3)= -20;	C2(5,5)= +16;
			case 4:	C2(4,0)=  +1;	C2(4,2)=  -8;	C2(4,4)=  +8;
			case 3:	C2(3,1)=  -3;	C2(3,3)=  +4;
			case 2:	C2(2,0)=  -1;	C2(2,2)=  +2;
			case 1:	C2(1,1)=  +1;
			case 0:	C2(0,0)=  +1;
		}
		
		// Make matrix of fit parameters
		DMatrix P(order2+1, order1+1);
		for(unsigned int i=0; i<=order1; i++){
			for(unsigned int j=0; j<=order2; j++){
				P(j, i) = section.pp[i][j];
			}
		}

		// Multiply out the matrices and store the final one with the section
		section.Q.ResizeTo(P);
		section.Q = C1_t*P*C2;

		// Copy Q matrix into more efficient format
		section.cc = section.pp;
		for(unsigned int i=0; i<=order1; i++){
			for(unsigned int j=0; j<=order2; j++){
				section.cc[i][j] = section.Q(j,i);
			}
		}

		// Add this section object to list
		string type = vals["Bi"];
		if(type=="Bx")sections_Bx.push_back(section);
		if(type=="Bz")sections_Bz.push_back(section);		
	}
}

//---------------------------------
// GetField
//---------------------------------
void DMagneticFieldMapParameterized::GetField(double x, double y, double z, double &Bx, double &By, double &Bz, int method) const
{
	// default 
	Bx = By = Bz = 0.0;
	double Br = 0.0;

	double r = sqrt(x*x + y*y);

	// Find Bx Dsection object for this z value
	for(unsigned int i=0; i<sections_Bx.size(); i++){
		if(sections_Bx[i].IsInRange(z)){
			Br = sections_Bx[i].Eval(r,z);
			break;
		}
	}

	// Find Bz Dsection object for this z value
	for(unsigned int i=0; i<sections_Bz.size(); i++){
		if(sections_Bz[i].IsInRange(z)){
			Bz = sections_Bz[i].Eval(r,z);
			break;
		}
	}

	// Convert Br to Bx and By
	if(r!=0.0){
		Bx = Br*x/r;
		By = Br*y/r;
	}
}
//---------------------------------
// Get the z-component of the magnetic field
//---------------------------------
double DMagneticFieldMapParameterized::GetBz(double x, double y, double z)const
{

  double r = sqrt(x*x + y*y);
  
  // Find Bz Dsection object for this z value
  for(unsigned int i=0; i<sections_Bz.size(); i++){
    if(sections_Bz[i].IsInRange(z)){
      return (sections_Bz[i].Eval(r,z));
    }
  }
  return 0.;
}


//---------------------------------
// Eval
//---------------------------------
double DMagneticFieldMapParameterized::Dsection::Eval(double &r, double &z) const
{
	// Evaluate the value (either Bz or Bx) at the specified r,z values

	// Initialize value.
	double Bi = 0.0;

	// Convert to "normalized" coordinates (where
	// values range from -1 to +1).
	double r_prime = (r-rmid)/rnorm;
	double z_prime = (z-zmid)/znorm;

#if 0
	// Calculate powers of "r" used by Chebyshev functions
	double &r_1 = r_prime;
	double r_2 = r_1*r_1;
	double r_3 = r_2*r_1;
	double r_4 = r_3*r_1;
	double r_5 = r_4*r_1;
	double r_6 = r_5*r_1;
	double r_7 = r_6*r_1;
	double r_8 = r_7*r_1;
	double r_9 = r_8*r_1;

	// Calculate powers of "z" used by Chebyshev functions
	double &z_1 = z_prime;
	double z_2 = z_1*z_1;
	double z_3 = z_2*z_1;
	double z_4 = z_3*z_1;
	double z_5 = z_4*z_1;
	double z_6 = z_5*z_1;
	double z_7 = z_6*z_1;
	double z_8 = z_7*z_1;
	double z_9 = z_8*z_1;

	// Calculate each parameter based on r_prime and then
	// add the term based on z_prime to Bi.
	for(unsigned int i=0; i<order1; i++){

		// Calculate the 1st level parameter from the 2nd
		// level parameters. Note that this funny switch actually
		// does the sum of 2nd level poly.
		const double *my_pp = &pp[i][0];
		double p = my_pp[0]*T0();
		switch(order2){
			case 9:	p += my_pp[9]*T9(r_9,r_7,r_5,r_3,r_1);
			case 8:	p += my_pp[8]*T8(r_8,r_6,r_4,r_2);
			case 7:	p += my_pp[7]*T7(r_7,r_5,r_3,r_1);
			case 6:	p += my_pp[6]*T6(r_6,r_4,r_2);
			case 5:	p += my_pp[5]*T5(r_5,r_3,r_1);
			case 4:	p += my_pp[4]*T4(r_4,r_2);
			case 3:	p += my_pp[3]*T3(r_3,r_1);
			case 2:	p += my_pp[2]*T2(r_2);
			case 1:	p += my_pp[1]*T1(r_1);
		}

		// Add the i-th poly term to Bi
		switch(i){
			case 0:	Bi += p*T0(); break;
			case 1:	Bi += p*T1(z_1); break;
			case 2:	Bi += p*T2(z_2); break;
			case 3:	Bi += p*T3(z_3,z_1); break;
			case 4:	Bi += p*T4(z_4,z_2); break;
			case 5:	Bi += p*T5(z_5,z_3,z_1); break;
			case 6:	Bi += p*T6(z_6,z_4,z_2); break;
			case 7:	Bi += p*T7(z_7,z_5,z_3,z_1); break;
			case 8:	Bi += p*T8(z_8,z_6,z_4,z_2); break;
			case 9:	Bi += p*T9(z_9,z_7,z_5,z_3,z_1); break;
		}
	}
#else

	double zz = 1;
	for(unsigned int i=1; i<=order2; i++){
		const double *col_ptr = &cc[i][0];
		double v = *col_ptr++;
		double rr = 1;
		for(unsigned int j=1; j<=order1; j++, col_ptr++){
			rr*=r_prime;
			v+=(*col_ptr)*rr;
		}
		
		zz*=z_prime;
		Bi += v*zz;
	}

#if 0
	// Form matrix of r_prime raised to powers from 0 to order2
	DMatrix R(order2+1, 1);
	double rr= R(0,0) = 1;
	for(unsigned int i=1; i<=order2; i++)R(i,0)=(rr*=r_prime); 

	// Form matrix of z_prime raised to powers from 0 to order1
	DMatrix Z(1,order1+1);
	double zz= Z(0,0) = 1;
	for(unsigned int i=1; i<=order1; i++)Z(0,i)=(zz*=z_prime);
	
	// Multiply R and Z vectors using parameters matrix
	DMatrix B = Q*R;
	Bi = B(0,0);
#endif
#endif
//_DBG_<<"Bi="<<Bi<<" B(0,0)="<<B(0,0)<<endl;

	return Bi;
}

//---------------------------------
// GetFieldGradient
//---------------------------------
void DMagneticFieldMapParameterized::GetFieldGradient(double x, double y, double z,
										  double &dBxdx, double &dBxdy,
										  double &dBxdz,
										  double &dBydx, double &dBydy,
										  double &dBydz,
										  double &dBzdx, double &dBzdy,
										  double &dBzdz) const
{

}

//---------------------------------
// GetFieldBicubic
//---------------------------------
void DMagneticFieldMapParameterized::GetFieldBicubic(double x,double y,double z, double &Bx,double &By,double &Bz) const
{
	// Since we calculate the gradient analytically here, a bicubic
	// interpolation will not give a better result than what is
	// returned by the GetFieldGradient method.
	GetField(x,y,z,Bx, By, Bz);
}

//---------------------------------
// GetFieldAndGradient
//---------------------------------
void DMagneticFieldMapParameterized::GetFieldAndGradient(double x,double y,double z,
				 double &Bx,
				 double &By,
				 double &Bz,
				 double &dBxdx, 
				 double &dBxdy,
				 double &dBxdz,
				 double &dBydx, 
				 double &dBydy,
				 double &dBydz,
				 double &dBzdx, 
				 double &dBzdy,
				 double &dBzdz) const
{
	// Just defer to the individual functions here.
	GetField(x,y,z,Bx, By, Bz);
	GetFieldGradient(x, y, z, dBxdx, dBxdy, dBxdz, dBydx, dBydy, dBydz, dBzdx, dBzdy, dBzdz);
}

