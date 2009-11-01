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
#define T0(x) (1)
#define T1(x) (x)
#define T2(x) (2*pow(x,2)-1)
#define T3(x) (4*pow(x,3)-3*x)
#define T4(x) (8*pow(x,4)-8*pow(x,2)+1)
#define T5(x) (16*pow(x,5)-20*pow(x,3)+5*x)
#define T6(x) (32*pow(x,6)-48*pow(x,4)+18*pow(x,2)-1)
#define T7(x) (64*pow(x,7)-112*pow(x,5)+56*pow(x,3)-7*x)
#define T8(x) (128*pow(x,8)-256*pow(x,6)+160*pow(x,4)-32*pow(x,2)+1)
#define T9(x) (256*pow(x,9)-576*pow(x,7)+432*pow(x,5)-120*pow(x,3)+9*x)

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
	for(unsigned int i=0; i<sections_map.size(); i++){
		map<string, string> &vals = sections_map[i];

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
		
		string type = vals["Bi"];
		if(type=="Bx")sections_Bx.push_back(section);
		if(type=="Bz")sections_Bz.push_back(section);
		
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

	// Calculate each parameter based on r_prime and then
	// add the term based on z_prime to Bi.
	for(unsigned int i=0; i<order1; i++){

		// Calculate the 1st level parameter from the 2nd
		// level parameters. Note that this funny switch actually
		// does the sum of 2nd level poly.
		const double *my_pp = &pp[i][0];
		double p = my_pp[0]*T0(r_prime);
		switch(order2){
			case 9:	p += my_pp[9]*T9(r_prime);
			case 8:	p += my_pp[8]*T8(r_prime);
			case 7:	p += my_pp[7]*T7(r_prime);
			case 6:	p += my_pp[6]*T6(r_prime);
			case 5:	p += my_pp[5]*T5(r_prime);
			case 4:	p += my_pp[4]*T4(r_prime);
			case 3:	p += my_pp[3]*T3(r_prime);
			case 2:	p += my_pp[2]*T2(r_prime);
			case 1:	p += my_pp[1]*T1(r_prime);
		}

		// Add the i-th poly term to Bi
		switch(i){
			case 0:	Bi += p*T0(z_prime); break;
			case 1:	Bi += p*T1(z_prime); break;
			case 2:	Bi += p*T2(z_prime); break;
			case 3:	Bi += p*T3(z_prime); break;
			case 4:	Bi += p*T4(z_prime); break;
			case 5:	Bi += p*T5(z_prime); break;
			case 6:	Bi += p*T6(z_prime); break;
			case 7:	Bi += p*T7(z_prime); break;
			case 8:	Bi += p*T8(z_prime); break;
			case 9:	Bi += p*T9(z_prime); break;
		}
	}

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

