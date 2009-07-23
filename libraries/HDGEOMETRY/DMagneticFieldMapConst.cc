// $Id$
//
//    File: DMagneticFieldMapConst.cc
// Created: Fri Nov  7 04:01:28 EST 2008
// Creator: davidl (on Darwin Amelia.local 8.11.1 i386)
//

#include <cmath>
using namespace std;

#include "DMagneticFieldMapConst.h"



//---------------------------------
// DMagneticFieldMapConst    (Constructor)
//---------------------------------
DMagneticFieldMapConst::DMagneticFieldMapConst(JApplication *japp, string namepath)
{
	int runnumber = 1;
	jcalib = japp->GetJCalibration(runnumber);
	if(GetValues(namepath, runnumber)==0){
		_DBG_<<"Error getting JCalibration object for magnetic field!"<<endl;
		japp->Quit();
	}
}

//---------------------------------
// DMagneticFieldMapConst    (Constructor)
//---------------------------------
DMagneticFieldMapConst::DMagneticFieldMapConst(JCalibration *jcalib, string namepath)
{
	this->jcalib = jcalib;
	if(GetValues(namepath)==0){
		_DBG_<<"Error getting JCalibration object for magnetic field!"<<endl;
		exit(-1);
	} 
}

//---------------------------------
// ~DMagneticFieldMapConst    (Destructor)
//---------------------------------
DMagneticFieldMapConst::~DMagneticFieldMapConst()
{

}

//---------------------------------
// GetValues
//---------------------------------
int DMagneticFieldMapConst::GetValues(string namepath, int runnumber, string context)
{
	/// Read the parameters for the constant magnetic field map from the calibration database.

	if(!jcalib)return 0;
	
	cout<<"Reading Constant Magnetic field values from "<<namepath<<" ..."<<endl;
	map<string,double> vals;
	jcalib->Get(namepath, vals);
	if(vals.size()==0)return 0;
	
	Br = vals["Br"];
	Bphi = vals["Bphi"];
	Bz = vals["Bz"];
	cout<<"   Br="<<Br<<"  Bphi="<<Bphi<<"  Bz="<<Bz<<endl;
	
	return vals.size();
}


//-------------
// GetFieldGradient
//-------------
void DMagneticFieldMapConst::GetFieldGradient(double x, double y, double z,
						double &dBxdx, double &dBxdy,
						double &dBxdz,
						double &dBydx, double &dBydy,
						double &dBydz,		
						double &dBzdx, double &dBzdy,
						double &dBzdz) const{
  

	// Constant field has zero gradient
	dBxdx = 0.0;
	dBxdy = 0.0;
	dBxdz = 0.0;
	dBydx = 0.0;
	dBydy = 0.0;
	dBydz = 0.0;
	dBzdx = 0.0;
	dBzdy = 0.0;
	dBzdz = 0.0;
}



//---------------------------------
// GetField
//---------------------------------
void DMagneticFieldMapConst::GetField(double x, double y, double z, double &Bx, double &By, double &Bz, int method) const
{
	/// This calculates the magnetic field at an arbitrary point
	/// in space using the constat field map parameters read from the calibaration
	/// database.

	if(Br!=0.0){
		double r = sqrt(x*x + y*y);
		double cos_phi = x/r;
		double sin_phi = y/r;
		if(r==0.0){
			cos_phi=1.0;
			sin_phi=0.0;
		}
		Bx = Br*cos_phi;
		By = Br*sin_phi;
	}else{
		Bx=By=0.0;
	}
	Bz = this->Bz;
}



void DMagneticFieldMapConst::GetFieldBicubic(double x,double y,double z,
					     double &Bx,double &By,double &Bz) 
  const{
  GetField(x,y,z,Bx,By,Bz);
}


void DMagneticFieldMapConst::GetFieldAndGradient(double x,double y,double z,
						  double &Bx,double &By,
						  double &Bz,
						  double &dBxdx, double &dBxdy,
						  double &dBxdz,
						  double &dBydx, double &dBydy,
						  double &dBydz,
						  double &dBzdx, double &dBzdy,
						  double &dBzdz) const{
  GetField(x,y,z,Bx,By,Bz);
  // Constant field has zero gradient
  dBxdx = 0.0;
  dBxdy = 0.0;
  dBxdz = 0.0;
  dBydx = 0.0;
  dBydy = 0.0;
  dBydz = 0.0;
  dBzdx = 0.0;
  dBzdy = 0.0;
  dBzdz = 0.0;
}
