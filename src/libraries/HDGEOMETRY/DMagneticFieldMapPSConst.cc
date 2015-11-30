// $Id$
//
//    File: DMagneticFieldMapPSConst.cc
//

#include <cmath>
using namespace std;

#include "DMagneticFieldMapPSConst.h"



//---------------------------------
// DMagneticFieldMapPSConst    (Constructor)
//---------------------------------
DMagneticFieldMapPSConst::DMagneticFieldMapPSConst(JApplication *japp, string namepath)
{
	int32_t runnumber = 1;
	jcalib = japp->GetJCalibration(runnumber);
	if(GetValues(namepath, runnumber)==0){
		_DBG_<<"Error getting JCalibration object for magnetic field!"<<endl;
		japp->Quit();
	}
}

//---------------------------------
// DMagneticFieldMapPSConst    (Constructor)
//---------------------------------
DMagneticFieldMapPSConst::DMagneticFieldMapPSConst(JCalibration *jcalib, string namepath)
{
	this->jcalib = jcalib;
	if(GetValues(namepath)==0){
		_DBG_<<"Error getting JCalibration object for magnetic field!"<<endl;
		exit(-1);
	} 
}

//---------------------------------
// DMagneticFieldMapPSConst    (Constructor)
//---------------------------------
DMagneticFieldMapPSConst::DMagneticFieldMapPSConst(double Bx, double By, double Bz)
{
	this->jcalib = NULL;

	this->Bx = Bx;
	this->By = By;
	this->Bz = Bz;
}

//---------------------------------
// ~DMagneticFieldMapPSConst    (Destructor)
//---------------------------------
DMagneticFieldMapPSConst::~DMagneticFieldMapPSConst()
{

}

//---------------------------------
// GetValues
//---------------------------------
int DMagneticFieldMapPSConst::GetValues(string namepath, int32_t runnumber, string context)
{
	/// Read the parameters for the constant magnetic field map from the calibration database.

	if(!jcalib)return 0;
	
	cout<<"Reading Constant Magnetic field values from "<<namepath<<" ..."<<endl;
	map<string,double> vals;
	jcalib->Get(namepath, vals);
	if(vals.size()==0)return 0;
	
	Bx = vals["Bx"];
	By = vals["By"];
	Bz = vals["Bz"];
	cout<<"   Bx="<<Bx<<"  By="<<By<<"  Bz="<<Bz<<endl;
	
	return vals.size();
}


//-------------
// GetFieldGradient
//-------------
void DMagneticFieldMapPSConst::GetFieldGradient(double x, double y, double z,
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
void DMagneticFieldMapPSConst::GetField(double x, double y, double z, double &Bx, double &By, double &Bz, int method) const
{
	/// This calculates the magnetic field at an arbitrary point
	/// in space using the constant field map parameters read from the calibaration
	/// database.

	Bx = this->Bx;
	By = this->By;
	Bz = this->Bz;
}

//---------------------------------
// GetField
//---------------------------------
void DMagneticFieldMapPSConst::GetField(const DVector3 &pos,
				      DVector3 &Bout) const
{
	/// This calculates the magnetic field at an arbitrary point
	/// in space using the constant field map parameters read from the calibaration
	/// database.

	Bout.SetXYZ(this->Bx,this->By,this->Bz);
}


void DMagneticFieldMapPSConst::GetFieldBicubic(double x,double y,double z,
					     double &Bx,double &By,double &Bz) 
  const{
  GetField(x,y,z,Bx,By,Bz);
}


void DMagneticFieldMapPSConst::GetFieldAndGradient(double x,double y,double z,
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
