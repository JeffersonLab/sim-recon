

#include <iostream>
#include <vector>
using namespace std;

#include <JANA/JCalibrationFile.h>
#include <HDGEOMETRY/DMagneticFieldMapCalibDB.h>

extern "C" {
#include "calibDB.h"
};

DMagneticFieldMap *Bmap=NULL;
static JCalibration *jcalib=NULL;

//----------------
// initcalibdb_
//----------------
void initcalibdb_(void)
{
	ios::sync_with_stdio(true);
	
	// Create a JCalibration object using the JANA_CALIB_URL environment variable
	// Right now, we hardwire this to use JCalibrationFile.
	const char *url = getenv("JANA_CALIB_URL");
	if(!url){
		_DBG_<<"JANA_CALIB_URL environment not set."<<endl;
		exit(-1);
	}
	jcalib = new JCalibrationFile(url, 1, "");
	
	// Read in the field map from the calibration DB
	Bmap = new DMagneticFieldMapCalibDB(jcalib);

	
}

//----------------
// gufld2_
//----------------
void gufld2_(float *r, float *B)
{
	/// Wrapper function to allow the FORTRAN gufld routine to
	/// use the C++ class DMagneticFieldMap to access the 
	/// B-field.

	if(!Bmap){
		_DBG_<<"Call to gufld2_ when Bmap not intialized! Exiting."<<endl;
		exit(-1);
	}
	
	double x = r[0];
	double y = r[1];
	double z = r[2];
	double Bx, By, Bz;
	
	Bmap->GetField(x, y, z, Bx, By, Bz);

	B[0] = Bx;
	B[1] = By;
	B[2] = Bz;
}

//----------------
// GetCalib
//----------------
int GetCalib(const char* namepath, unsigned int *Nvals, float* vals)
{
	/// C-callable routine for accessing calibraion constants.
	/// The values specified by "namepath" will be read into the array
	/// "vals". The "vals" array should have enough memory allocated
	/// to hold *Nvals elements. If not, only the first *Nvals elements
	/// will be copied and a non-zero value returned. If the number
	/// of values in the database are less than *Nvals, then all values
	/// are copied, *Nvals is updated to reflect the number of valid
	/// elements in "vals", and a value of 0 is returned.
	
	if(!jcalib){
		_DBG_<<"ERROR - GetCalib() called when jcalib not set!"<<endl;
		_DBG_<<"ERROR - request for \""<<namepath<<"\""<<endl;
		_DBG_<<"ERROR - Exiting ..."<<endl;
		exit(-1);
	}

	vector<float> vvals;
	jcalib->Get(namepath, vvals);
	if(vvals.size()<*Nvals)*Nvals = vvals.size();
	for(unsigned int i=0; i<*Nvals; i++)vals[i] = vvals[i];
	
	return vvals.size()>*Nvals; // return 0 if OK, 1 if not
}

//----------------
// GetLorentzDefelections
//----------------
void GetLorentzDefelections(float *lorentz_x, float *lorentz_z, float **lorentz_nx, float **lorentz_nz
	, const unsigned int Nxpoints, const unsigned int Nzpoints)
{
	/// C-callable routine for accessing calibraion constants.
	/// The values specified by "namepath" will be read into the array
	/// "vals". The "vals" array should have enough memory allocated
	/// to hold *Nvals elements. If not, only the first *Nvals elements
	/// will be copied and a non-zero value returned. If the number
	/// of values in the database are less than *Nvals, then all values
	/// are copied, *Nvals is updated to reflect the number of valid
	/// elements in "vals", and a value of 0 is returned.
	
	// Make sure jcalib is set
	if(!jcalib){
		_DBG_<<"ERROR - GetLorentzDefelections() called when jcalib not set!"<<endl;
		_DBG_<<"ERROR - Exiting ..."<<endl;
		exit(-1);
	}

	// Get constants and do basic check on number of elements
	vector< map<string, float> > tvals;
	jcalib->Get("FDC/lorentz_deflections", tvals);
	if(tvals.size() != Nxpoints*Nzpoints){
		_DBG_<<"ERROR - GetLorentzDefelections() number of elements in calib DB"<<endl;
		_DBG_<<"ERROR - not the same as expected. DB="<<tvals.size()<<" expected"<<Nxpoints*Nzpoints<<endl;
		_DBG_<<"ERROR - Exiting ..."<<endl;
		exit(-1);
	}
	
	// Notify user
	cout<<"Read "<<tvals.size()<<" values from FDC/lorentz_deflections in calibDB"<<endl;
	cout<<"   lorentz_deflections columns (alphabetical): ";
	map<string,float>::iterator iter;
	for(iter=tvals[0].begin(); iter!=tvals[0].end(); iter++)cout<<iter->first<<" ";
	cout<<endl;
	
	// Copy values into tables. We preserve the order since that is how it was
	// originally done in hitFDC.c
	for(unsigned int i=0; i<tvals.size(); i++){
		map<string, float> &row = tvals[i];
		unsigned int xindex = i/Nzpoints;
		unsigned int zindex = i%Nzpoints;
		lorentz_x[xindex] = row["x"];
		lorentz_z[zindex] = row["z"];
		lorentz_nx[xindex][zindex] = row["nx"];
		lorentz_nz[xindex][zindex] = row["nz"];
	}	
}

