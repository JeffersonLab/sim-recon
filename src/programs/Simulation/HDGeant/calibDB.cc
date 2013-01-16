
#include <stdlib.h>
#include <dlfcn.h>
#include <unistd.h>

#include <iostream>
#include <vector>
#include <string>
using namespace std;

#include <JANA/JCalibrationFile.h>
#include <HDGEOMETRY/DMagneticFieldMapCalibDB.h>
#include <HDGEOMETRY/DMagneticFieldMapConst.h>
#include "HDGEOMETRY/DMagneticFieldMapSpoiled.h"
#include "HDGEOMETRY/DMagneticFieldMapParameterized.h"

extern "C" {
#include "calibDB.h"
};
#include "controlparams.h"


extern "C" int hddsgeant3_runtime_(void);  // called from uginit.F. defined in calibDB.cc
extern "C" void md5geom_(char *md5);
void init_runtime_xml(void);
void md5geom_runtime(char *md5);
extern "C" const char* GetMD5Geom(void);

DMagneticFieldMap *Bmap=NULL;
static JCalibration *jcalib=NULL;
//static void *dlgeom_handle=NULL;
//string HDDS_XML = "$HDDS_HOME/main_HDDS.xml";

extern "C" {
	void md5geom_wrapper_(char *md5);
}


//----------------
// initcalibdb_
//----------------
void initcalibdb_(char *bfield_type, char *bfield_map)
{
	ios::sync_with_stdio(true);

	if(japp) jcalib = japp->GetJCalibration(1);

	if(!jcalib){
		jerr << "----- Couldn't create JCalibration object!!!" << endl;
	}

#if 0
	// Create a JCalibration object using the JANA_CALIB_URL environment variable
	// Right now, we hardwire this to use JCalibrationFile.
	const char *url = getenv("JANA_CALIB_URL");
	if(!url){
		_DBG_<<"JANA_CALIB_URL environment not set."<<endl;
		exit(-1);
	}
	jcalib = new JCalibrationFile(url, 1, "");
#endif
	
	// The actual DMagneticFieldMap subclass can be specified in
	// the control.in file. Since it is read in as integers of
	// "MIXED" format through ffkey though (who knows what that
	// means!) then there can be trailing white spaceat the end
	// of the string. Here, we replace terminate the string with
	// a null to eliminate that white space.
	while(strlen(bfield_type)>0 && bfield_type[strlen(bfield_type)-1]==' ')bfield_type[strlen(bfield_type)-1] = 0;
	while(strlen(bfield_map)>0 && bfield_map[strlen(bfield_map)-1]==' ')bfield_map[strlen(bfield_map)-1] = 0;
	
	// Read in the field map from the appropriate source
	if(bfield_type[0] == 0)strcpy(bfield_type, "CalibDB");
	string bfield_type_str(bfield_type);
	if(bfield_type_str=="CalibDB"){
		if(strlen(bfield_map))
			Bmap = new DMagneticFieldMapCalibDB(jcalib, bfield_map);
		else
			Bmap = new DMagneticFieldMapCalibDB(jcalib);
	}else if(bfield_type_str=="Const"){
		if(strlen(bfield_map))
			Bmap = new DMagneticFieldMapConst(jcalib, bfield_map);
		else
			Bmap = new DMagneticFieldMapConst(jcalib);
	}else if(bfield_type_str=="Spoiled"){
		if(strlen(bfield_map))
			Bmap = new DMagneticFieldMapSpoiled(jcalib, bfield_map);
		else
			Bmap = new DMagneticFieldMapSpoiled(jcalib);
	}else if(bfield_type_str=="Parameterized"){
		if(strlen(bfield_map))
			Bmap = new DMagneticFieldMapParameterized(jcalib, bfield_map);
		else
			Bmap = new DMagneticFieldMapParameterized(jcalib);
	}else{
		_DBG_<<" Unknown DMagneticFieldMap subclass \"DMagneticFieldMap"<<bfield_type_str<<"\" !!"<<endl;
		exit(-1);
	}	
}

//----------------
// gufld_db_
//----------------
void gufld_db_(float *r, float *B)
{
	/// Wrapper function to allow the FORTRAN gufld routine to
	/// use the C++ class DMagneticFieldMap to access the 
	/// B-field.

	if(!Bmap){
		_DBG_<<"Call to gufld_db when Bmap not intialized! Exiting."<<endl;
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
// GetLorentzDeflections
//----------------
void GetLorentzDeflections(float *lorentz_x, float *lorentz_z, float **lorentz_nx, float **lorentz_nz
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
		_DBG_<<"ERROR - GetLorentzDeflections() called when jcalib not set!"<<endl;
		_DBG_<<"ERROR - Exiting ..."<<endl;
		exit(-1);
	}

	// Get constants and do basic check on number of elements
	vector< map<string, float> > tvals;
	jcalib->Get("FDC/lorentz_deflections", tvals);
	if(tvals.size() != Nxpoints*Nzpoints){
		_DBG_<<"ERROR - GetLorentzDeflections() number of elements in calib DB"<<endl;
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
//----------------
// GetConstants
//----------------
int GetConstants(const char* namepath, int *Nvals, float* vals, mystr_t *strings)
{
	/// C-callable routine for accessing calibraion constants.
	/// The values specified by "namepath" will be read into the array
	/// "vals". The "vals" array should have enough memory allocated
	/// to hold *Nvals elements. If not, only the first *Nvals elements
	/// will be copied and a non-zero value returned. If the number
	/// of values in the database are less than *Nvals, then all values
	/// are copied, *Nvals is updated to reflect the number of valid
	/// elements in "vals", and a value of 0 is returned.
        /// Similar the variable names are stored in the array strings.
	
	if(!jcalib){
		_DBG_<<"ERROR - GetConstants() called when jcalib not set!"<<endl;
		_DBG_<<"ERROR - request for \""<<namepath<<"\""<<endl;
		_DBG_<<"ERROR - Exiting ..."<<endl;
		exit(-1);
	}

	map <string, float> detparms;
	jcalib->Get(namepath, detparms);

	if((int)detparms.size()<*Nvals)
	  *Nvals = (int)detparms.size();
	int i=0;
	for( map<string, float>::iterator ii=detparms.begin(); ii!=detparms.end(); ++ii){
	  if (i<*Nvals){
	    strcpy (strings[i].str, (*ii).first.c_str());
	    vals[i++] = (*ii).second;
	  }
	}
	return (int)detparms.size()>*Nvals; // return 0 if OK, 1 if not
}


//----------------
// GetArrayConstants
//----------------
int GetArrayConstants(const char* namepath, int *Nvals, float* vals, mystr_t *strings)
{
	/// C-callable routine for accessing calibration constants.
	/// The values specified by "namepath" will be read into the array
	/// "vals". The "vals" array should have enough memory allocated
	/// to hold *Nvals elements. If not, only the first *Nvals elements
	/// will be copied and a non-zero value returned. If the number
	/// of values in the database are less than *Nvals, then all values
	/// are copied, *Nvals is updated to reflect the number of valid
	/// elements in "vals", and a value of 0 is returned.
        /// Similar the variable names are stored in the array strings.
	
	if(!jcalib){
		_DBG_<<"ERROR - GetArrayConstants() called when jcalib not set!"<<endl;
		_DBG_<<"ERROR - request for \""<<namepath<<"\""<<endl;
		_DBG_<<"ERROR - Exiting ..."<<endl;
		exit(-1);
	}

	vector<map <string, float> >detparms;
	jcalib->Get(namepath, detparms);

	unsigned int i=0;
	int j=0;
	for (i=0;i<detparms.size();i++){
	  for( map<string, float>::iterator ii=detparms[i].begin(); ii!=detparms[i].end(); ++ii){
	    if (j<*Nvals){
	      strcpy (strings[j].str, (*ii).first.c_str());
	      vals[j] = (*ii).second;
	      j++;
	    }
	    else return 1;
	  }
	}
	*Nvals=j;

	return 0; // return 0 if OK, 1 if not
}

//------------------
// GetMD5Geom
//------------------
const char* GetMD5Geom(void)
{
	// Get the MD5 checksum of the geometry that will be
	// used for the simulation. This will retrieve the
	// geometry checksum from either what has been statically
	// linked in, or dynamically, whichever is being used.

	// This is a little odd since the string originates
	// in a FORTRAN routine.
	static char md5[256];
	memset(md5, 0, 256);
	md5geom_wrapper_(md5);
	
	md5[32] = 0; // truncate string at 32 characters (FORTRAN adds a space)
	
	return md5;
}

