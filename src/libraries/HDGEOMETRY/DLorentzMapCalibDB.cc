#include "DLorentzMapCalibDB.h"
#include "DLorentzDeflections.h"
#include <map>
#include <vector>
using namespace std;


//---------------------------------
// DLorentzMapCalibDB    (Constructor)
//---------------------------------
DLorentzMapCalibDB::DLorentzMapCalibDB(JApplication *japp)
{
  int runnumber = 1;
  jcalib = japp->GetJCalibration(runnumber);

  int Npoints = GetLorentzDeflections(runnumber); 
  if(Npoints==0){
    _DBG_<<"Error getting JCalibration object for Lorentz corrections!"<<
      endl;
    japp->Quit();
  }
}

//---------------------------------
// DLorentzMapCalibDB    (Constructor)
//---------------------------------
DLorentzMapCalibDB::DLorentzMapCalibDB(JCalibration *jcalib)
{
  this->jcalib = jcalib;
  if(GetLorentzDeflections()==0){
    _DBG_<<"Error getting JCalibration object for Lorentz corrections!"<<
      endl;
    exit(1);
  } 
}



unsigned int DLorentzMapCalibDB::GetLorentzDeflections(unsigned int runno){  
  /// Routine for accessing calibration constants adapted from code written by 
  /// David Lawrence.  
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
    return 0;
  }
  
  // Get constants and do basic check on number of elements
  vector< map<string, float> > tvals;
  jcalib->Get("FDC/lorentz_deflections", tvals);
  if(tvals.size() != LORENTZ_X_POINTS*LORENTZ_Z_POINTS){
    _DBG_<<"ERROR - GetLorentzDefelections() number of elements in calib DB"<<endl;
    _DBG_<<"ERROR - not the same as expected. DB="<<tvals.size()<<" expected"<<LORENTZ_X_POINTS*LORENTZ_Z_POINTS<<endl;
    _DBG_<<"ERROR - Exiting ..."<<endl;
    return 0;
  }
  
  // Notify user
  jout<<"Read "<<tvals.size()<<" values from FDC/lorentz_deflections in calibDB"<<endl;
  jout<<"   lorentz_deflections columns (alphabetical): ";
  map<string,float>::iterator iter;
  for(iter=tvals[0].begin(); iter!=tvals[0].end(); iter++)jout<<iter->first<<" ";
  jout<<endl;

  // Copy values into tables. We preserve the order since that is how it was
  // originally done in hitFDC.c
  for(unsigned int i=0; i<tvals.size(); i++){
    map<string, float> &row = tvals[i];
    unsigned int xindex = i/LORENTZ_Z_POINTS;
    unsigned int zindex = i%LORENTZ_Z_POINTS;
    lorentz_x[xindex] = row["x"];
    lorentz_z[zindex] = row["z"];
    lorentz_nx[xindex][zindex] = row["nx"];
    lorentz_nz[xindex][zindex] = row["nz"];
  }
  return tvals.size();
}
