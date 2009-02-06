#include "DMaterialMapCalibDB.h"

#include <map>
#include <vector>
using namespace std;


//---------------------------------
// DMaterialMapCalibDB    (Constructor)
//---------------------------------
DMaterialMapCalibDB::DMaterialMapCalibDB(JApplication *japp)
{
  int runnumber = 1;
  jcalib = japp->GetJCalibration(runnumber);

  int Npoints = GetMaterialMap(runnumber); 
  if(Npoints==0){
    _DBG_<<"Error getting JCalibration object for Material properties"<<
      endl;
    japp->Quit();
  }
}

//-----------------------------------
// GetMaterialMap - routine for retreiving material information from the 
// calibration database
//-----------------------------------
unsigned int DMaterialMapCalibDB::GetMaterialMap(unsigned int runno){
  // Make sure jcalib is set
  if(!jcalib){
    _DBG_<<"ERROR - GetMaterialMap() called when jcalib not set!"<<endl;
    _DBG_<<"ERROR - Exiting ..."<<endl;
    return 0;
  }
  
  // Get constants and do basic check on number of elements
  vector< map<string, float> > tvals;
  jcalib->Get("Material/radlen", tvals);
  if(tvals.size() != NUM_X_POINTS*NUM_Z_POINTS){
    _DBG_<<"ERROR - GetMaterialMap() number of elements in calib DB"<<endl;
    _DBG_<<"ERROR - not the same as expected. DB="<<tvals.size()<<", expected "<<NUM_X_POINTS*NUM_Z_POINTS<<endl;
    _DBG_<<"ERROR - Exiting ..."<<endl;
    return 0;
  }
  
  // Notify user
  cout<<"Read "<<tvals.size()<<" values from Material/radlen in calibDB"<<endl;
  cout<<"   radlen columns (alphabetical): ";
  map<string,float>::iterator iter;
  for(iter=tvals[0].begin(); iter!=tvals[0].end(); iter++)cout<<iter->first<<" ";
  cout<<endl;

  // Copy values into tables. 
  for(unsigned int i=0; i<tvals.size(); i++){
    map<string, float> &row = tvals[i];
    unsigned int zindex = i/NUM_X_POINTS;
    unsigned int rindex = i%NUM_X_POINTS;
    material_x[rindex] = row["x"];
    material_z[zindex] = row["z"];
    radlen[zindex][rindex] = row["radlen"];// radiation lengths of the material
    atomic_Z[zindex][rindex]=row["Z"]; // atomic number
    atomic_A[zindex][rindex]=row["A"]; // atomic weight
    density[zindex][rindex]=row["dens"];
  }
  return tvals.size();

}
