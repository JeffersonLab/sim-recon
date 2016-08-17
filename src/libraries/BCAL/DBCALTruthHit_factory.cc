//most of the code was orginally written by Matt Shepherd in
//DBCALCluster_factory.cc

#include <vector>
using namespace std;

#include <JANA/JApplication.h>
using namespace jana;

#include "BCAL/DBCALTruthHit_factory.h"
#include "BCAL/DBCALTruthCell.h"
#include "BCAL/DBCALGeometry.h"

#include "DANA/DApplication.h"

#include "units.h"

/// temp
//static const int BCAL_NUM_MODULES  = 48;
//static const int BCAL_NUM_LAYERS   =  4;
//static const int BCAL_NUM_SECTORS  =  4;


//----------------
// brun
//----------------
jerror_t DBCALTruthHit_factory::brun(JEventLoop *loop, int32_t runnumber) {
  // Only print messages for one thread whenever run number changes
  static pthread_mutex_t print_mutex = PTHREAD_MUTEX_INITIALIZER;
  static set<int> runs_announced;
  pthread_mutex_lock(&print_mutex);
  bool print_messages = false;
  if(runs_announced.find(runnumber) == runs_announced.end()){
    print_messages = true;
    runs_announced.insert(runnumber);
  }
  pthread_mutex_unlock(&print_mutex);

  DApplication* app = dynamic_cast<DApplication*>(loop->GetJApplication());
  DGeometry* geom = app->GetDGeometry(runnumber);
  //geom->GetTargetZ(m_z_target_center);

  if(print_messages) jout << "in DBCALTruthHit_factory, loading constants ..." << endl;

  // load attenuation correction parameters 
  attenuation_parameters.clear();

  vector< vector<double> > attenuation_parameters_temp;
  loop->GetCalib("/BCAL/attenuation_parameters", attenuation_parameters_temp);

  // avoid potential crash ...
  for (unsigned int i = 0; i < attenuation_parameters_temp.size(); i++){
      attenuation_parameters.push_back(attenuation_parameters_temp.at(i));
  }

   if (PRINTCALIBRATION) {
       int channel = 0;
       for (int module=1; module<=BCAL_NUM_MODULES; module++) {
           for (int layer=1; layer<=BCAL_NUM_LAYERS; layer++) {
               for (int sector=1; sector<=BCAL_NUM_SECTORS; sector++) {
                   printf("%2i  %2i  %2i %12.4f %12.4f %12.4f\n",
                          module,layer,sector,
                          attenuation_parameters[channel][0],
                          attenuation_parameters[channel][1],
                          attenuation_parameters[channel][2]);
               }
               channel++;
           }
       }
   }


  // load effective velocities
  effective_velocities.clear();

  // Passing in the member vector directly was sometimes causing a crash...
  vector <double> effective_velocities_temp;
  loop->GetCalib("/BCAL/effective_velocities", effective_velocities_temp);

  for (unsigned int i = 0; i < effective_velocities_temp.size(); i++){
    effective_velocities.push_back(effective_velocities_temp.at(i));
  }

  // load track parameters (the parameters of the quadratic fit in histograms of z_track = f(tUp - tDown)  )
  track_parameters.clear();

  vector< vector<double> > track_parameters_temp;
  loop->GetCalib("/BCAL/z_track_parms", track_parameters_temp);

  // Passing in the member vector directly was sometimes causing a crash...
  for (unsigned int i = 0; i < track_parameters_temp.size(); i++){
      track_parameters.push_back(track_parameters_temp.at(i));
  }

  return NOERROR;
}

//----------------
// evnt
//----------------
jerror_t DBCALTruthHit_factory::evnt(JEventLoop *loop, uint64_t eventnumber) {

  vector<const DBCALTruthCell*> cell;
  loop->Get(cell);
  if (cell.size() <= 0) return NOERROR;

  vector<pair < int ,double > > E_up;
  vector<pair < int, double > > E_down;
  E_up.clear();
  E_down.clear();

  double summed_E_up[768];
  double summed_E_down[768];

  double atten_E_up = 0.;
  double atten_E_down = 0;

  for(unsigned int i = 0 ; i < 768; i++){
    summed_E_up[i]=0;
    summed_E_down[i]=0;
  }

  double summed_E = 0.;
  int chan_indx = 0;
  int fADC_layer = 0;  
 
  int table_id = GetCalibIndex(1, 1, 1 );
  float fibLen = DBCALGeometry::GetBCAL_length();
  double attenuation_length = DBCALGeometry::ATTEN_LENGTH;
  double attenuation_L1=-1., attenuation_L2=-1.;
  double cEff = GetEffectiveVelocity(table_id);
  double track_p0 = -1.0; // will be updated from GetTrackParameters (dimensions: cm)
  double track_p1 = -1.0; // will be updated from GetTrackParameters (dimensions: cm/ns)
  double track_p2 = -100.0; // will be updated from GetTrackParameters (dimensions: cm/ns^2)

  for( vector< const DBCALTruthCell*>::const_iterator celltr = cell.begin();
       celltr != cell.end();
       ++celltr ){

       const DBCALTruthCell& cells = (**celltr);
   
       GetAttenuationParameters(table_id, attenuation_length, attenuation_L1, attenuation_L2);
       GetTrackParameters(table_id, track_p0, track_p1, track_p2);

       if( (**celltr).layer == 1) fADC_layer = 1;
       if( (**celltr).layer > 1 && (**celltr).layer < 4) fADC_layer = 2;
       if( (**celltr).layer > 3 && (**celltr).layer < 7) fADC_layer = 3;
       if( (**celltr).layer > 6 && (**celltr).layer < 11) fADC_layer = 4;

       int chan_num = 16*((**celltr).module-1) + 4*(fADC_layer-1) + ((**celltr).sector-1);

       double dUp = fibLen/2. + (**celltr).zLocal;
       double dDown = fibLen/2. - (**celltr).zLocal;

       double time_up = dUp/cEff + (**celltr).t;
       double time_down = dDown/cEff + (**celltr).t;

       atten_E_up = (**celltr).E*exp(-dUp/attenuation_length);
       atten_E_down = (**celltr).E*exp(-dDown/attenuation_length);

       E_up.push_back(make_pair(chan_num,atten_E_up));
       E_down.push_back(make_pair(chan_num,atten_E_down));
    }

    for(unsigned int i = 0 ; i < E_up.size(); i++){
      summed_E_up[E_up[i].first] += E_up[i].second;
      summed_E_down[E_down[i].first] += E_down[i].second;
    }
	
    for(int j = 0 ; j < 768 ; ++j){
       for(int i = 0 ; i < 2 ; ++i){
         GetAttenuationParameters(table_id, attenuation_length, attenuation_L1, attenuation_L2);
         GetTrackParameters(table_id, track_p0, track_p1, track_p2);
	 int end_temp = i;
         if(end_temp == 0){
		summed_E = summed_E_up[j];
	 }
	 if(end_temp == 1){
		summed_E = summed_E_down[j];
	 }
         chan_indx = j;
         if((summed_E_up[j] > 0 && end_temp==0)|| (summed_E_down[j]>0 && end_temp==1)){
	 DBCALTruthHit *truthhit = new DBCALTruthHit(attenuation_length,cEff,track_p0,track_p1,track_p2,end_temp,chan_indx,summed_E);

          _data.push_back(truthhit);
	 }
       }
     }
  
  return NOERROR;
}

bool DBCALTruthHit_factory::GetAttenuationParameters(int id, double &attenuation_length, 
						    double &attenuation_L1, double &attenuation_L2)
{
	vector<double> &parms = attenuation_parameters.at(id);

	attenuation_length = parms[0];
	attenuation_L1 = parms[1];
	attenuation_L2 = parms[2];

	return true;
}

double DBCALTruthHit_factory::GetEffectiveVelocity(int id)
{
	return effective_velocities.at(id);
}

bool DBCALTruthHit_factory::GetTrackParameters(int id, double &track_p0, 
						    double &track_p1, double &track_p2)
{
	vector<double> &z_parms = track_parameters.at(id);

	if(!z_parms.empty()){
		track_p0 = z_parms[0];
		track_p1 = z_parms[1];
		track_p2 = z_parms[2];
		return true;
	}
  	else{
  		jerr<<"Failed to retrieve the z_track parameters from CCDB!!!" << endl;
  		exit(-1);
  	}
}
