//most of the code was orginally written by Matt Shepherd in
//DBCALCluster_factory.cc

#include <vector>
using namespace std;

#include <JANA/JApplication.h>
using namespace jana;

#include "BCAL/DBCALPoint_factory.h"
#include "BCAL/DBCALHit.h"

#include "DANA/DApplication.h"

#include "units.h"

/// temp
static const int BCAL_NUM_MODULES  = 48;
static const int BCAL_NUM_LAYERS   =  4;
static const int BCAL_NUM_SECTORS  =  4;


//----------------
// brun
//----------------
jerror_t DBCALPoint_factory::brun(JEventLoop *loop, int runnumber) {
  DApplication* app = dynamic_cast<DApplication*>(loop->GetJApplication());
  DGeometry* geom = app->GetDGeometry(runnumber);
  geom->GetTargetZ(m_z_target_center);

  cout << "in DBCALPoint_factory, loading constants ..." << endl;

  // load attenuation correction parameters 
  vector< vector<double> > in_atten_parameters;
  loop->GetCalib("/BCAL/attenuation_parameters", in_atten_parameters);

  // if (PRINTCALIBRATION) {
  //   jout << "DBCALPoint_factory::brun >>Printing " << endl;
  // }
  attenuation_parameters.clear();
  int channel = 0;
  for (int module=1; module<=BCAL_NUM_MODULES; module++) {
	  for (int layer=1; layer<=BCAL_NUM_LAYERS; layer++) {
		  for (int sector=1; sector<=BCAL_NUM_SECTORS; sector++) {
			  //int cell_id = GetCalibIndex(module,layer,sector);

			  vector<double> new_params(3,0.);
			  //attenuation_parameters[cell_id][0] = in_atten_parameters[channel][0];
			  //attenuation_parameters[cell_id][1] = in_atten_parameters[channel][1];
			  //attenuation_parameters[cell_id][2] = in_atten_parameters[channel][2];
			  // hack to workaround odd CCDB behavior
			  //attenuation_parameters[cell_id][0] = in_atten_parameters[channel][1];
			  //attenuation_parameters[cell_id][1] = in_atten_parameters[channel][2];
			  //attenuation_parameters[cell_id][2] = in_atten_parameters[channel][0];
			  
			  new_params[0] = in_atten_parameters[channel][1];
			  new_params[1] = in_atten_parameters[channel][2];
			  new_params[2] = in_atten_parameters[channel][0];
			  attenuation_parameters.push_back( new_params );

			  if (PRINTCALIBRATION) {
			    printf("%2i  %2i  %2i %12.4f %12.4f %12.4f\n",
				   module,layer,sector,
				   attenuation_parameters[channel][0],
				   attenuation_parameters[channel][1],
				   attenuation_parameters[channel][2]);
			  }
/*
			  cerr << " loaded " << cell_id << " = " << attenuation_parameters[cell_id][0] << ", "
			       << attenuation_parameters[cell_id][1] << ", "
			       << attenuation_parameters[cell_id][2] << endl;
*/

			  channel++;
		  }
	  }
  }


  // load effective velocities
  effective_velocities.clear();
  loop->GetCalib("/BCAL/effective_velocities", effective_velocities);

    return NOERROR;
}

//----------------
// evnt
//----------------
jerror_t DBCALPoint_factory::evnt(JEventLoop *loop, int eventnumber) {

  vector<const DBCALUnifiedHit*> hits;
  loop->Get(hits);
  if (hits.size() <= 0) return NOERROR;

  // first arrange the list of hits so they are grouped by cell
  map< int, cellHits > cellHitMap;
  for( vector< const DBCALUnifiedHit* >::const_iterator hitPtr = hits.begin();
       hitPtr != hits.end();
       ++hitPtr ){
    
    const DBCALUnifiedHit& hit = (**hitPtr);
    
    int id = DBCALGeometry::cellId( hit.module, hit.layer, hit.sector );

    // Add hit to appropriate list for this cell
    if(hit.end == DBCALGeometry::kUpstream){
      cellHitMap[id].uphits.push_back( *hitPtr );
    }else{
      cellHitMap[id].dnhits.push_back( *hitPtr );
    }
  }

  // now go through this list and group hits into BCAL points
  // this combines information from both ends
  for( map< int, cellHits >::const_iterator mapItr = cellHitMap.begin();
       mapItr != cellHitMap.end();
       ++mapItr ){
    
    const vector<const DBCALUnifiedHit*> &uphits = mapItr->second.uphits;
    const vector<const DBCALUnifiedHit*> &dnhits = mapItr->second.dnhits;

    //require double-ended hits
    if(uphits.size()==0 || dnhits.size()==0) continue;

    // Each SiPM sum can have multiple hits, some caused purely by
    // dark hits. A more sophisticated algorithm may be needed here
    // to decipher the multi-hit events. For now, we just take the
    // most energetic hit from each end. (Single ended hits are
    // ignored.

    const DBCALUnifiedHit *uphit=uphits[0];
    const DBCALUnifiedHit *dnhit=dnhits[0];

    for(unsigned int i=1; i<uphits.size(); i++){
      if(uphits[i]->E > uphit->E) uphit = uphits[i];
    }

    for(unsigned int i=1; i<dnhits.size(); i++){
      if(dnhits[i]->E > dnhit->E) dnhit = dnhits[i];
    }

    // first check that the hits don't have absurd timing information

    //int id = DBCALGeometry::cellId( uphit->module, uphit->layer, uphit->sector );  // key the cell identification off of the upstream cell
    int table_id = GetCalibIndex( uphit->module, uphit->layer, uphit->sector );  // key the cell identification off of the upstream cell

    float fibLen = DBCALGeometry::BCALFIBERLENGTH;
    //float cEff = DBCALGeometry::C_EFFECTIVE;    
    float cEff = GetEffectiveVelocity(table_id);

    // get the position with respect to the center of the module -- positive
    // z in the downstream direction
    double zLocal = 0.5 * cEff * ( uphit->t - dnhit->t );

    // if the timing information indicates that the z position is more than 60 cm outside the BCAL, likely the hit is contamined by noise or entirely noise, skip this cell
    double tol = 60*k_cm;

    if (zLocal > (0.5*fibLen + tol) || zLocal < (-0.5*fibLen - tol)) continue;

    // pass attenuation length parameters to the DBCALPoint constructor, since
    // many of the calculations are implemented there
    double attenuation_length = DBCALGeometry::ATTEN_LENGTH;
    double attenuation_L1=-1., attenuation_L2=-1.;  // these parameters are ignored for now
    GetAttenuationParameters(table_id, attenuation_length, attenuation_L1, attenuation_L2);
    // if (GetAttenuationParameters(id, attenuation_length, attenuation_L1, attenuation_L2)) {
    //   printf("got new att length %f\n",attenuation_length);
    // } else {
    //   printf("default att length %f\n",attenuation_length);
    // }

    DBCALPoint *point = new DBCALPoint(*uphit,*dnhit,m_z_target_center,attenuation_length);

    point->AddAssociatedObject(uphit);
    point->AddAssociatedObject(dnhit);

    _data.push_back(point);
  }

  //Possibly we should also construct points from single-ended hits here.
  //The code for this is currently (commented out) in
  //DBCALCluster_factory.cc

  return NOERROR;
}

bool DBCALPoint_factory::GetAttenuationParameters(int id, double &attenuation_length, 
						    double &attenuation_L1, double &attenuation_L2)
{
	vector<double> &parms = attenuation_parameters.at(id);

	attenuation_length = parms[0];
	attenuation_L1 = parms[1];
	attenuation_L2 = parms[2];

	return true;
}

double DBCALPoint_factory::GetEffectiveVelocity(int id)
{
	return effective_velocities.at(id);
}
