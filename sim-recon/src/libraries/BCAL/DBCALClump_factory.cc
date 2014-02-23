/*
 *  DBCALClump.h
 *
 *  Created by Beni Zihlmann Tue Mar 12 2013
 *
 */


#include <JANA/JApplication.h>
using namespace jana;

#include <DBCALClump_factory.h>


//----------------
// init
//----------------
jerror_t DBCALClump_factory::init(void)
{
  return NOERROR;
}

//----------------
// init
//----------------
jerror_t DBCALClump_factory::brun(JEventLoop *loop)
{
  map<string, double> bcalparms;

  if ( !loop->GetCalib("BCAL/mc_parms", bcalparms)){
    cout<<"DBCALClump_factory: loading values from TOF data base"<<endl;
  } else {
    cout << "DBCALClumpo_factory: Error loading values from BCAL MC data base" <<endl;

    VELOCITY = 15.;  // set to some reasonable value

    return NOERROR;
  }

  VELOCITY = bcalparms["C_EFFECTIVE"];


  return NOERROR;
}



//----------------
// evnt
//----------------
jerror_t DBCALClump_factory::evnt(JEventLoop *loop, int eventnumber) {

  const DBCALHit *BcalMatrixU[48*4][4]; // Up stream
  const DBCALHit *BcalMatrixD[48*4][4]; // Down stream

  int MatrixUD[48*4][4];

  float EU[48*4];
  float ED[48*4];
  int EUi[48*4];
  int EDi[48*4];

  vector<const DBCALHit*> AllBcalHits;
  loop->Get(AllBcalHits);

  for (int i=0;i<48*4;i++){
    EU[i] = 0.;
    ED[i] = 0.;
    EUi[i] = 0;
    EDi[i] = 0;
    for (int m=0;m<4;m++){
      BcalMatrixU[i][m] = NULL;
      BcalMatrixD[i][m] = NULL;
      MatrixUD[i][m] = 0;
    }
  }

  // setup BCAL Matrix
  
  for (int i=0; i< (int)AllBcalHits.size();i++){
    const DBCALHit* hit = AllBcalHits[i];
    //cout<<"Module: "<<hit->module<<"   Sector: "<<hit->sector<<"   Layer: "<<hit->layer<<endl;
    int chan = hit->sector-1 + (hit->module-1)*4;
    if ((hit->t<200.) && (hit->t>0.))  { // HARD CODED
      if (hit->end == DBCALGeometry::kUpstream){
	BcalMatrixU[chan][hit->layer-1] = hit;
      } else{
	BcalMatrixD[chan][hit->layer-1] = hit;
      }
    }
  }

  // Build the up/down stream coincidence Matrix
  for (int i=0;i<48*4;i++){
    for (int m=0;m<4;m++){
      if ((BcalMatrixU[i][m]) && (BcalMatrixD[i][m])){
	float meant = (BcalMatrixU[i][m]->t + BcalMatrixD[i][m]->t - 390./16.75)/2.; // HARD CODED VALUES
	// this mean time has to be positive!!
	// and the cut could be made stronger.
	if (meant>0){
	  MatrixUD[i][m] = 1;
	} else {
	  BcalMatrixU[i][m] = NULL;
	  BcalMatrixD[i][m] = NULL;
	}
      }
    }
  }


  // first get rid of single-ended isolated hits
  for (int i=0;i<48*4;i++){
    int bef = i-1;
    int aft = i+1;
    if (bef<0){
      bef+=48*4;
    }
    if (aft>=48*4){
      aft -=48*4;
    }
    // layer 1;
    if (!MatrixUD[i][0]){
      if (BcalMatrixU[i][0]){
	if ((!BcalMatrixU[bef][0]) &&
	    (!BcalMatrixU[aft][0]) &&
	    (!BcalMatrixU[i][1]) &&
	    (!BcalMatrixU[aft][1]) &&
	    (!BcalMatrixU[aft][1])){
	  BcalMatrixU[i][0] = NULL;
	}
      }
      if (BcalMatrixD[i][0]){
	if ((!BcalMatrixD[bef][0]) &&
	    (!BcalMatrixD[aft][0]) &&
	    (!BcalMatrixD[i][1]) &&
	    (!BcalMatrixD[aft][1]) &&
	    (!BcalMatrixD[aft][1])){
	  BcalMatrixD[i][0] = NULL;
	}		
      }
    }
    // layer 2;
    if (!MatrixUD[i][1]){
      if (BcalMatrixU[i][1]){
	if ((!BcalMatrixU[bef][1]) &&
	    (!BcalMatrixU[aft][1]) &&
	    (!BcalMatrixU[i][2]) &&
	    (!BcalMatrixU[aft][2]) &&
	    (!BcalMatrixU[aft][2])){
	  BcalMatrixU[i][1] = NULL;
	}
      }
      if (BcalMatrixD[i][1]){
	if ((!BcalMatrixD[bef][1]) &&
	    (!BcalMatrixD[aft][1]) &&
	    (!BcalMatrixD[i][2]) &&
	    (!BcalMatrixD[aft][2]) &&
	    (!BcalMatrixD[aft][2])){
	  BcalMatrixD[i][1] = NULL;
	}		
      }
    }
    // layer 3;
    if (!MatrixUD[i][2]){
      if (BcalMatrixU[i][2]){
	if ((!BcalMatrixU[bef][2]) &&
	    (!BcalMatrixU[aft][2]) &&
	    (!BcalMatrixU[i][3]) &&
	    (!BcalMatrixU[aft][3]) &&
	    (!BcalMatrixU[aft][3])){
	  BcalMatrixU[i][2] = NULL;
	}
      }
      if (BcalMatrixD[i][2]){
	if ((!BcalMatrixD[bef][2]) &&
	    (!BcalMatrixD[aft][2]) &&
	    (!BcalMatrixD[i][3]) &&
	    (!BcalMatrixD[aft][3]) &&
	    (!BcalMatrixD[aft][3])){
	  BcalMatrixD[i][2] = NULL;
	}		
      }
    }
     // layer 4;
    if (!MatrixUD[i][3]){
      if (BcalMatrixU[i][3]){
	if ((!BcalMatrixU[bef][3]) &&
	    (!BcalMatrixU[aft][3]) &&
	    (!BcalMatrixU[i][2]) &&
	    (!BcalMatrixU[aft][2]) &&
	    (!BcalMatrixU[aft][2])){
	  BcalMatrixU[i][3] = NULL;
	}
      }
      if (BcalMatrixD[i][3]){
	if ((!BcalMatrixD[bef][3]) &&
	    (!BcalMatrixD[aft][3]) &&
	    (!BcalMatrixD[i][2]) &&
	    (!BcalMatrixD[aft][2]) &&
	    (!BcalMatrixD[aft][2])){
	  BcalMatrixD[i][3] = NULL;
	}		
      }
    }      
  }

  for (int i=0;i<48*4;i++){
    for (int m=0;m<4;m++){
      
      if (BcalMatrixU[i][m]){
	EU[i] += BcalMatrixU[i][m]->E;
	EUi[i] += 1;
      }
      if (BcalMatrixD[i][m]){
	ED[i] += BcalMatrixD[i][m]->E;
	EDi[i] += 1;
      }
      
    }
  }

  // clean up Matrix from obvious stray hits
  for (int i=1;i<48*4-1;i++) {   
    int a = EDi[i-1]+EUi[i-1];
    int b = EDi[i] + EUi[i];
    int c = EDi[i+1]+EUi[i+1];
    if ((b<2)&&(a==0)&&(c==0)){
      EDi[i] = 0;
      EUi[i] = 0;
      EU[i] = 0;
      ED[i] = 0;
    }
  }

  int CNT=1;
  while (CNT) {

    //cout<<"Find Next Clump"<<endl;
    int MaxU=999;
    int MaxD=999;
    
    float mU=0.;
    float mD=0.;
    // find largest energy deposition in the matrix
    for (int i=0;i<48*4;i++){
      if ((EU[i]>0.0)&&(ED[i]>0.0)){
	if (EU[i]>mU){
	  mU = EU[i];
	  MaxU = i;
	}
	if (ED[i]>mD){
	  mD = ED[i];
	  MaxD = i;
	}
      }
    }

    if (MaxU==999 || MaxD==999){
      return NOERROR;  
    }

    if (mU>mD){
      MaxD = MaxU;
    } else {
      MaxU = MaxD;
    }

    //cout<<MaxD<<" / "<<MaxU<<endl;

    int mod = MaxD;
    float Mt = 0;
    float cnt = 0;
    float EmaxU = 0.;
    float EmaxD = 0.;
    int idxMaxU[2] = {999,999};
    int idxMaxD[2] = {999,999};

    // loop over all layers in the seed sector
    // and its neighouring sectors and use all
    // instances with hits on both ends to determine 
    // the mean time    
    for (int k=mod-1; k<mod+2; k++) { 
      int idx = k;
      if (k<0){
	idx += 48*4;
      } else if (k>=4*48){
	idx -= 48*4;
      }
      for (int n=0;n<3;n++){
	const DBCALHit* hitD = BcalMatrixD[idx][n];
	const DBCALHit* hitU = BcalMatrixU[idx][n];
	if (hitD && hitU){
	  if (hitD->E > EmaxD){
	    EmaxD = hitD->E; 
	    idxMaxD[0] = idx;  // index for seed
	    idxMaxD[1] = n;
	  }
	  if (hitU->E > EmaxU){
	    EmaxU = hitU->E; 
	    idxMaxU[0] = idx;  // index for seed
	    idxMaxU[1] = n;
	  }
	  float mt = hitD->t + hitU->t;
	  Mt += mt;
	  cnt +=1.;
	}
      }
    }
    if (cnt>0){
      Mt /= cnt;  // need improve to use Energy weight 
    }

    //cout<<Mt<<"  "<<cnt<<endl;

    if (EmaxU>EmaxD){
      idxMaxD[0] = idxMaxU[0];
      idxMaxD[1] = idxMaxU[1];
    } else {
      idxMaxU[0] = idxMaxD[0];
      idxMaxU[1] = idxMaxD[1];
    }

    // Note: Mt is the flight time of the patricle creating this shower!
    //       This means Mt is the same for up and down stream!.
    //
    // Now collect all hits related to the seed hit for up and down stream
    // separately. Only hits that appear in time (+/- 2ns) with the seed hit 
    // are added to the lists for layers 1,2 and 3. Layer 4 is always added.
    // This may be changed if the timing information from the FADC is good enough.
    // All hits in the list regardless if they contribut or not are then removed 
    // from the list of hits so they can not be used twice.

    CNT = (int)cnt;

    if (cnt>0) {
      vector <const DBCALHit*> HitListU;
      vector <const DBCALHit*> HitListD;
      vector <float> MeanTime;
      vector <float> DeltaTime;
      vector <int> sector;
      vector <int> layer;
      
      float seedTime = BcalMatrixU[idxMaxU[0]][idxMaxU[1]]->t;
      int Idx = MaxU;
      // collect hits for up stream shower
      // first move upwards in sectors from the seed
      while(EUi[Idx]) {
	//cout<<"EUi> "<<Idx<<"  "<<EUi[Idx]<<endl;
	for (int i=0;i<4;i++){
	  const DBCALHit* hitU = BcalMatrixU[Idx][i];
	  const DBCALHit* hitD = BcalMatrixD[Idx][i];
	  if ((hitU)&&(hitD)){
	    float mt =  hitD->t + hitU->t;
	    //mt -= (390./VELOCITY);
	    mt -= (390./16.75); // subract detector internal path HARD CODED VALUES!!!!!!
	    mt /= 2.;
	    MeanTime.push_back(mt);

	    // Note below D-U means positive dt is upstream
	    // from bcal center
	    float dt = hitD->t - hitU->t;
	    //dt = dt/2./VELOCITY;
	    dt = dt/2.*16.75; // convert to cm HARD CODED VALUE!!!!!
	    DeltaTime.push_back(dt);
	    sector.push_back(Idx);
	    layer.push_back(i);
	  }

	  if (hitU){
	    //if ((hitU->t - seedTime)<2.){
	    HitListU.push_back(hitU);
	    //}
	    BcalMatrixU[Idx][i] = NULL; 
	  }
	}
	EUi[Idx] = 0;
	EU[Idx] = 0.;
	Idx++;
	if (Idx>=4*48){
	  Idx -= 4*48;
	}
      }
      // now move downwards from the seed center
      Idx = MaxU-1;
      if (Idx<0){
	Idx += 4*48;
      } 
      while(EUi[Idx]) {
	//cout<<"EUi< "<<Idx<<"  "<<EUi[Idx]<<endl;
	for (int i=0;i<4;i++){
	  const DBCALHit* hitU = BcalMatrixU[Idx][i];
	  const DBCALHit* hitD = BcalMatrixD[Idx][i];
	  if ((hitU)&&(hitD)){
	    float mt =  hitD->t + hitU->t;
	    //mt -= (390./VELOCITY);
	    mt -= (390./16.5); // subract detector internal drift
	    mt /= 2.;
	    MeanTime.push_back(mt);
	    float dt = hitD->t - hitU->t;
	    //dt = dt/2./VELOCITY;
	    dt = dt/2.*16.5;   // convert to cm
	    DeltaTime.push_back(dt);
	    sector.push_back(Idx);
	    layer.push_back(i);
	  }

	  if (hitU){
	    //if ((hitU->t - seedTime)<2.){
	    HitListU.push_back(hitU);
	    //}
	    BcalMatrixU[Idx][i] = NULL; 
	  }
	}

	EUi[Idx] = 0;
	EU[Idx] = 0.;
	Idx--;
	if (Idx<0){
	  Idx += 4*48;
	}
      }
      
      seedTime = BcalMatrixD[idxMaxD[0]][idxMaxD[1]]->t;
      // collect hits for down stream shower
      // first move upwards in sectors
      Idx = MaxD;
      while(EDi[Idx]) {
	//cout<<"EDi> "<<Idx<<"  "<<EUi[Idx]<<endl;
	for (int i=0;i<4;i++){
	  const DBCALHit* hit = BcalMatrixD[Idx][i];
	  if (hit){
	    //if ((hit->t - seedTime)<2.){
	    HitListD.push_back(hit);
	    //}
	    BcalMatrixD[Idx][i] = NULL; 
	  }
	}
	EDi[Idx] = 0;
	ED[Idx] = 0.;
	Idx++;
	if (Idx>=4*48){
	  Idx -= 4*48;
	}
      }
      // now move downwrads from center
      Idx = MaxD-1;
      if (Idx<0){
	Idx += 4*48;
      }
      while(EDi[Idx]) {
	//cout<<"EDi< "<<Idx<<"  "<<EUi[Idx]<<endl;
	for (int i=0;i<4;i++){
	  const DBCALHit* hit = BcalMatrixD[Idx][i];
	  if (hit){
	    //if ((hit->t - seedTime)<2.){
	    HitListD.push_back(hit);
	    //}
	    BcalMatrixD[Idx][i] = NULL; 
	  }
	}
	EDi[Idx] = 0;
	ED[Idx] = 0.;
	Idx--;
	if (Idx<0){
	  Idx += 4*48;
	}
      }

      if ( (HitListU.size()>0) && (HitListD.size()>0) ) {
	//cout<<"Make New Clump:"<<cnt<<endl;
	DBCALClump* myClump = new DBCALClump(HitListU, HitListD);
	myClump->MeanTime = MeanTime;
	myClump->DeltaTime = DeltaTime;
	myClump->Sector = sector;
	myClump->Layer = layer;
	myClump->AnalyzeClump();
	
	int oK = 1;
	if ((HitListU.size()==1) || (HitListD.size()==1)){
	  if (myClump->ClumpE[0]<50.) { // HARD CODED
	    oK = 0;
	  }
	}
	if (oK){
	  _data.push_back(myClump);
	} else {
	  delete myClump;
	}
      } else {
	cout<<"Error no hits in this Clump!!!! Event: "<<eventnumber<<"      cnt= "<<cnt <<endl;
      }
    }
        
  }

  return NOERROR;
}

