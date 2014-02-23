/*
 *  DBCALClump.cc
 *
 *  Created by Beni Zihlmann Tue Mar 12 2013
 *  version 0.1
 *
 */

#include "BCAL/DBCALClump.h"

#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMath.h>

DBCALClump::DBCALClump(vector<const DBCALHit*> U, vector<const DBCALHit*> D){

  HitsU = U;
  HitsD = D;

}

void DBCALClump::AnalyzeClump(){

  // now do some work on the Clumps up stream and down stream
  // create clump profiles and determine Energy Sums

  //const DBCALHit *BcalMatrixU[48*4][4]; // Up stream
  //const DBCALHit *BcalMatrixD[48*4][4]; // Down stream
  float EU[48*4];
  float ED[48*4];

  fillArrays(EU,ED);

  float EmaxU=0.;
  float EmaxD=0.;

  // find Maximum of energy in Clump array
  int MaxU=999;
  int MaxD=999;
  for (int i=0;i<48*4;i++){
    if ((EU[i]>0) && (ED[i]>0)) {
      if (EU[i]>EmaxU){
	EmaxU = EU[i];
	MaxU = i;
      }
      if (ED[i]>EmaxD){
	EmaxD = ED[i];
	MaxD = i;
      }
    }
  }

  //cout<<MaxU<<" / "<<MaxD<<endl;
  // use largest energy to define central bin
  if (MaxU>MaxD){
    MaxD = MaxU;
  } else {
    MaxU = MaxD;
  }

  resetProfiles();

  int UCounter=0;
  int DCounter=0;
  
  // set up a shower profile such that the maximum energy of the shower
  // ends up in bin 20;
  int CenterOffsetU = 20-MaxU; 
  int CenterOffsetD = 20-MaxD;

  // create shower profiles for up stream and down stream 
  int n=MaxU;
  // start from center upwards
  while (EU[n]>0.){
    int idx = n+CenterOffsetU;
    if (idx<0){
      idx += MaxU;
    }
    //cout<<"U> "<<n<<" "<<EU[n]<<" "<<idx<<endl;
    if ((idx>=60)||(idx<0)){
      break;
    }
    ProfileU[idx] = EU[n]; // limit shower upper limit
    UCounter++;
    // section to get average mean time and average time difference for sector
    int loc=999;
    for (unsigned int k=0;k<Sector.size();k++){
      if (Sector[k] == n){
	loc = k;
	break;
      }
    }
    if (loc<999){
      float Counter=0.;
      while (Sector[loc]==n){
	ProfileMT[idx] += MeanTime[loc];
	ProfileTD[idx] += DeltaTime[loc];
	loc++;
	Counter += 1.0;
      }
      ProfileMT[idx] /= Counter;
      ProfileTD[idx] /= Counter;
    }

    n++;
    if (n>=48*4){
      n -= 48*4;
    }
  }
  // now move from center down wards
  n = MaxU-1;
  if (n<0){
    n += 48*4;
  }
  while (EU[n]>0.){
    int idx = n+CenterOffsetU;
    //cout<<"U< "<<n<<" "<<EU[n]<<"  "<<idx<<endl;
    if (idx>20){
      idx -= (48*4);
    }
    if (idx<0){
      break;
    }
    ProfileU[idx] = EU[n];
    UCounter++;
    // section to get average mean time and average time difference for sector
    int loc=999;
    for (unsigned int k=0;k<Sector.size();k++){
      if (Sector[k] == n){
	loc = k;
	break;
      }
    }
    if (loc<999){
      float Counter=0.;
      while (Sector[loc]==n){
	ProfileMT[idx] += MeanTime[loc];
	ProfileTD[idx] += DeltaTime[loc];
	loc++;
	Counter += 1.0;
      }
      ProfileMT[idx] /= Counter;
      ProfileTD[idx] /= Counter;
    }

    n--;
    if (n<0){
      n += 48*4;
    }
  }
  n=MaxD;
  // start from center upwards
  while (ED[n]>0.){

    int idx = n+CenterOffsetD;
    //cout<<"D> "<<n<<" "<<ED[n]<<"  "<<idx<<endl;
    if (idx<0){
      idx += MaxD;
    }
    if ((idx>=60)||(idx<0)){
      break;
    }
    ProfileD[idx] = ED[n];
    DCounter++;
    n++;
    if (n>=48*4){
      n -= 48*4;
    }
  }
  // now move from center down wards
  n = MaxD-1;
  while (ED[n]>0.){

    int idx = n+CenterOffsetD;
    //cout<<"D< "<<n<<" "<<ED[n]<<"   "<<idx<<endl;
    if (idx>20){
      idx -= (48*4);
    }
    if (idx<0){
      break;
    }
    ProfileD[idx] = ED[n];
    DCounter++;
    n--;
    if (n<0){
      n += 48*4;
    }
  }

  // at this point we have a shower profile for the upstream and down stream Clump

  if (0){

    for (int k=0; k<60; k++){
      cout<<ProfileU[k]<<" ";
    }
    cout<<endl;
    for (int k=0; k<60; k++){
      cout<<ProfileD[k]<<" ";
    }
    cout<<endl;

  }


  // AT THIS POINT CODE SHOULD GO THAT TESTS IF THERE IS MORE THAN ONE
  // PARTICLE-SHOWER IN THE CLUMP!!!!!!
  // THIS IS NOT YET DONE IN VERSION 0.1 each Clump is just one shower

  double x[60];
  for (int k=0;k<60;k++){
    x[k] = (double) k;
  }
  int low = 0;
  while (ProfileU[low]==0.){
    low++;
  }
  int high = low;
  while (ProfileU[high]!=0){
    high++;
  }
  low--;
  high++;
  
  //double chi2U = 999.;
  double posU = 20.;
  if (UCounter>1){
    TGraph* grU = new TGraph(60,x,ProfileU);
    //grU->Fit("gaus","Q","r",(Double_t)low,(Double_t)high);
    
    grU->Fit("gaus","Q","r",15.,25.);
    TF1 *f1U = grU->GetFunction("gaus");

    if (f1U) {
      //chi2U = f1U->GetChisquare();  
      posU = f1U->GetParameter(1);
    } else {
      posU = 20.;
      //chi2U = 999.;
    }
    if (TMath::Abs(posU-20.)>2.){ // this shower has a complicated shape
      f1U->SetParameter(1,20.);
      f1U->FixParameter(2,1.0);
      //grU->Fit("gaus","Q","r",(Double_t)low,(Double_t)high);
      grU->Fit("gaus","Q","r",15.,25.);
      f1U = grU->GetFunction("gaus");
      posU = f1U->GetParameter(1);
      if (TMath::Abs(posU-20.)>2.){ // still no better result
	posU = 20.; // fix it.
      }
    }
  } else if (UCounter){
    posU = 20.;
    //chi2U = 999;
  } else {
    cout<<"Error no data to fit U: this should never happen!!!"<<endl;
    posU = 0;
    //chi2U=999999;
  }

  low = 0;
  while (ProfileU[low]==0.){
    low++;
  }
  high = low;
  while (ProfileU[high]!=0){
    high++;
  }
  low--;
  high++;
  
  //double chi2D = 999.;
  double posD = 20.;
  if (DCounter>1){
    TGraph* grD = new TGraph(60,x,ProfileD);
    //grD->Fit("gaus","Q","r",(Double_t)low,(Double_t)high);
    grD->Fit("gaus","Q","r",15.,25.);
    
    TF1 *f1D = grD->GetFunction("gaus");
    if (f1D){
      //chi2D = f1D->GetChisquare();  
      posD = f1D->GetParameter(1);
    } else {
      posD = 20.;
      //chi2D = 999.;
    }
    if (TMath::Abs(posD-20.)>2.){ // this shower has a complicated shape
      f1D->SetParameter(1,20.);
      f1D->FixParameter(2,1.0);
      //grD->Fit("gaus","Q","r",(Double_t)low,(Double_t)high);
      grD->Fit("gaus","Q","r",15.,25.);
      f1D = grD->GetFunction("gaus");
      posD = f1D->GetParameter(1);
      if (TMath::Abs(posD-20.)>2.){ // still no better result
	posD = 20.; // fix it.
      }
    }
  } else if (DCounter){
      posD = 20.;
      //chi2D = 999.;    
  } else {
    cout<<"Error no data to fit D: this should never happen!!!"<<endl;
    posU = 0;
    //chi2U=999999;
  }


  //cout<<chi2U<<"  "<<chi2D<<"         ";
  //cout<<posU<<"  "<<posD<<endl;

  //double UD_asymmetry = (posD-posU)/posD; // Note this quantity may be used for multi shower searches. 

  // The following assumes that the first BCAL module (0) starts  at 9 o'clock
  // looking down stream and the counting is goind clock wise.
  float phiD = (posD-20.+MaxD+0.5) * 2.*3.1415926/(48.*4.) ;
  float phiU = (posU-20.+MaxU+0.5) * 2.*3.1415926/(48.*4.) ;
  ClumpPhi.push_back((phiD+phiU)/2.);
  //cout<<phid*180./3.1415926<<endl;


  // calculate total energy of the shower from up and down stream clump energy
  // using the geometri mean approach
  double ESumU = 0.0;
  double ESumD = 0.0;
  for (unsigned int k=0;k<60;k++){
    ESumU += ProfileU[k];
    ESumD += ProfileD[k];
  }

  // ATTENTION HARD CODED VALUES 195/300 = BcalLength/2/ATTENUATIONLENGTH !!!!!!! 
  double E = TMath::Sqrt(ESumU*ESumD)/TMath::Exp(-195./300.);
  ClumpE.push_back(E);
  //  cout<<E<<endl;

  float MTaverage=0;
  double Position[4] = {0.,0.,0.,0.};
  double err[4] = {200.,200.,200.,200.};
  double ex[4] = {0.,0.,0.,0.};
  double Pcnt[4] = {0.,0.,0.,0.};
  for (unsigned int k=0;k<MeanTime.size();k++){
    MTaverage += MeanTime[k]; 
    int idx = Layer[k];
    Position[idx] += DeltaTime[k];
    Pcnt[idx] += 1.;
  }
  if (MeanTime.size()){
    MTaverage /= (float)MeanTime.size();
  } else{
    MTaverage = 0.;
  }
  ClumpMT.push_back(MTaverage);

  low=4;
  high=0;
  for (int k=0;k<4;k++){
    if (Pcnt[k]>0.){
      Position[k] /= Pcnt[k];
      err[k] = 2.5; // Hard coded error
      if (k<low){
	low = k;
      }
      if (k>high){
	high = k;
      }
    }
  }
  
  double Pos=0;
  double cnt = 0;
  // don't do a fit with only two or less points
  // fit approach seem to be too instable 
  // do at the moment do only average
  if ((high-low)<6){ 
    for (int k=0;k<4;k++){
      if (Pcnt[k]>0){
	Pos += Position[k];
	cnt += 1.;
      }
    }
    if (cnt>0){
      Pos /= cnt;
    }

    //cout<<"           "<<Pos<<endl;
    //for (int i=0;i<4;i++){
    //  cout<<Position[i]<<"   ";
    //}
    //cout<<endl;
    
  } else { // do a fit with 3 or more points
    double xx[4] = {1.,3.,8.,15.};
    TGraphErrors* grT = new TGraphErrors(4,xx,Position,ex,err);
    grT->Fit("pol1","Q","R",low,high);
    TF1 * lfit = grT->GetFunction("pol1");
    double offset = lfit->GetParameter(0);
    //double slope = lfit->GetParameter(1);
    //cout<<slope<<" * x + "<<offset<<endl;
    //for (int i=0;i<4;i++){
    //  cout<<Position[i]<<"   ";
    // }
    //cout<<endl;
    Pos = offset;
  }

  // adjust postion to be zero at the up stream end
  // and full length at the down stream end
  Pos = 195.-Pos; // hard code  half length of module
  ClumpPos.push_back(Pos);
  //cout<<Pos<<endl;

}

void DBCALClump::fillArrays(float* EU, float* ED){

  // reset array
  for (int i=0;i<48*4;i++){
    EU[i] = 0.;
    ED[i] = 0.;
  }

  for (unsigned int i=0;i<HitsU.size();i++){
    const DBCALHit *hit = HitsU[i];
    int idx = (hit->sector-1) + (hit->module-1)*4;
    EU[idx] += hit->E;
  }
  for (unsigned int i=0;i<HitsD.size();i++){
    const DBCALHit *hit = HitsD[i];
    int idx = (hit->sector-1) + (hit->module-1)*4;
    ED[idx] += hit->E; 
  }
}

void DBCALClump::resetProfiles(void){

  for (int k=0;k<60;k++){
    ProfileU[k] = 0.0;
    ProfileD[k] = 0.0;
    ProfileMT[k] = 0.0;
    ProfileTD[k] = 0.0;
  }

}

