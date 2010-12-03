// $Id$
//
//    File: DTOFMCResponse_factory.cc
// Created: Mon Aug 15 11:33:45 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8 i686)
//
// changes:  Thu Jun 21 12:16:55 EDT 2007 add digitization for energy and 
//                                        smearing for timeing
//

#include <iostream>
using namespace std;

#include "DTOFMCResponse_factory.h"
#include "DTOFHitRaw.h"
#include "DTOFGeometry.h"

// root specific stuff
#include "TRandom.h"  // random generators

#define TWO_HIT_RESOL   25.
#define MAX_HITS        100
#define THRESH_MEV      0.8
#define ECHARGE         1.602e-7 // electric charge in C

TRandom Ran;
float nph[2], R, A1, A2;
float charge[2], energy[2],tof[2];
float dist,loca[2],sigm[2];
int adc[2],t1,tdc[2];


//------------------
// evnt
//------------------
jerror_t DTOFMCResponse_factory::evnt(JEventLoop *loop, int eventnumber)
{

  vector<const DTOFHitRaw*> hddmhits;
  eventLoop->Get(hddmhits);

  vector<const DTOFGeometry*> tofGeomVect;
  eventLoop->Get(tofGeomVect);
 
  const DTOFGeometry& tofGeom = (*(tofGeomVect[0]));

  //  cout<<"TOF Response Factory: Number of hits in FTOF: "
  //    <<hddmhits.size()<< endl;

  for (unsigned int i = 0; i < hddmhits.size(); i++){

    const DTOFHitRaw *hddmhit = hddmhits[i];
    DTOFMCResponse *response = new DTOFMCResponse;

    // do any run-dependent smearing here

    // calculate total number of generated photons in 4pi

    nph[0] = hddmhit->dE_north*1000.* (PHOTONS_PERMEV);
    nph[1] = hddmhit->dE_south*1000.* (PHOTONS_PERMEV);

    //    R = tofGeom.LONGBARLENGTH/2.;

    tof[0] = hddmhit->t_north;  // TOF at PMT
    tof[1] = hddmhit->t_south; 

    if (hddmhit->plane==0){  // 0 is vertical plane, 1 is horizontal plane
      dist = hddmhit->y;
    } else {
      dist = hddmhit->x;
    }
    for (int i=0;i<2;i++){
      loca[i] = pow(-1.,double(i))*dist;
      // position dependent timing resolution see NIM A525 (2004)183
      sigm[i] = loca[i]*(-0.0008926) + TOF_CENT_TRES; 
      sigm[i] = (double)pow(10.,double(sigm[i]));
      if (tof[i]>0){
	tof[i] += Ran.Gaus(0.,sigm[i]); // add time resolution of PMT
	tdc[i] = int(tof[i]/TDC_RES);
      } else{
	tdc[i] = 0;
      }
    }

    double esmear =  Ran.Gaus(0.,sqrt(nph[0]));
    for (int i=0;i<2;i++){

      energy[i] = nph[i] += esmear;   // smear the number of generated photons
      loca[i] -= HALFPADDLE ;                            // distance to PMT

      // fraction of 4pi seen by PMT directly
      // surace area ratio between PMT and photon emitted in shpere that reach the PMT (total interal reflection)
      A1 = PMT_SURFACE/4./loca[i]/loca[i]/3.1415926;   
                                                         
      // fraction of 4pi below critical angle
      A2 = (1.-cos(THETA_MAX/180.*3.1415926))/2. - A1;
      nph[i] *= (A1+A2*REFLECT_EFF)*exp(loca[i]/ATTEN_LENGTH);  // number of photons hitting the photo cathode
      nph[i] *= PHE_EFF;             // number of photo electrons
      charge[i] = nph[i]*PMT_GAIN*ECHARGE; // charge in pC at PMT Anode
      adc[i] = (int)(charge[i]/ADC_RES<1024.0 ? charge[i]/ADC_RES:1024.0);
      //adc[i] = (int)fmin(charge[i]/ADC_RES,1024.); // Replaced with above line because fmin is not available on all compilers 6/24/07 D.L.


      energy[i] = energy[i]/(PHOTONS_PERMEV)/1000.; // energy in GeV
    }

    response->id          = hddmhit->id;
    response->orientation = hddmhit->plane;
    response->bar         = hddmhit->bar;
    response->y           = tofGeom.bar2y( hddmhit->bar, hddmhit->plane );
    response->t_north     = tof[0];
    response->t_south     = tof[1];
    response->ADC_north   = adc[0];
    response->ADC_south   = adc[1];
    response->TDC_north   = tdc[0];
    response->TDC_south   = tdc[1];
    response->E_north     = energy[0];
    response->E_south     = energy[1];
    response->ptype       = hddmhit->ptype;
    _data.push_back(response);

  }

  return NOERROR;
}


//------------------
// init
//------------------
jerror_t DTOFMCResponse_factory::brun(JEventLoop *loop, int runnumber)
{

  map<string, double> tofparms;
 
  if ( !loop->GetCalib("TOF/tof_parms", tofparms)){
    cout<<"DTOFMCResponse_factory: loading values from TOF data base"<<endl;
  } else {
    cout << "DTOFMCResponse_factory: Error loading values from TOF data base" <<endl;
  }

  ATTEN_LENGTH   =    tofparms["TOF_ATTEN_LENGTH"]; 
  C_EFFECTIVE    =    tofparms["TOF_C_EFFECTIVE"];
  PHOTONS_PERMEV =    tofparms["TOF_PHOTONS_PERMEV"];
  THETA_MAX      =    tofparms["TOF_THETA_MAX"];
  PMT_SURFACE    =    tofparms["TOF_PMT_SURFACE"];
  REFLECT_EFF    =    tofparms["TOF_REFLECT_EFF"];
  PHE_EFF        =    tofparms["TOF_PHE_EFF"];
  PMT_GAIN       =    tofparms["TOF_PMT_GAIN_MC"];
  ADC_RES        =    tofparms["TOF_ADC_RES_MC"];
  TOF_CENT_TRES  =    tofparms["TOF_CENT_TRES"];
  TDC_RES        =    tofparms["TOF_TDC_RES_MC"];
  HALFPADDLE     =    tofparms["TOF_HALFPADDLE"];

  return NOERROR;

}
