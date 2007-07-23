// $Id$
//
//    File: DTOFMCResponse_factory.cc
// Created: Mon Aug 15 11:33:45 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8 i686)
//
// changes:  Thu Jun 21 12:16:55 EDT 2007 add digitization for energy and 
//                                        smearing for timeing
//

#include "DTOFMCResponse_factory.h"
#include "DHDDMTOFHit.h"
#include "DTOFGeometry.h"

// root specific stuff
#include "TRandom.h"  // random generators

#define ATTEN_LENGTH    150
#define C_EFFECTIVE     15
#define TWO_HIT_RESOL   25.
#define MAX_HITS        100
#define THRESH_MEV      0.8
#define PHOTONS_PERMEV  10000. // scintillation photon generated per MeV energy deposition 
#define THETA_MAX       50.74 // total internal reflection
#define PMT_SURFACE     19.63 // PMT surface in cm^2
#define REFLECT_EFF     0.5   // reflection efficiency of light 
#define PHE_EFF         0.2   // efficiency to create photo elecrons
#define PMT_GAIN        40000000. // PMT gain factor 4*10^7
#define ECHARGE         1.602e-7 // electric charge in C
#define ADC_RES         50.      // adc resolution pC/count
#define TOF_CENT_TRES   -0.698970004 // time resolution of tof paddel for hit in the center log10(0.2ns)
#define TDC_MC_RES         0.06 //TDC resolution in ns

TRandom Ran;
float nph[2], R, A1, A2;
float charge[2], energy[2],tof[2];
float dist,loca[2],sigm[2];
int adc[2],t1,tdc[2];

void debugfuncMCResponse(){
  //do nothing
}


//------------------
// evnt
//------------------
jerror_t DTOFMCResponse_factory::evnt(JEventLoop *loop, int eventnumber)
{

  vector<const DHDDMTOFHit*> hddmhits;
  eventLoop->Get(hddmhits);

  vector<const DTOFGeometry*> tofGeomVect;
  eventLoop->Get(tofGeomVect);
 
  const DTOFGeometry& tofGeom = (*(tofGeomVect[0]));

  //  cout<<"TOF Response Factory: Number of hits in FTOF: "
  //    <<hddmhits.size()<< endl;

  for (unsigned int i = 0; i < hddmhits.size(); i++){

    const DHDDMTOFHit *hddmhit = hddmhits[i];
    DTOFMCResponse *response = new DTOFMCResponse;

    // do any run-dependent smearing here

    // calculate total number of generated photons in 4pi


    debugfuncMCResponse();

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
	tdc[i] = int(tof[i]/TDC_MC_RES);
      } else{
	tdc[i] = 0;
      }
    }

    for (int i=0;i<2;i++){
      energy[i] = nph[i] += Ran.Gaus(0.,sqrt(nph[i]));   // smear the number of generated photons
      // fraction of 4pi seen by PMT directly
      A1 = PMT_SURFACE/4./loca[i]/loca[i];
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
// toString
//------------------
const string DTOFMCResponse_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
  Get();
  if(_data.size()<=0)return string(); // don't print anything if we have no data!

  // Put the class specific code to produce nicely formatted ASCII here.
  // The JFactory_base class has several methods defined to help. They
  // rely on positions of colons (:) in the header. Here's an example:
  
  printheader( "id: orientation: end:    t [ns]:    x/y (orth.):   dE [MeV]:" );

  for(unsigned int i=0; i<_data.size(); i++ ){

    DTOFMCResponse *myTOF = _data[i];
    
    printnewrow();
    printcol("%d",	myTOF->id );
    printcol("%d",	myTOF->orientation );
    printcol("%2.3f",	myTOF->y );
    printcol("%1.3f",	myTOF->t_north );
    printcol("%1.3f",	myTOF->E_north );
    printcol("%1.3f",	myTOF->t_south );
    printcol("%1.3f",	myTOF->E_south );
    printcol("%1.3f",	myTOF->ADC_north );
    printcol("%1.3f",	myTOF->ADC_south );
    printrow();
  }
  
  return _table;

}
