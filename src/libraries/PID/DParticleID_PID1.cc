// $Id$
//
//    File: DParticleID_PID1.cc
// Created: Mon Feb 28 15:25:35 EST 2011
// Creator: staylor (on Linux ifarml1 2.6.18-128.el5 x86_64)
//

#include "DParticleID_PID1.h"

//---------------------------------
// DParticleID_PID1    (Constructor)
//---------------------------------
DParticleID_PID1::DParticleID_PID1(JEventLoop *loop):DParticleID(loop)
{

}

//---------------------------------
// ~DParticleID_PID1    (Destructor)
//---------------------------------
DParticleID_PID1::~DParticleID_PID1()
{

}


// Compute a metric for how well the measured dEdx agrees with the expected
// dEdx for a given particle type. Also returns the dEdx itself and the number
// of hits that went into the calculation.
jerror_t DParticleID_PID1::GetdEdxChiSq(const DTrackTimeBased *track,
					double &dEdx,
					unsigned int &num,double &chi2){
  //initialization
  dEdx=chi2=0.;
  num=0;

  // Get the dEdx info from the CDC/FDC hits on the track
  vector<DParticleID::dedx_t>dEdx_list;
  GetdEdx(track,dEdx_list);
  
  // if the track does not intersect with any of the hit wires, then the 
  // parameters are clearly wrong for the set of hits!
  if (dEdx_list.size()==0) return VALUE_OUT_OF_RANGE;
  
  // Truncated mean:  loop over a subset of this list, throwing away a
  // number of the highest dE/dx values.  Since the list is sorted according 
  // to dEdx values, with smaller values being earlier in the list, we need 
  // only loop over a fraction of the total number of hits.
  double dEdx_diff=0.;
  double p_avg=0.;
  double mean_path_length=0.;
  num=dEdx_list.size()/2+1;
  double N=double(num);
  for (unsigned int i=0;i<num;i++){
    double p=dEdx_list[i].p;
    double dx=dEdx_list[i].dx;
    double dE=dEdx_list[i].dE;						
    double my_dedx=dE/dx;
    
    // Get the expected (most probable) dE/dx for a particle with this mass
    // and momentum for this hit
    double dEdx_mp=GetMostProbabledEdx(p,track->mass(),dx);

    dEdx+=my_dedx;
    dEdx_diff+=my_dedx-dEdx_mp;
    p_avg+=p;
    mean_path_length+=dx;
  }
  dEdx/=N; 
  dEdx_diff/=N;
  mean_path_length/=N;
  p_avg/=N;
 
  // Estimated uncertainty in the dEdx measurement
  double dEdx_sigma=GetdEdxSigma(N,p_avg,track->mass(),mean_path_length);
  
  // Chi2 for dedx measurement
  chi2=dEdx_diff*dEdx_diff/dEdx_sigma/dEdx_sigma;

  return NOERROR;
}
