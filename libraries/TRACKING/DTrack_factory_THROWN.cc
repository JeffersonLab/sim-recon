// $Id$
//
//    File: DTrack_factory_THROWN.cc
// Created: Mon Sep  3 19:57:11 EDT 2007
// Creator: davidl (on Darwin Amelia.local 8.10.1 i386)
//

#include <cmath>
using namespace std;

#include "DANA/DApplication.h"

#include "DTrack_factory_THROWN.h"
#include "DMCThrown.h"
#include "DReferenceTrajectory.h"
#include "DRandom.h"
#include "DMatrix.h"


//------------------
// evnt
//------------------
jerror_t DTrack_factory_THROWN::evnt(JEventLoop *loop, int eventnumber)
{
  vector<const DMCThrown*> mcthrowns;
  loop->Get(mcthrowns);

  for(unsigned int i=0; i< mcthrowns.size(); i++){
    const DMCThrown *thrown = mcthrowns[i];

    if(fabs(thrown->q)<1)continue;

    // Create new DTrack object and initialize parameters with those
    // from track candidate
    /*
     * Dave Lawrence's preliminary smearing algorithm.
     DTrack *track = new DTrack;
     track->q			= thrown->q;
     track->p			= thrown->p*(1.0+SampleGaussian(0.02));
     track->theta	= thrown->theta + SampleGaussian(0.001);
     if(track->theta<0.0)track->theta = 0.0;
     if(track->theta>=M_PI)track->theta = M_PI;
     track->phi		= thrown->phi + SampleGaussian(0.002);
     if(track->phi<0.0)track->phi+=2.0*M_PI;
     if(track->phi>=2.0*M_PI)track->phi-=2.0*M_PI;
     track->x			= thrown->x + SampleGaussian(0.01);
     track->y			= thrown->y + SampleGaussian(0.01);
     track->z			= thrown->z + SampleGaussian(1.00);
     track->candidateid = 0;
     track->chisq	= 1.0;
     */
    DTrack *track = new DTrack;
    track->q			= thrown->q;
    track->candidateid = 0;
    track->chisq	= 1.0;
    track->x			= thrown->x;
    track->y			= thrown->y;
    track->z			= thrown->z;
    track->p			= thrown->p;
    track->theta	= thrown->theta;
    track->phi		= thrown->phi;

    // Adapted from Alex's fortran code to paramatrize the resolution 
    // of GlueX. He developed this from a study Dave Lawrence did. 

    // Both the reference trajectory and the kinematic data section below
    // use DVector3 objects for position and momentum.
    DVector3 mom, pos;
    pos.SetXYZ(track->x, track->y, track->z);
    mom.SetMagThetaPhi(track->p, track->theta, track->phi);

    // We need to swim a reference trajectory here. To avoid the overhead
    // of allocating/deallocating them every event, we keep a pool and
    // re-use them. If the pool is not big enough, then add one to the
    // pool.
    if(rt.size()<=_data.size()){
      // This is a little ugly, but only gets called a few times throughout the life of the process
      // Note: these never get deleted, even at the end of process.
      DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
      rt.push_back(new DReferenceTrajectory(dapp->GetBfield()));
    }
    rt[_data.size()]->Swim(pos, mom, track->q);
    track->rt = rt[_data.size()];

    // Create and fill the covariance matrix for the track.
    // We need to fill this using errors estimated from the thrown
    // momentum and angle. 
    DMatrixDSym errMatrix(1,7);

    // Fill in DKinematicData part
    track->setMass(0.0);
    track->setMomentum(mom);
    track->setPosition(pos);
    //track->setCharge(track->q);
    track->setErrorMatrix(errMatrix);

    this->SmearMomentum(track);

    _data.push_back(track);
  }

  return NOERROR;
}

//------------------
// SampleGaussian
//------------------
double DTrack_factory_THROWN::SampleGaussian(double sigma)
{
  // We loop to ensure not to return values greater than 3sigma away
  double val;
  do{
    double epsilon = 1.0E-10;
    double r1 = epsilon+((double)random()/(double)RAND_MAX);
    double r2 = (double)random()/(double)RAND_MAX;
    val = sigma*sqrt(-2.0*log(r1))*cos(2.0*M_PI*r2);
  }while(fabs(val/sigma) > 3.0);

  return val;
}

//------------------
// toString
//------------------
const string DTrack_factory_THROWN::toString(void)
{
  // Ensure our Get method has been called so _data is up to date
  GetNrows();
  if(_data.size()<=0)return string(); // don't print anything if we have no data!

  printheader("row: q:       p:   theta:   phi:    x:    y:    z:");

  for(unsigned int i=0; i<_data.size(); i++){

    DTrack *track = _data[i];

    printnewrow();

    printcol("%x", i);
    printcol("%+d", (int)track->q);
    printcol("%3.3f", track->p);
    printcol("%1.3f", track->theta);
    printcol("%1.3f", track->phi);
    printcol("%2.2f", track->x);
    printcol("%2.2f", track->y);
    printcol("%2.2f", track->z);

    printrow();
  }

  return _table;
}

//------------------
// The various smearing resolutions taken
// directly from Alex's smear.f code.
//------------------
void DTrack_factory_THROWN::SmearMomentum(DTrack *trk)
{
  DRandom rnd;

  // Track phi, theta, mag
  //double pmag = trk->p;
  //double phi = trk->phi;
  //double theta = trk->theta;
  double pmag = trk->lorentzMomentum().Rho();
  double theta = trk->lorentzMomentum().Theta();
  double phi = trk->lorentzMomentum().Phi();

  // Alex's study used phi in degrees
  double theta_deg = theta*180.0/TMath::Pi();

  // Set the random seed based on phi
  rnd.SetSeed((uint)(10000*phi));

  double dp, dphi, dtheta, dtheta_deg;

  // Grabbed this from 
  if(theta_deg<15) 
    {dp=res_p1(theta_deg); dtheta_deg=res_p10(pmag); dphi=res_p9(pmag);}
  else if(theta_deg>=15 && theta_deg<25) 
    {dp=res_p2(pmag); dtheta_deg=res_p10(pmag); dphi=res_p9(pmag);}
  else if(theta_deg>=25 && theta_deg<35) 
    {dp=res_p3(pmag); dtheta_deg=res_p10(pmag); dphi=res_p9(pmag);}
  else if(theta_deg>=35 && theta_deg<45) 
    {dp=res_p4(pmag); dtheta_deg=res_p11(pmag); dphi=res_p9(pmag);}
  else if(theta_deg>=45 && theta_deg<55) 
    {dp=res_p4(pmag); dtheta_deg=res_p12(pmag); dphi=res_p9(pmag);}
  else if(theta_deg>=55)
  {
    dtheta_deg=res_p12(pmag); 
    dphi=res_p9(pmag);
    if(pmag<2.0) 
      dp=res_p8(theta_deg);  /// Note that this is different from the FORTRAN code due to mismatch from orginal study.
    else if(pmag>=2.0 && pmag<3.0) 
      dp=res_p7(theta_deg); 
    else if(pmag>=3.0) /// Note that this is different from the FORTRAN code due to mismatch from orginal study.
      dp=res_p6(theta_deg); 
  }

  // Both dphi and dtheta are returned in degrees, so we convert back to radians
  dphi *= TMath::Pi()/180.0;
  dtheta = dtheta_deg*TMath::Pi()/180.0;

  // Generate the spearing
  double deltap = rnd.Gaus(0.0, dp);
  double deltatheta = rnd.Gaus(0.0, dtheta);
  double deltaphi = rnd.Gaus(0.0, dphi);
  pmag += deltap;
  theta += deltatheta;
  phi += deltaphi;

  // Set the 3momentum of the track with the new parameters
  DVector3 mom;
  mom.SetMagThetaPhi(pmag, theta, phi);
  trk->setMomentum(mom);
  trk->setMass(0.0);

  // Set up the error matrix properly
  DMatrixDSym errMat(7);
  DMatrix sphErrMat(3,3);
  DMatrix transformMat(3,3);
  DMatrix transformMat_T(3,3);
  DMatrix cartErrMat(3,3);

  // Start off with a diagonal matrix in the spherical coordinates.
  sphErrMat[0][0] = dp*dp;
  sphErrMat[1][1] = dtheta*dtheta;
  sphErrMat[2][2] = dphi*dphi;

  // Transformation matrix
  // T_ij = dr_i/dx_j
  transformMat[0][0] = sin(theta)*cos(phi);
  transformMat[1][0] = pmag*cos(theta)*cos(phi);
  transformMat[2][0] = -pmag*sin(theta)*sin(phi);
  transformMat[0][1] = sin(theta)*sin(phi);
  transformMat[1][1] = pmag*cos(theta)*sin(phi);
  transformMat[2][1] = pmag*sin(theta)*cos(phi);
  transformMat[0][2] = cos(theta);
  transformMat[1][2] = -pmag*sin(theta);
  transformMat[2][2] = 0.0;

  // Carry out the conversion to cartesian coordinates
  transformMat_T.Transpose(transformMat);
  cartErrMat = transformMat_T * sphErrMat * transformMat;

  // Set up the full error matrix.
  // Note that we have not smeared the energy [3] nor
  // the vertex [4,5,6].
  errMat[0][0] = cartErrMat[0][0];
  errMat[0][1] = cartErrMat[0][1];
  errMat[0][2] = cartErrMat[0][2];
  errMat[1][0] = cartErrMat[1][0];
  errMat[1][1] = cartErrMat[1][1];
  errMat[1][2] = cartErrMat[1][2];
  errMat[2][0] = cartErrMat[2][0];
  errMat[2][1] = cartErrMat[2][1];
  errMat[2][2] = cartErrMat[2][2];

  // Set the error matrix
  trk->setErrorMatrix(errMat);
}

double DTrack_factory_THROWN::res_p1(double x)
{
  // Return dpmag/pmag(theta) 
  // 0<theta<15 
  // 0<pmag<inf
  double k0=0.11088;
  double k1=-0.023485;
  double k2=0.0021123;
  double k3=-0.000058412;

  double res = k0 + k1*x + k2*pow(x,2)+ k3*pow(x,3);
  return res;
}

double DTrack_factory_THROWN::res_p2(double x)
{
  // Return dpmag/pmag(theta) 
  // 15<theta<25 
  // 0<pmag<inf
  double a=0.013814;
  double b=0.010706;
  double s=0.342;

  double res=sqrt(pow(a*s*x,2)+pow(b,2)*(1+0.02/pow(x,2))/s);
  return res;
}

double DTrack_factory_THROWN::res_p3(double x)
{
  // Return dpmag/pmag(theta) 
  // 25<theta<35 
  // 0<pmag<inf
  double a=0.0073285;
  double b=0.0083813;
  double s=0.5;

  double res=sqrt(pow(a*s*x,2)+pow(b,2)*(1+0.02/pow(x,2))/s);
  return res;
}

double DTrack_factory_THROWN::res_p4(double x)
{
  // Return dpmag/pmag(theta) 
  // 35<theta<45 
  // 0<pmag<inf
  double a=0.007551;
  double b=0.0087766;
  double s=0.643;

  double res=sqrt(pow(a*s*x,2)+pow(b,2)*(1+0.02/pow(x,2))/s);
  return res;
}

double DTrack_factory_THROWN::res_p5(double x)
{
  // Return dpmag/pmag(theta) 
  // 45<theta<55 
  // 0<pmag<inf
  double a=0.0077688;
  double b=0.009204;
  double s=0.766;

  double res=sqrt(pow(a*s*x,2)+pow(b,2)*(1+0.02/pow(x,2))/s);
  return res;
}

double DTrack_factory_THROWN::res_p6(double x)
{
  // Return dpmag/pmag(theta) 
  // 55<theta<180 
  // 0<pmag<2 GeV
  double k0=0.37978;
  double k1=-0.0178;
  double k2=0.0003155;
  double k3=-.0000024105;
  double k4=0.000000006728;

  double res = k0 + k1*x + k2*pow(x,2) + k3*pow(x,3)+ k4*pow(x,4);;
  return res;
}

double DTrack_factory_THROWN::res_p7(double x)
{
  // Return dpmag/pmag(theta) 
  // 55<theta<180 
  // 2<pmag<3 GeV
  double k0=0.36913;
  double k1=-0.017553;
  double k2=0.00031951;
  double k3=-.0000024915;
  double k4=0.000000007067;

  double res = k0 + k1*x + k2*pow(x,2) + k3*pow(x,3)+ k4*pow(x,4);;
  return res;
}

double DTrack_factory_THROWN::res_p8(double x)
{
  // Return dpmag/pmag(theta) 
  // 55<theta<180 
  // 3<pmag<inf
  double k0=0.5213;
  double k1=-0.024712;
  double k2=0.00044629;
  double k3=-.0000034504;
  double k4=0.0000000096985;

  double res = k0 + k1*x + k2*pow(x,2) + k3*pow(x,3)+ k4*pow(x,4);;
  return res;
}

double DTrack_factory_THROWN::res_p9(double x)
{
  // Return dphi(pmag)
  // 0<theta<180 
  // 0<pmag<inf
  double k0=0.95004;
  double k1=-0.71694;
  double k2=0.26776;
  double k3=-0.044906;
  double k4=0.0027684;

  double res = k0 + k1*x + k2*pow(x,2) + k3*pow(x,3)+ k4*pow(x,4);;
  if(x>=5.5768) 
  {
    res = 0.16847;
  }
  return res;
}

double DTrack_factory_THROWN::res_p10(double x)
{
  // Return dtheta/theta(p) 
  // 0<theta<35 
  // 0<pmag<inf
  double xo=0.28323;
  double yo=0.092655;
  double a1=0.13457;
  double a2=0.16388;
  double t1=0.27009;
  double t2=0.95528;

  double res=yo+ a1*exp(-(x-xo)/t1) + a2*exp(-(x-xo)/t2);
  if(x>=4)
  {
    res=0.096163;
  }
  return res;
}

double DTrack_factory_THROWN::res_p11(double x)
{
  // Return dtheta/theta(p) 
  // 35<theta<45 
  // 0<pmag<inf
  double xo=0.243846;
  double yo=0.14657;
  double a1=0.27375;
  double a2=0.13852;
  double t1=0.15888;
  double t2=0.76164;

  double res=yo+ a1*exp(-(x-xo)/t1) + a2*exp(-(x-xo)/t2);
  if(x>=4)
  {
    res=0.1478;
  }
  return res;
}

double DTrack_factory_THROWN::res_p12(double x)
{
  // Return dtheta/theta(p) 
  // 45<theta<180 
  // 0<pmag<inf
  double xo=0.235728;
  double yo=0.20553;
  double a1=0.16369;
  double a2=0.062797;
  double t1=0.28976;
  double t2=1.056;

  double res=yo+ a1*exp(-(x-xo)/t1) + a2*exp(-(x-xo)/t2);
  if(x>=4)
  {
    res=0.20749;
  }
  return res;
}

