#include "TTCS_crs.hh"
#include <cmath>
#include "kin_funcs.h"
#include <TMath.h>
#include "GPDs.hh"
#include <TH2D.h>
#include <TCanvas.h>

using namespace std;


//typedef double (TTCS_crs::*BH_crs_section_mf)(double *, double *);

TTCS_crs::TTCS_crs()
{
  f_BH = new TF2("f_BH", BH_crs_section, 0, 360, 0, 180, 4);
  f_INT = new TF2("f_INT", INT_crs_section, 0., 360., 0., 180., 12);
  gp = new GPDs("CFFs_DD_Feb2012.dat",  17, 17, 9, 1.49, -0.20, 0.072);
}

TTCS_crs::TTCS_crs( double a_s, double a_Q2, double a_t )// s(GeV)^2, Q2(GeV)^2, t(GeV)^2
{
  Set_SQ2t(a_s, a_Q2, a_t);
  iweight = -1;
  //  BH_crs_section_mf p = &TTCS_crs::BH_crs_section;
  //cout<<"p = "<<p<<endl<<"  Address of p = "<<&p<<endl;
  f_BH = new TF2("f_BH", BH_crs_section, 0, 360, 0, 180, 4);
  f_BH->SetParameters(is, iQ2, it);
  
  f_INT = new TF2("f_INT", INT_crs_section, 0., 360., 0., 180., 12);
  double eta = iQ2/( 2*(is - M_p*M_p) - iQ2 );
  gp = new GPDs("CFFs_DD_Feb2012.dat",  17, 17, 9, iQ2, it, eta);

  gp->Set_q2_t_eta(iQ2, it, eta);
}

double TTCS_crs::BH_crs_section( double *x, double *par)
{
  double phi = x[0]/radian;
  double theta = x[1]/radian;
  
  double s = par[0];
  double Q2 = par[1];
  double t = par[2];

  //cout<<" s Q2 t phi, theta ="<<s<<"  "<<Q2<<"  "<<t<<"  "<<phi<<"  "<<theta<<endl;

  double weight;
  
  double beta = sqrt( 1 - (4*m_e*m_e)/Q2 );
  double r = sqrt( (s - Q2 - M_p*M_p)*(s - Q2 - M_p*M_p) - 4*Q2*M_p*M_p );
  double tau = Q2/( s - M_p*M_p );
  double cos_TH_Cm = (2*s*(t - 2*M_p*M_p) + (s + M_p*M_p)*(s + M_p*M_p - Q2))/sqrt(LambdaFunc(s, M_p*M_p, 0)*LambdaFunc(s, M_p*M_p, Q2));
  double sin_TH_Cm = sqrt( 1 - cos_TH_Cm*cos_TH_Cm );
  double Delta_Perp = sin_TH_Cm*r/(2*sqrt(s));
  //iDelta_Perp = Delta_Perp;
  //double Delta_Perp = sqrt((-t)*(1 - tau) - tau*tau*M_p*M_p );
  double a = beta*r*cos(theta);
 
  double b = beta*((Q2*(s - M_p*M_p - Q2) + t*(s - M_p*M_p + Q2))/r)*cos(theta) -
    beta*( (2*(s-M_p*M_p)*sqrt(Q2)*Delta_Perp)/r )*sin(theta)*cos(phi);

  double L_BH = ((Q2 - t)*(Q2 - t) - b*b)/4.;
  double L0_BH = Q2*Q2*sin(theta)*sin(theta)/4.;
  double A_BH = (s - M_p*M_p)*(s - M_p*M_p)*Delta_Perp*Delta_Perp - t*a*(a + b) - M_p*M_p*b*b - t*(4*M_p*M_p - t)*Q2 +
    (m_e*m_e/L_BH)*(  TMath::Power( (Q2 - t)*(a + b) - (s - M_p*M_p)*b, 2) + t*(4*M_p*M_p - t)*(Q2 - t)*(Q2 - t) );
  
  double B_BH = (Q2 + t)*(Q2 + t) + b*b + 8*m_e*m_e*Q2 - 4*m_e*m_e*(t + 2*m_e*m_e)*(Q2 - t)*(Q2 - t)/L_BH;
    
  double F1p = (1./((1 - t/0.71)*(1 - t/0.71)))*(1/(1 - t/(4.*M_p*M_p) ))*( 1 - 2.79*t/( 4*M_p*M_p ) );
  double F2p = (1./((1 - t/0.71)*(1 - t/0.71)))*(1/(1 - t/(4.*M_p*M_p) ))*(ammp - 1);
  //double F1n = (1./((1 - t/0.71)*(1 - t/0.71)))*(1/(1 - t/(4.*M_p*M_p) ))*(-ammn*t/(4*M_p*M_p));
  //double F2n = (1./((1 - t/0.71)*(1 - t/0.71)))*(1/(1 - t/(4.*M_p*M_p) ))*ammn;

  if( par[3] == 1 )
    { weight = L_BH/L0_BH; }
  else
    {
      weight = 1.;
    }
  
  double crs = weight*sin(theta)*(TMath::Power(alpha_em, 3)/(4*PI*(s - M_p*M_p)*(s - M_p*M_p)) )*(beta/(-t*L_BH))*
    //double crs = weight*(TMath::Power(alpha_em, 3)/(4*PI*(s - M_p*M_p)*(s - M_p*M_p)) )*(beta/(-t*L_BH))*
    ( (A_BH/(-t))*(F1p*F1p - (t/(4*M_p*M_p))*F2p*F2p ) + (F1p + F2p)*(F1p + F2p)*B_BH/2 )*
    0.389379*1e9;

  return crs;
}

//==================The Interference term ======================
double TTCS_crs::INT_crs_section(double *x, double *par)
{
  double phi = x[0]/radian;
  double theta = x[1]/radian;
  
  double s = par[0];
  double Q2 = par[1];
  double t = par[2];
  double sc_D = par[3]; // Scale of Dterm;
  double weight;
  
  double ImH = par[5];
  double ReH = par[6];
  double ImE = par[7];
  double ReE = par[8];
  double ImHtild = par[9];
  double ReHtild = par[10];
  double Dterm = par[11];

  double sigma = 1.;
    
  double beta = sqrt( 1 - (4*m_e*m_e)/Q2 );
  double r = sqrt( (s - Q2 - M_p*M_p)*(s - Q2 - M_p*M_p) - 4*Q2*M_p*M_p );
  double tau = Q2/( s - M_p*M_p );
  double eta = tau/(2 - tau);
  double cos_TH_Cm = (2*s*(t - 2*M_p*M_p) + (s + M_p*M_p)*(s + M_p*M_p - Q2))/sqrt(LambdaFunc(s, M_p*M_p, 0)*LambdaFunc(s, M_p*M_p, Q2));
  double sin_TH_Cm = sqrt( 1 - cos_TH_Cm*cos_TH_Cm );
  double Delta_Perp = sin_TH_Cm*r/(2*sqrt(s));
  //double Delta_Perp = sqrt((-t)*(1 - tau) - tau*tau*M_p*M_p );
  double a = beta*r*cos(theta);
  double b = sigma*beta*sqrt( (Q2-t)*(Q2-t) - TMath::Power((2*(s-M_p*M_p)*sqrt(Q2)*Delta_Perp)/r , 2) )*cos(theta) - 
    beta*( (2*(s-M_p*M_p)*sqrt(Q2)*Delta_Perp)/r )*sin(theta)*cos(phi);
  double L_BH = ((Q2 - t)*(Q2 - t) - b*b)/4.;
  double L0_BH = Q2*Q2*sin(theta)*sin(theta)/4.;

  double F1p = (1./((1 - t/0.71)*(1 - t/0.71)))*(1/(1 - t/(4.*M_p*M_p) ))*( 1 - 2.79*t/( 4*M_p*M_p ) );
  double F2p = (1./((1 - t/0.71)*(1 - t/0.71)))*(1/(1 - t/(4.*M_p*M_p) ))*(ammp - 1);
  //double F1n = (1./((1 - t/0.71)*(1 - t/0.71)))*(1/(1 - t/(4.*M_p*M_p) ))*(-ammn*t/(4*M_p*M_p));
  //double F2n = (1./((1 - t/0.71)*(1 - t/0.71)))*(1/(1 - t/(4.*M_p*M_p) ))*ammn;
  
  double t_min = -4*eta*eta*M_p*M_p/(1 - eta*eta);
  double M2_int = 2*sqrt( t_min - t )/M_p*(1-eta)/(1+eta)*(F1p*(ReH + sc_D*Dterm) - eta*(F1p + F2p)*ReHtild - 
							   (t/(4*M_p*M_p))*F2p*(ReE - sc_D*Dterm));
  
  if( par[4] == 1 )
    { weight = L_BH/L0_BH; }
  else
    {
      weight = 1;
    }
  //cout<<"M2_int = "<<M2_int<<endl;
  double sigma_int = sin(theta)*weight*(-TMath::Power(alpha_em, 3))/(4*PI*s*s)*(-1/t)*(M_p/sqrt(Q2))*(1/(tau*sqrt(1 - tau)))*L0_BH/L_BH*
  //double sigma_int = weight*(-TMath::Power(alpha_em, 3))/(4*PI*s*s)*(-1/t)*(M_p/sqrt(Q2))*(1/(tau*sqrt(1 - tau)))*L0_BH/L_BH*
		       cos(phi)*(1 + cos(theta)*cos(theta))/sin(theta)*M2_int*
		       0.389379*1e9;
  return sigma_int; // This is dsigma/dQ2dts\thetad\phi pbn*GeV-4*rad-1
}



// double TTCS_crs::Get_Delta_Perp() const
// {
//   return this->iDelta_Perp;
// }

double TTCS_crs::Eval_BH(double a_phi, double a_th) const
{
  f_BH->SetParameters(is, iQ2, it, iweight);
  
  return f_BH->Eval(a_phi, a_th)/sin(a_th/radian); // 1/sin(theta) is for getting d\sigma/dcos(\theta)
}

double TTCS_crs::Eval_BH(double a_s, double a_Q2, double a_t, double a_weight, double a_phi, double a_th) const// s, Q2, t, weight, phi, theta
{
  f_BH->SetParameters(a_s, a_Q2, a_t, a_weight);
  return f_BH->Eval(a_phi, a_th)/sin(a_th/radian); // 1/sin(theta) is for getting d\sigma/dcos(\theta)
}

double TTCS_crs::Eval_INT(double a_phi, double a_th, double a_sc_D) const
{
  double eta = iQ2/( 2*(is - M_p*M_p) - iQ2 );
  gp->Set_q2_t_eta(iQ2, it, eta);
  
  double ImH = gp->GetImH();
  double ReH = gp->GetReH();
  double ImE = gp->GetImE();
  double ReE = gp->GetReE();
  double ImHtild = gp->GetImHtild();
  double ReHtild = gp->GetReHtild();
  double Dterm = gp->GetDterm();
  
  f_INT->SetParameters(is, iQ2, it, a_sc_D, iweight, ImH, ReH, ImE, ReE, ImHtild, ReHtild);
  f_INT->SetParameter(11, Dterm);
  return f_INT->Eval(a_phi, a_th)/sin(a_th/radian); // 1/sin(theta) is for getting d\sigma/dcos(\theta)
}

double TTCS_crs::Eval_INT(double a_s, double a_Q2, double a_t, double a_weight, double a_phi, double a_th, double a_sc_D) const
{
  double eta = a_Q2/( 2*(a_s - M_p*M_p) - a_Q2 );
  gp->Set_q2_t_eta(a_Q2, a_t, eta);

  double ImH = gp->GetImH();
  double ReH = gp->GetReH();
  double ImE = gp->GetImE();
  double ReE = gp->GetReE();
  double ImHtild = gp->GetImHtild();
  double ReHtild = gp->GetReHtild();
  double Dterm = gp->GetDterm();
  //  cout<<"ImH   ReH  ImE  ReE  ImHtild ReHtild Dterm  "<<ImH<<"   "<<ReH<<"   "<<ImE<<"   "<<ReE<<"   "<<ImHtild<<"   "<<ReHtild<<"   "<<Dterm<<endl;

  f_INT->SetParameters(a_s, a_Q2, a_t, a_sc_D, a_weight, ImH, ReH, ImE, ReE, ImHtild, ReHtild);
  f_INT->SetParameter(11, Dterm);
  return f_INT->Eval(a_phi, a_th)/sin(a_th/radian); // 1/sin(theta) is for getting d\sigma/dcos(\theta)
}

void TTCS_crs::Set_SQ2t(double a_s, double a_Q2, double a_t)
{
  is = a_s;
  it = a_t;
  iQ2 = a_Q2;
}

// double TTCS_crs::Get_M2int() const
// {
//   return M2_int;
// }

// double TTCS_crs::Get_M2int(double a_s, double a_Q2, double a_t, double a_sc_D)
// {
//   Set_SQ2t( a_s,  a_Q2, a_t); // phi and theta 10, are just random angles (no matter), 
//   Eval_INT(10, 10, a_sc_D);
//   return M2_int;
// }

double TTCS_crs::Integral_BH_phi_th( double a_phi_min, double a_phi_max, double a_th_min, double a_th_max)
{
  f_BH->SetParameters(is, iQ2, it, iweight);
  return f_BH->Integral(a_phi_min, a_phi_max, a_th_min, a_th_max);
}

void TTCS_crs::Set_Weight( double a_weight)
{
  iweight = a_weight;
}

void TTCS_crs::Set_sc_D( double a_sc_D)
{
  isc_D = a_sc_D;
}

void TTCS_crs::Draw_BH(const char* option)
{
  f_BH->SetParameters(is, iQ2, it, iweight);
  f_BH->Draw(option);
}

void TTCS_crs::Draw_INT(const char* option, double a_sc_D)
{
  double eta = iQ2/( 2*(is - M_p*M_p) - iQ2 );
  gp->Set_q2_t_eta(iQ2, it, eta);
  
  double ImH = gp->GetImH();
  double ReH = gp->GetReH();
  double ImE = gp->GetImE();
  double ReE = gp->GetReE();
  double ImHtild = gp->GetImHtild();
  double ReHtild = gp->GetReHtild();
  double Dterm = gp->GetDterm();

  f_INT->SetParameters(is, iQ2, it, a_sc_D, iweight, ImH, ReH, ImE, ReE, ImHtild, ReHtild);
  f_INT->SetParameter(11, Dterm);
  f_INT->SetNpx(500);
  f_INT->SetNpy(500);
  f_INT->Draw(option);
}

TH2D *TTCS_crs::Get_BH_Crs_Histogream_ThPhi(const char *name)
{
  TCanvas *ctmp = new TCanvas();
  f_BH->SetParameters(is, iQ2, it, iweight);
  f_BH->SetNpx(500);
  f_BH->SetNpy(500);
  f_BH->Draw("colz");
  TH2D *h_tmp = (TH2D*)f_BH->GetHistogram();
  h_tmp->SetName(name);
  delete ctmp;
  return h_tmp;
}

TH2D *TTCS_crs::Get_INT_Crs_Histogream_ThPhi(const char *name)
{
  TCanvas *ctmp = new TCanvas(); // Without this I was getting crazy results in utility programs :-)
  double eta = iQ2/( 2*(is - M_p*M_p) - iQ2 );
  gp->Set_q2_t_eta(iQ2, it, eta);

  double ImH = gp->GetImH();
  double ReH = gp->GetReH();
  double ImE = gp->GetImE();
  double ReE = gp->GetReE();
  double ImHtild = gp->GetImHtild();
  double ReHtild = gp->GetReHtild();
  double Dterm = gp->GetDterm();

  f_INT->SetParameters(is, iQ2, it, isc_D, iweight, ImH, ReH, ImE, ReE, ImHtild, ReHtild);
  f_INT->SetParameter(11, Dterm);
  f_INT->SetNpx(500);
  f_INT->SetNpy(500);
  f_INT->Draw("colz");
  TH2D *h_tmp = (TH2D*)f_INT->GetHistogram();
  h_tmp->SetName(name);
  delete ctmp;
  return h_tmp;
}
