#ifndef TTCS_crs_HH
#define TTCS_crs_HH

#include <TF1.h>
#include <TF2.h>
#include "GPDs.hh"

class TH2D;

class TTCS_crs: public TObject
{
public:
  TTCS_crs();
  TTCS_crs(double, double, double); // s(GeV)^2, Q2(GeV)^2, t(GeV)^2
  double Eval_BH(double, double) const;// phi and theta in radians
  double Eval_BH(double, double, double, double, double, double) const; // s, Q2, t, weight, phi, theta
  double Eval_INT(double, double, double) const;// (phi, theta, scale of Dterm) phi and theta in radians
  double Eval_INT(double, double, double, double, double, double, double) const; // s, Q2, t, weight, phi, theta, sc_D
  //double Get_Delta_Perp() const;
  // double Get_M2int() const; // M^{--}
  // double Get_M2int(double, double, double, double ); // M^{--}
  void Set_SQ2t(double, double, double);
  double Integral_BH_phi_th( double phi_min = 0, double phi_max = 360, double th_min = 0, double th_max = 180);
  void Set_Weight( double weight = -1.); // with +1 it will weight with L/L0, otherwise it will not
  void Set_sc_D(double); // Set magnitude of D-term
  void Draw_BH(const char* option);
  void Draw_INT(const char* option, double sc_D = 1.);
  TH2D *Get_BH_Crs_Histogream_ThPhi(const char *name);
  TH2D *Get_INT_Crs_Histogream_ThPhi(const char *name);

private:
  //  double M2_int;
  double is, iQ2, it;
  double iweight;
  double isc_D;
  GPDs *gp;
  
  //static double iDelta_Perp; // This 'static' is not a good idea, but it didn't work without this...

  static constexpr double radian = 57.2957795130823229;
  static constexpr double m_e = 0.00051;
  static constexpr double M_p = 0.938272;
  static constexpr double alpha_em = 1./137.;
  static constexpr double PI = 3.14159265358979312;
  static constexpr double ammp = 2.793;
  static constexpr double ammn = -1.913;

  static double BH_crs_section( double *, double *); // Hmmm why it worked with static ?, and didn't work without static?
  static double INT_crs_section( double *, double *); // Hmmm why it worked with static ?, and didn't work without static?
  TF2 *f_BH;
  TF2 *f_INT;
};

#endif
