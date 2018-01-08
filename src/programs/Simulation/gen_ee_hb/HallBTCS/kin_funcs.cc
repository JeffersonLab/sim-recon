//#include "kin_funcs.h"
#include <math.h>


// This is from Byukling Kayanti Formula (6.3)
double LambdaFunc( double x, double y, double z )
{
  return (x - y - z)*(x - y - z) - 4*y*z;
}

//From Byukling Kayanti Formula (5.14) Page 86
double T_min( double ma_2, double mb_2, double m1_2, double m2_2, double s) // arguments are squares of masses of particles in the reaction a+b->1+2, and s is the square of the total c.m. energy i.e. (a+b)^2
{
  return ma_2 + m1_2 - (1/(2*s))*( (s + ma_2 - mb_2)*(s + m1_2 - m2_2) - sqrt( LambdaFunc(s, ma_2, mb_2)*LambdaFunc(s, m1_2, m2_2) ) );
}

//From Byukling Kayanti Formula (5.14) page 86
double T_max( double ma_2, double mb_2, double m1_2, double m2_2, double s)
{
  return ma_2 + m1_2 - (1/(2*s))*( (s + ma_2 - mb_2)*(s + m1_2 - m2_2) + sqrt( LambdaFunc(s, ma_2, mb_2)*LambdaFunc(s, m1_2, m2_2) ) );
}

double Q2_min( double s, double Eb, double M )
{
  // M is the target mass;
  double me = 0.00051;
  double Eg = (s - M*M)/(2*M);
  double E_pr = Eb - Eg;
  double P0 = sqrt(Eb*Eb - me*me);
  double P_pr = sqrt(E_pr*E_pr - me*me);
  double Q2min = 2*(Eb*E_pr - P0*P_pr - me*me);

  return Q2min;
}
