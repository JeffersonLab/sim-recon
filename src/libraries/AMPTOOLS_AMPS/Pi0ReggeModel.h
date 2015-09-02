/*
 ============================================================================
 Name        : Pi0ReggeModel.h
 Author      : Vincent Mathieu (Adapted for C++ and GlueX sim-recon by Justin Stevens)
 Version     : v1.0 June 2015
 Copyright   : MyLab
 Publication : arXiv:1505.02321
 Description : Photoproduction of a pseudoscalar : gamma + N --> 0- + N'
               Applied to gamma + p --> pi0 + p
               To change the model, change 'CGLN_Ai'
 ============================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>

using std::complex;
using namespace std;

// Structure with all the kinematical quantities
struct Kin
{
  double m1;			// Mass particle 1
  double m2;			// Mass particle 2
  double m3;			// Mass particle 3
  double m4;			// Mass particle 4
  complex<double> s;		// Mandelstam s variable
  complex<double> t;		// Mandelstam t variable
  complex<double> phi;		// phi = stu + ...
  complex<double> ks;		// photon  momentum in s-channel frame
  complex<double> qs;		// pion    momentum in s-channel frame
  complex<double> kt;		// photon  momentum in t-channel frame
  complex<double> pt;		// nucleon momentum in t-channel frame
  complex<double> E2s;		// initial nucleon energy in s-channel frame
  complex<double> E4s;		// final   nucleon energy in s-channel frame
  complex<double> E2t;		// initial nucleon energy in t-channel frame
  complex<double> E4t;		// final   nucleon energy in t-channel frame
  complex<double> coss;	        // cos \theta_s
  complex<double> sins;	        // sin \theta_s
  complex<double> cossh;	// cos \theta_s/2
  complex<double> sinsh;	// sin \theta_s/2
  complex<double> cost;	        // cos \theta_t
  complex<double> sint;	        // sin \theta_t
  complex<double> costh;	// cos \theta_t/2
  complex<double> sinth;	// sin \theta_t/2
};

/*
 ============================================================================
 Description of the subroutines:
 Pi0PhotCS_S
 	 differential cross section from s-channel amplitudes
 Pi0PhotAmpS
 	 s-channel helicity amplitudes with inputs = four vectors and helicities
 CGLN_Ai
 	 CGLN Scalar amplitudes Ai - the model has to be specified there
 CGLNA2F
 	 transformation between CGLN Fi (ouputs) and Ai (inputs)
 kin2to2
 	 four vectors of a + b --> c + d for given Ecm and theta_s
 kinematics
 	 compute all kinematical quantities (s- and t-channel)
 	 and store then in a struct Kin

 functions:
 ---------
 Cos2T
 	 return Mandelstam t from cos(theta_s)
 cgamma
 	 return gamma(z) or log( gamma(z) )
 lambda
 	 return a*a + b*b + c*c - 2*(a*b + b*c + c*a)
 snorm
 	 return
 ============================================================================
 */

void kin2to2(double Ecm, double theta, double mass[], double pa[],double pb[],double pc[], double pd[]);
complex<double> lambda(complex<double> a, double b, double c);
double snorm(double p[]);
complex<double> cgamma(complex<double> z,int OPT);
complex<double> Pi0PhotAmpS(double pa[],double pb[],double pc[], int hel[]);
void CGLN_Ai(complex<double> s, complex<double> t, complex<double> CGLNA[]);
double Pi0PhotCS_S(double E,double theta, double &BeamSigma);
void CGLNA2F(struct Kin var, complex<double> Ai[], complex<double> Fi[]);
void kinematics(complex<double> s, complex<double> t, double mass[], struct Kin *var);
