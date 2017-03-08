/*
 ============================================================================
 Name        : Pi0ReggeModel.cc
 Author      : Vincent Mathieu (Adapted for C++ and GlueX sim-recon by Justin Stevens)
 Version     : v1.0 June 2015
 Copyright   : MyLab
 Publication : arXiv:1505.02321
 Description : Photoproduction of a pseudoscalar : gamma + N --> 0- + N'
               Applied to gamma + p --> pi0 + p
               To change the model, change 'CGLN_Ai'
 ============================================================================
*/

#include "Pi0ReggeModel.h"

// *********************************************************************************
double Pi0PhotCS_S(double E, double theta, double &BeamSigma){
/*
 * compute the differential cross section from s-channel helicities
 * in micro barns/Gev^2
 */

  double MP = 0.938272046;
  double MPI = 0.1349766;
  double mass[5] = {0.0, 0.0, MP, MPI, MP};
  double pa[4], pb[4], pc[4], pd[4];
  int hel[4] = {2,1,0,1};						// hel = {1,+,0,+} (x2)
  complex<double>  res[4] ;
  double sig;
  kin2to2(E, theta, mass, pa, pb, pc, pd);
  
  res[0] = Pi0PhotAmpS(pa, pb, pc, hel);		                // hel = {1,+,0,+} (x2)
  
  hel[1] = -1;								// hel = {1,-,0,+} (x2)
  res[1] = Pi0PhotAmpS(pa, pb, pc, hel);
  
  hel[3] = -1;								// hel = {1,-,0,-} (x2)
  res[2] = Pi0PhotAmpS(pa, pb, pc, hel);
  
  hel[1] = +1;								// hel = {1,+,0,-} (x2)
  res[3] = Pi0PhotAmpS(pa, pb, pc, hel);
  
  sig = pow(abs(res[0]),2) + pow(abs(res[1]),2) + pow(abs(res[2]),2) + pow(abs(res[3]),2);
  sig = sig/32 * 389.3;
  //sig = sig/32/M_PI/(E*E-MP*MP)/(E*E-MP*MP) * 389.3; 
 
  // Sigma beam asymmetry from equation B3b from paper
  BeamSigma = real(res[0]*conj(res[2]) - res[1]*conj(res[3]));
  BeamSigma = BeamSigma/16 * 389.3;
  BeamSigma = BeamSigma/sig;
  
  return sig;
}

// *********************************************************************************
complex<double> Pi0PhotAmpS(double pa[],double pb[],double pc[], int hel[]){
  /* Regge amplitudes for gamma + p --> pi0  + p
   * See arXiv:1505.02321 for the formulas
   * Amplitudes in the S-CHANNEL !
   * Helicities are defined in the c.o.m. of (gamma,p)
   * Inputs:
   * 		pa, pb, pc 					: four vectors in the c.o.m. frame.
   * 		hel={mua, mub, muc, mud} 	: 2 x helicities of a,b,c,d
   * Outputs:
   * 		ampl						: the amplitude ; complex number
   */
  complex<double> amp = 0.0;
  complex<double> Ai[5] = {0}, Fi[5] = {0};
  double pab[4], pbc[4], pca[4], pd[4];
  struct Kin var;
  double mass[5]={0};
  
  // Only valid for Q^2 = 0 --> photon helicity should be +1 or -1
  if ( abs(hel[0]) != 2 ) return 0.0;
  
  int i;
  for(i=0;i<4;i++){							// compute the Mandelstam invariant
    pab[i] = pa[i] + pb[i];						// from four-vectors
    pbc[i] = pb[i] - pc[i];						//
    pca[i] = pc[i] - pa[i];						//
    pd[i]  = pa[i] + pb[i] - pc[i];				//
  }
  
  double EPS = 0.0;
  complex<double> I(0,1);
  var.s = snorm( pab ) + I*EPS;
  var.t = snorm( pca ) - I*EPS;
  // photon should be real
  mass[1]   = 0.0;
  mass[2]  = sqrt( snorm( pb ) );
  mass[3]  = sqrt( snorm( pc ) );
  mass[4]  = sqrt( snorm( pd ) );
  
  kinematics(var.s, var.t, mass, &var);	// fill the structure var with all kin. quantities
  CGLN_Ai(var.s,var.t,Ai);				// call the model : CHANGE HERE FOR ANOTHER MODEL
  
  CGLNA2F(var, Ai, Fi);					// Convert CGLN Ai to CGLN Fi

  // Test nucleon helicities hel[1] and hel[3]
  if ( hel[1] == 1 && hel[3] == 1 ) {
    // hel = {1,+,0,+}
    amp = sqrt(2) * var.sinsh * (Fi[2] + Fi[1] )
      + 1/sqrt(2) * var.sins * var.cossh * (Fi[3] + Fi[4]);
  }
  else if  ( hel[1] == -1 && hel[3] == -1 ){
    // hel = {1,-,0,-}
    amp = -1/sqrt(2) * var.sins * var.cossh * (Fi[3] + Fi[4]);
  }
  else if ( hel[1] == -1 && hel[3] == 1 ){
    // hel = {1,-,0,+}
    amp = sqrt(2) * var.cossh * (Fi[2] - Fi[1] )
      + 1/sqrt(2) * var.sins * var.sinsh * (Fi[3] - Fi[4]);
  }
  else if ( hel[1] == 1 && hel[3] == -1 ){
    // hel = {1,+,0,-}
    amp = 1/sqrt(2) * var.sins * var.sinsh * (Fi[3] - Fi[4]);
  }
  else return 0.0 ;
  
  // if negative photon helicity, nucleon flip amplitude changes sign
  if (hel[0] == -2 && hel[1] != hel[3] ) amp = -1.*amp;
  
  // there is a factor 8*PI*W between Hi and helicity amplitudes
  return amp * 8. * M_PI * sqrt(var.s);
}

// *********************************************************************************
void CGLN_Ai(complex<double> s, complex<double> t, complex<double> CGLNA[]){
  /*
   * Model for gamma p --> pi0 p from arXiv:1505.02321
   * Vincent Mathieu June 2015
   */

  // Model parameters
  double alpV[3] = { 0.442457, 1.099910, 0.0};		// vector trajectory
  double alpC[3] = { 0.461811, 0.166543, 0.0};		// vector cut trajectory
  double alpA[3] = {-0.193332, 1.021300, 0.0};		// axial-vector trajectory
  double g1 = 3.8873, g4 = -10.1403, g1c = -1.76187, g4c = -3.58089, g2 = -8.887; // residues
  
  complex<double> Rv, Rc, Ra;
  complex<double> avec, acut, aaxi;
  // trajectories:
  avec = alpV[0] + t*alpV[1] + t*t*alpV[2];
  acut = alpC[0] + t*alpC[1] + t*t*alpC[2];
  aaxi = alpA[0] + t*alpA[1] + t*t*alpA[2];
  
  // Regge factors:
  complex<double> I(0,1);
  Rv = cgamma( 1.0 - avec, 0)/2. * ( 1.-exp(-1.*I*M_PI*avec) ) * pow(s,avec-1.);
  Rc = cgamma( 1.0 - acut, 0)/2. * ( 1.-exp(-1.*I*M_PI*acut) ) * pow(s,acut-1.);
  Ra = cgamma( 1.0 - aaxi, 0)/2. * ( 1.-exp(-1.*I*M_PI*aaxi) ) * pow(s,aaxi-1.);
  Rc = Rc / log(s);
  // IF ONLY VECTOR POLE
  //Rc = 0; Ra = 0;
  
  // Scalar amplitudes:
  CGLNA[1] = t* ( g1*Rv + g1c*Rc);
  CGLNA[2] = g2 * Ra - ( g1*Rv + g1c*Rc);
  CGLNA[3] = 0;
  CGLNA[4] = g4*Rv + g4c*Rc;
  
  return ;

}

// *********************************************************************************
void CGLNA2F(struct Kin var, complex<double> Ai[], complex<double> Fi[]){
  /*
   * Compute CGLN Fi(s,t) from CGLN Ai(s,t)
   */
  complex<double> w;
  complex<double> E2p, E2m, E4p, E4m;
  
  w = sqrt(var.s);
  E2p = var.E2s + var.m2 ;
  E2m = var.E2s - var.m2 ;
  E4p = var.E4s + var.m4 ;
  E4m = var.E4s - var.m4 ;
  
  Fi[1] = (w - var.m2) * Ai[1] + (var.m3*var.m3 - var.t)/2. * (Ai[3] - Ai[4]);
  Fi[1] = Fi[1]  + ( w - var.m2 ) * ( w - var.m4 ) * Ai[4];
  Fi[1] = Fi[1] * sqrt( E2p * E4p ) / ( 8. * M_PI * w);

  Fi[2] = -1.*(w + var.m2) * Ai[1] + (var.m3*var.m3 - var.t)/2. * (Ai[3] - Ai[4]);
  Fi[2] = Fi[2]  + ( w + var.m2 ) * ( w + var.m4 ) * Ai[4];
  Fi[2] = Fi[2] * sqrt( E2m * E4m ) / ( 8. * M_PI * w);
  
  Fi[3] = (w + var.m2) * (  (w - var.m2)*Ai[2] + Ai[3] - Ai[4] );
  Fi[3] = Fi[3] * sqrt( E2m * E4p) * var.qs  / ( 8. * M_PI * w);
  
  Fi[4] = (w - var.m2) * ( -1.*(w + var.m2)*Ai[2] + Ai[3] - Ai[4] );
  Fi[4] = Fi[4] * sqrt( E2p * E4m) * var.qs  / ( 8. * M_PI * w);
  
  return;
}

// *********************************************************************************
void kin2to2(double Ecm, double theta, double mass[], double pa[],double pb[],double pc[], double pd[]){
/*
 * Kinematics of the 2-to-2 reaction, a + b --> c + d
 * Inputs:
 * 		Ecm center of mass energy
 * 		cos cosine of the scattering angle in the center of mass
 * 		mass = {ma,mb,mc,md} vector with the masses of external particles
 * 	Outputs:
 * 		pa, pb, pc, pd are the momenta of the particles
 */
  double Ea, Eb, Ec, Ed; 			// Energies of the particles
  double pi, pf;					// initial and final breakup momenta
  double ma, mb, mc, md;			// masses of the particles

  ma = mass[1]; mb = mass[2]; mc = mass[3]; md = mass[4];
  
  // Check that inputs are valid
  if( ( ma<0 ) || ( mb<0 ) || ( mc<0 ) || ( md<0 )  ) {
    printf("\n*** Wrong masses in kin2to2 ! *** \n\n");
    return ;
  }
  if( ( Ecm < ma + mb ) || ( Ecm < mc + md )  ) {
    printf("\n*** Wrong total energy in kin2to2 ! *** \n\n");
    return ;
  }
  
  Ea = ( Ecm*Ecm + ma*ma - mb*mb )/2./Ecm;
  Eb = ( Ecm*Ecm - ma*ma + mb*mb )/2./Ecm;
  Ec = ( Ecm*Ecm + mc*mc - md*md )/2./Ecm;
  Ed = ( Ecm*Ecm - mc*mc + md*md )/2./Ecm;
  
  pi = sqrt(Ea*Ea - ma*ma);
  pa[0] = Ea; pa[1] = 0; pa[2] = 0; pa[3] = +pi;
  pb[0] = Eb; pb[1] = 0; pb[2] = 0; pb[3] = -pi;

  pf = sqrt(Ec*Ec - mc*mc);
  pc[0] = Ec; pc[1] = +pf*sin(theta); pc[2] = 0; pc[3] = +pf*cos(theta);
  pd[0] = Ed; pd[1] = -pf*sin(theta); pd[2] = 0; pd[3] = -pf*cos(theta);
  
  return;
}

// *********************************************************************************
void kinematics(complex<double> s, complex<double> t, double mass[], struct Kin *var)
{
  double m12, m22, m32, m42; 	// masses squared
  complex<double> t0, t1, u ;
  
  var->s = s;
  var->t = t;
  
  var->m1 = mass[1];
  var->m2 = mass[2];
  var->m3 = mass[3];
  var->m4 = mass[4];
  
  m12 = mass[1] * mass[1] ;
  m22 = mass[2] * mass[2] ;
  m32 = mass[3] * mass[3] ;
  m42 = mass[4] * mass[4] ;
  
  var->ks = sqrt( lambda(s, m12, m22) / 4. / s );
  var->qs = sqrt( lambda(s, m32, m42) / 4. / s );
  var->kt = sqrt( lambda(t, m12, m32) / 4. / t );
  var->pt = sqrt( lambda(t, m22, m42) / 4. / t );
  
  // nuclei energies in s- and t-channel frames
  var->E2s = ( s + m22 - m12 ) / 2. / sqrt(s);
  var->E4s = ( s + m42 - m32 ) / 2. / sqrt(s);
  var->E2t = ( t + m22 - m42 ) / 2. / sqrt(t);
  var->E4t = ( t + m42 - m22 ) / 2. / sqrt(t);
  
  t1 = pow(m12 - m32 - m22 + m42,2)/(4.*s) - (var->ks - var->qs) * (var->ks - var->qs);
  t0 = t1 - 4. * var->ks *var->qs ;
  u = -1.*s - t + m12 + m22 + m32 + m42 ;	// Mandelstam s variable
  var->phi = s * (t - t1) * (t0 - t) ;	// Kibble function
  
  var->coss = 1. + (t - t1)/(2. * var->qs * var->ks) ;
  var->cost = (t*(s-u) + (m12-m32) * (m22-m42) )/(sqrt(lambda(var->t,m12,m32) * lambda(var->t,m22,m42) ) ) ;
  var->sins = sqrt( var->phi / s ) / ( 2. * var->qs * var->ks) ;
  var->sint = sqrt( var->phi / t ) / ( 2. * var->kt * var->pt) ;
  var->cossh = sqrt( (1. + var->coss ) / 2. );
  var->sinsh = sqrt( (1. - var->coss ) / 2. );
  var->costh = sqrt( (1. + var->cost ) / 2. );
  var->sinth = sqrt( (1. - var->cost ) / 2. );

  return ;
}

// *********************************************************************************
complex<double> lambda(complex<double> a, double b, double c){
// triangle function
  return a*a + b*b + c*c - 2.*(a*b + b*c + c*a);
}

// *********************************************************************************
double snorm(double p[]){
// Norm squared of a quadri-vector ; p[0] is the energy ; p[1-3] are x,y,z components
  return p[0]*p[0] - ( p[1]*p[1] + p[2]*p[2] + p[3]*p[3] );
}

// *********************************************************************************
complex<double> cgamma(complex<double> z,int OPT)
{
  complex<double> I(0,1);
  complex<double> g, infini= 1e308+ 0.*I; // z0,z1
  double x0,q1,q2,x,y,th,th1,th2,g0,gr,gi,gr1,gi1;
  double na,t,x1,y1,sr,si;
  int j,k;
  x1=9e9;
  na=9e9;
  
  static double a[] = {
    8.333333333333333e-02,
    -2.777777777777778e-03,
    7.936507936507937e-04,
    -5.952380952380952e-04,
    8.417508417508418e-04,
    -1.917526917526918e-03,
    6.410256410256410e-03,
    -2.955065359477124e-02,
    1.796443723688307e-01,
    -1.39243221690590};
  
  x = real(z);
  y = imag(z);
  if (x > 171) return infini;
  if ((y == 0.0) && (x == (int)x) && (x <= 0.0))
    return infini;
  else if (x < 0.0) {
    x1 = x;
    y1 = y;
    x = -x;
    y = -y;
  }
  x0 = x;
  if (x <= 7.0) {
    na = (int)(7.0-x);
    x0 = x+na;
  }
  q1 = sqrt(x0*x0+y*y);
  th = atan(y/x0);
  gr = (x0-0.5)*log(q1)-th*y-x0+0.5*log(2.0*M_PI);
  gi = th*(x0-0.5)+y*log(q1)-y;
  for (k=0;k<10;k++){
    t = pow(q1,-1.0-2.0*k);
    gr += (a[k]*t*cos((2.0*k+1.0)*th));
    gi -= (a[k]*t*sin((2.0*k+1.0)*th));
  }
  if (x <= 7.0) {
    gr1 = 0.0;
    gi1 = 0.0;
    for (j=0;j<na;j++) {
      gr1 += (0.5*log((x+j)*(x+j)+y*y));
      gi1 += atan(y/(x+j));
    }
    gr -= gr1;
    gi -= gi1;
  }
  if (x1 <= 0.0) {
    q1 = sqrt(x*x+y*y);
    th1 = atan(y/x);
    sr = -sin(M_PI*x)*cosh(M_PI*y);
    si = -cos(M_PI*x)*sinh(M_PI*y);
    q2 = sqrt(sr*sr+si*si);
    th2 = atan(si/sr);
    if (sr < 0.0) th2 += M_PI;
    gr = log(M_PI/(q1*q2))-gr;
    gi = -th1-th2-gi;
    x = x1;
    y = y1;
  }
  if (OPT == 0) {
    g0 = exp(gr);
    gr = g0*cos(gi);
    gi = g0*sin(gi);
  }
  g = gr + I*gi;
  return g;
}
