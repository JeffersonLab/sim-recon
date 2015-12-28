//
// CobremsGenerator - class implementation
//
// author: richard.t.jones at uconn.edu
// version: july 27, 2015
//
// notes:
//
// This class computes the spectrum of bremsstrahlung radiation from a
// crystal radiator.  The formalism is that described in the following paper.
//
//  W. Kaune, G. Miller, W. Oliver, R.W. Williams, and K.K. Young,
//  "Inclusive cross sections for pion and proton production by photons
//  using collimated coherent bremsstrahlung", Phys Rev D, vol 11,
//  no 3 (1975) pp. 478-494.
//
// The model for the photon beam contains the following parameters.
// 1. electron beam
//    * beam energy: mean and rms spread
//    * spot on radiator: gaussian model, cylindrical symmetry
//    * emittance: gaussian model, cylindrical symmetry
// 2. crystal target
//    * implemented for diamond, silicon
//    * uniform thickness across beam spot
//    * mosaic spread: gaussian model
//    * dipole atomic form factor
//    * Debye-Waller factor: defines coherent domain in q sum
// 3. downstream collimator
//    * fixed distance from radiator
//    * sharp cutoff at collimator radius
//    * perfect alignment with beam axis assumed
//
// The crystal orientation is computed based on the requested coherent
// edge position requested by the user. For a high-energy electron beam
// this fixes one of the angles of the crystal with respect to the
// electron beam. The other angle must be chosen based on other
// considerations. A default value for this secondary angle parameter
// is assigned below, based on the observation that
//  a) it is significantly larger than the primary edge-defining
//     angle, so that additional peaks from reciprocal lattice sites
//     (hkl) with the same h and different k values do not contribute
//     significantly to the spectrum below the endpoint, and
//  b) it is small enough, to render it unlikely that a random lattice
//     vector from a distant region in q-space will cross through the
//     coherent enhancement region as the primary (220) peak is moved
//     through its full range in x from 30% to 90% of the endpoint
// for beamline parameters similar to those describing Hall D and GlueX.
// Should the user wish to try other values for this angle, a public
// method is provided for this purpose.

#define COBREMS_GENERATOR_VERBOSITY 1

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <CobremsGenerator.hh>
using namespace std;

#if BOOST_PYTHON_WRAPPING
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/special_functions/erf.hpp>
#endif  // BOOST_PYTHON_WRAPPING

#ifndef _DBG_
#define _DBG_ cout<<__FILE__<<":"<<__LINE__<<" "
#define _DBG__ cout<<__FILE__<<":"<<__LINE__<<endl
#endif

const double CobremsGenerator::dpi = 3.1415926535897;
const double CobremsGenerator::me = 0.510998910e-3;
const double CobremsGenerator::alpha = 7.2973525698e-3;
const double CobremsGenerator::hbarc = 0.1973269718e-15;

// Replaced boost::math::expint(1,z) with this to remove dependency on BOOST.
// This uses the Gamma function as described here:
// http://mathworld.wolfram.com/ExponentialIntegral.html
#include <TMath.h>
#include <boost/math/special_functions/expint.hpp>
static double MyExpInt(double z)
{
	return boost::math::expint(1, z);
//	return TMath::Gamma(0,z);
}

CobremsGenerator::CobremsGenerator(double Emax_GeV, double Epeak_GeV)
{
   // Unique constructor for this class, initialize for the given
   // endpoint energy and peak position but these can be changed.

   fBeamEnergy = Emax_GeV;
   fBeamErms = 6.0e-4; // GeV
   fBeamEmittance = 2.5e-9; // m r
   fCollimatorSpotrms = 0.0005; // m
   fCollimatorDistance = 76.0; // m
   fCollimatorDiameter = 0.0034; // m
   fTargetThickness = 20e-6; // m
   fTargetThetay = 0.050; // radians
   fTargetThetaz = 0; // radians
   setTargetCrystal("diamond");
   setCoherentEdge(Epeak_GeV);
   fPhotonEnergyMin = 0.120; // GeV
   setPolarizedFlag(false);
   setCollimatedFlag(true);

#if COBREMS_GENERATOR_VERBOSITY > 0
   std::cout << std::endl
             << "Initialization for coherent bremsstralung calculation"
             << std::endl;
   printBeamlineInfo();
#endif
}

void CobremsGenerator::setTargetCrystal(std::string crystal)
{
   // declare the radiator target crystal type by name

   if (crystal == "diamond") {
      fTargetCrystal.name = "diamond";
      fTargetCrystal.Z = 6;
      fTargetCrystal.A = 12.01;
      fTargetCrystal.density = 3.534; // g/cm^3
      fTargetCrystal.lattice_constant = 3.56e-10; // m
      fTargetCrystal.Debye_Waller_const = 0.40e9; // 1/GeV^2
   }
   else if (crystal == "silicon") {
      fTargetCrystal.name = "silicon";
      fTargetCrystal.Z = 14;
      fTargetCrystal.A = 28.09;
      fTargetCrystal.density = 2.320; // g/cm^3
      fTargetCrystal.lattice_constant = 5.431e-10; // m
      fTargetCrystal.Debye_Waller_const = 1.5e9; // 1/GeV^2
   }
   else {
      std::cerr << "Error in CobremsGenerator::setTargetCrystal - "
                << "unknown crystal " << crystal << " requested, "
                << "cannot continue." << std::endl;
      exit(1);
   }

   // define the stanard unit cell of the diamond Bravais lattice
   fTargetCrystal.nsites = 8;
   fTargetCrystal.ucell_site.clear();
   fTargetCrystal.ucell_site.push_back(lattice_vector(0.0, 0.0, 0.0));
   fTargetCrystal.ucell_site.push_back(lattice_vector(0.0, 0.5, 0.5));
   fTargetCrystal.ucell_site.push_back(lattice_vector(0.5, 0.0, 0.5));
   fTargetCrystal.ucell_site.push_back(lattice_vector(0.5, 0.5, 0.0));
   fTargetCrystal.ucell_site.push_back(lattice_vector(0.25, 0.25, 0.25));
   fTargetCrystal.ucell_site.push_back(lattice_vector(0.25, 0.75, 0.75));
   fTargetCrystal.ucell_site.push_back(lattice_vector(0.75, 0.25, 0.75));
   fTargetCrystal.ucell_site.push_back(lattice_vector(0.75, 0.75, 0.25));
   fTargetCrystal.primaryHKL = lattice_vector(2,2,0);

   // approximate formula for atomic form factor beta
   fTargetCrystal.betaFF = 111 * pow(fTargetCrystal.Z, -1/3.) / me;

   // set mosaic spread to GlueX specification
   fTargetCrystal.mosaic_spread = 20e-6;

   // compute the radiation length
   fTargetCrystal.radiation_length = getTargetRadiationLength_Schiff();
}

double CobremsGenerator::getTargetDebyeWallerConstant(double DebyeT_K, 
                                                      double T_K)
{
   // Computes the Debye-Waller constant A for a simple model
   // assuming an isotropic crystal -- see Kaune et.al.
   //
   //       A(T) = A0 f(T)
   //  where
   //       A0 = 3 / (4 * atomicMass_GeV * DebyeTemperature_GeV)
   //  and
   //       f(T) = (2 / DebyeTemperature_GeV^2) *
   //              Integral_dw[0,DebyeTemperature_GeV] 
   //              {w (1 + 2 / (exp(w/T) - 1))}
   //
   //  T is the crystal temperature in GeV and A0 is the limiting
   //  value of A as T->0.

   double kBoltzmann = 8.617e-14; // GeV/K
   double amassGeV = fTargetCrystal.A * 0.932; // GeV
   double A0 = 3 / (2 * amassGeV * kBoltzmann * DebyeT_K); // /GeV^2
   double Tnormal = (T_K + 0.1) / DebyeT_K;
   int niter = 50;
   double f = 0;
   for (int iter=0; iter < niter; ++iter) {
      double x = (iter + 0.5) / niter;
      f += x * (1 + 2 / (exp(x / Tnormal) - 1)) / niter;
   }
   return A0 * f;
}

void CobremsGenerator::printBeamlineInfo()
{
   // Print a summary of the target crystal model parameters

   std::cout << " electron beam energy: " << fBeamEnergy << " GeV"
             << std::endl
             << " electron beam emittance: " 
             << fBeamEmittance * 1e9 << " mm.urad"
             << std::endl
             << " radiator crystal: " << fTargetCrystal.name
             << ", thickness " << fTargetThickness * 1e6 << " um"
             << std::endl
             << " radiation length: "
             << fTargetCrystal.radiation_length * 100 << " cm,"
             << " mosaic spread: " 
             << fTargetCrystal.mosaic_spread * 1e6 << " urad"
             << std::endl
             << " photon beam collimator half-angle: "
             << fCollimatorDiameter / (2 * fCollimatorDistance)
                * fBeamEnergy / me << " (m/E)"
             << std::endl
             << " collimator diameter: " 
             << fCollimatorDiameter * 100 << " cm"
             << std::endl
             << " crystal orientation: theta_x " 
             << fTargetThetax * 1e3 << " mrad"
             << std::endl
             << "                      theta_y " 
             << fTargetThetay * 1e3 << " mrad"
             << std::endl << std::endl;
}

void CobremsGenerator::printTargetCrystalInfo()
{
   // Print a summary of the target crystal model parameters

   double kBoltzmann = 8.617e-14; // GeV/K
   double amass = fTargetCrystal.A * 0.932; // GeV
   double DebyeTheta0 = 3 / (4 * amass * fTargetCrystal.Debye_Waller_const);
   double DebyeTheta300 = DebyeTheta0;
   for (int i=0; i < 5; i++) {
      DebyeTheta300 = DebyeTheta0 * (1 + 
                      pow(2 * dpi * 300 * kBoltzmann / DebyeTheta300, 2) / 6);
   }

   std::cout << "CobremsGenerator crystal type is " << fTargetCrystal.name
             << std::endl
             << "  atomic number Z=" << fTargetCrystal.Z
             << ", atomic weight A=" << fTargetCrystal.A << " amu"
             << std::endl
             << "  mass density: " << fTargetCrystal.density << " g/cm^3"
             << std::endl
             << "  radiation length: " << fTargetCrystal.radiation_length * 100
             << " cm" << std::endl
             << "  Debye-Waller constant: " << fTargetCrystal.Debye_Waller_const
             << " /GeV^2 (" << DebyeTheta300 / kBoltzmann << " K)"
             << std::endl
             << "  mosaic spread: " << fTargetCrystal.mosaic_spread * 1e6
             << " urad" << std::endl
             << "  atomic form-factor cutoff momentum: "
             << sqrt(1 / fTargetCrystal.betaFF) * 1e6 << " keV"
             << std::endl
             << "  primary lattice reflection h,k,l="
             << fTargetCrystal.primaryHKL.x << ","
             << fTargetCrystal.primaryHKL.y << ","
             << fTargetCrystal.primaryHKL.z
             << std::endl
             << "  lattice constant: " << fTargetCrystal.lattice_constant * 1e9
             << " nm" << std::endl
             << "  occupied sites of the crystal lattice unit cell are:"
             << std::endl;
   for (unsigned int i=0; i < fTargetCrystal.ucell_site.size(); ++i) {
      char s[100];
      snprintf(s, 100, "%4.2f %4.2f %4.2f", fTargetCrystal.ucell_site[i].x,
                                            fTargetCrystal.ucell_site[i].y,
                                            fTargetCrystal.ucell_site[i].z);
      std::cout << "    " << i + 1 << ":  " << s << std::endl;
   }
   std::cout << "  Crystal orientation matrix is:" << std::endl;
   for (int i=0; i < 3; ++i) {
      char s[100];
      snprintf(s, 100, "%15.12f %15.12f %15.12f", fTargetRmatrix[i][0],
                                                  fTargetRmatrix[i][1],
                                                  fTargetRmatrix[i][2]);
      std::cout << "    " << s << std::endl;
   }
}

void CobremsGenerator::applyBeamCrystalConvolution(int nbins, double *xvalues,
                                                              double *yvalues)
{
   // Electron beam emittance produces two effects in the coherent
   // bremsstrahlung spectrum:
   //   1) smears out the collimation acceptance function in production
   //      angle, so it varies smoothly to zero instead of being a step;
   //   2) combines with mosaic spread of the target crystal to smear out
   //      the relation between photon energy fraction x and production
   //      angle theta for a given lattice reflection.
   // The first one affects the left-hand (low energy) side of the
   // coherent peaks in the coherent bremsstrahlung spectrum, while the
   // second affects the right-hand (high energy) side, limiting the
   // sharpness of the edges in either case. 
   //
   // Effect (1) is taken into account in the way the acceptance function
   // is computed on the final state, but effect (2) is more difficult to
   // treat analytically because it involves smearing directions in the
   // initial state. This is only relevant to the coherent part of the
   // spectrum, where it acts by broadening the step on the high side of
   // the coherent edge from a sharp drop to a gradually sloping curve.
   // One way to take this into account in an effective manner is to treat
   // the coherent spectrum at every photon energy bin as being dominated
   // by a single reciprocal lattice vector, and smearing out the relation
   // x=x(theta) that follows from the two-body nature of the scattering
   // from that plane. Considering a 1D spectrum of beam intensity vs x,
   // this leads to a convolution of the distribution with an x-dependent
   // smearing function. The applyBeamConvolution method computes that
   // smearing function for each value of x in the input xvalues array
   // and applies it to the input spectrum represented by the yvalues
   // array. The yvalues array is overwritten with the convoluted spectrum.
   // For simplicity, the xvalues are assumed to be equally spaced.

   double x0 = xvalues[0];
   double x1 = xvalues[nbins - 1];
   double var0 = pow(fTargetCrystal.mosaic_spread, 2) +
                 pow(fBeamEmittance / fCollimatorSpotrms, 2);
   double varMS = Sigma2MS(fTargetThickness);

   // Here we have to guess which characteristic angle alph inside the crystal
   // is dominantly responsible for the coherent photons in each bin in x.
   // I just use the smallest of the two angles, but this does not work when
   // both angles are small, and you have to be more clever -- BEWARE!!!
   double alph = (fabs(fTargetThetax) < fabs(fTargetThetay))?
                  fabs(fTargetThetax) : fabs(fTargetThetay);
   if (alph == 0) {
      alph = (fabs(fTargetThetax) > fabs(fTargetThetay))?
              fabs(fTargetThetax) : fabs(fTargetThetay);
   }

   // In any case, fine-tuning below the mosaic spread limit makes no sense.
   else {
      alph = (alph > fTargetCrystal.mosaic_spread)?
              alph : fTargetCrystal.mosaic_spread;
   }

   double *norm = new double[nbins];
   double *result = new double[nbins];
   for (int j=0; j < nbins; ++j) {
      norm[j] = 0;
      result[j] = 0;
      for (int i=-nbins; i <= nbins; ++i) {
         double dx = (x1 - x0) * (j - i) / nbins;
         double x = x0 + (x1 - x0) * (j + 0.5) / nbins;
         double dalph = dx * alph / (x * (1 - x));
         double term;
         if (varMS / var0 > 1e-4) {
            term = dalph / varMS *
                          (erf(dalph / sqrt(2 * (var0 + varMS))) -
                           erf(dalph / sqrt(2 * var0))) +
                   sqrt(2 / dpi) / varMS *
                          (exp(-dalph*dalph / (2 * (var0 + varMS))) *
                                               sqrt(var0 + varMS) -
                           exp(-dalph*dalph / (2 * var0)) * sqrt(var0));
         }
         else {
            term = exp(-dalph*dalph / (2 * var0)) / sqrt(2 * dpi * var0);
         }
         term *= alph / x;
         norm[j] += term;
      }
   }

   for (int i=-nbins; i <= nbins; ++i) {
      int ii = abs(i);
      for (int j=0; j < nbins; ++j) {
         double dx = (x1 - x0) * (j - i) / nbins;
         double x = x0 + (x1 - x0) * (j + 0.5) / nbins;
         double dalph = dx * alph / (x * (1 - x));
         double term;
         if (varMS / var0 > 1e-4) {
            term = dalph / varMS *
                          (erf(dalph / sqrt(2 * (var0 + varMS))) -
                           erf(dalph / sqrt(2 * var0))) +
                   sqrt(2 / dpi) / varMS *
                          (exp(-dalph*dalph / (2 * (var0 + varMS))) *
                                                    sqrt(var0 + varMS) -
                           exp(-dalph*dalph/ (2 * var0)) * sqrt(var0));
         }
         else {
            term = exp(-dalph*dalph / (2 * var0)) / sqrt(2 * dpi * var0);
         }
         term *= alph / x;
         result[ii] += term * yvalues[j] / norm[j];
      }
   }

   for (int i=0; i < nbins; ++i) {
      if (fabs(result[i]) > 1e-35) {
         yvalues[i] = result[i];
      }
      else {
         yvalues[i] = 0;
      }
   }

   delete [] norm;
   delete [] result;
}

double CobremsGenerator::getTargetRadiationLength_PDG()
{
   // PDG formula for radiation length, converted to meters

   double Z = fTargetCrystal.Z;
   double N = fTargetCrystal.nsites;
   double a = fTargetCrystal.lattice_constant;
   double c = alpha * Z;
   double s = 4 * N * pow(alpha, 3) * pow(hbarc/(a*me), 2) / a *
              (Z*Z * (log(184.15 * pow(Z, -1/3.)) -
                      c*c * (1 / (1 + c*c) + 0.20206 - 0.0369 * c*c +
                             0.0083 * pow(c, 4) - 0.002 * pow(c, 6))) +
               Z * log(1194 * pow(Z, -2/3.)));
   return 1/s;
}

double CobremsGenerator::getTargetRadiationLength_Schiff()
{
   // Schiff formula for radiation length, converted to meters
 
   double Z = fTargetCrystal.Z;
   double N = fTargetCrystal.nsites;
   double a = fTargetCrystal.lattice_constant;
   double zeta = log(1440 * pow(Z, -2/3.)) / log(183 * pow(Z, -1/3.));
   double s = 4 * N * pow(alpha, 3) * pow(hbarc/(a*me), 2) / a *
              Z * (Z + zeta) * log(183 * pow(Z, -1/3.));
   return 1/s;
}

void CobremsGenerator::setCoherentEdge(double Epeak_GeV)
{
   // Adjust theta_x of the target to align the coherent edge at
   // energy Epeak_GeV in the photon spectrum, then orient the
   // crystal according to theta_x, theta_y, theta_z tip angles.

   double edge = Epeak_GeV;
   double qtotal = hbarc * (2 * dpi / fTargetCrystal.lattice_constant);
   lattice_vector hkl = fTargetCrystal.primaryHKL;
   qtotal *= sqrt(hkl.x * hkl.x + hkl.y * hkl.y + hkl.z * hkl.z);
   double qlong = edge * me*me / (2 * fBeamEnergy * (fBeamEnergy - edge));
   fTargetThetax = -qlong / qtotal;
   resetTargetOrientation();
   RotateTarget(0, dpi/2, 0);      // point (1,0,0) along beam
   RotateTarget(0, 0, dpi/4);      // point (0,1,1) vertically
   RotateTarget(-fTargetThetax, 0, 0);
   RotateTarget(0, -fTargetThetay, 0);
   RotateTarget(0, 0, -fTargetThetaz);
}

CobremsGenerator::CobremsGenerator(const CobremsGenerator &src)
{
   // copy constructor

   fTargetCrystal = src.fTargetCrystal;
   fTargetThickness = src.fTargetThickness;
   fTargetThetax = src.fTargetThetax;
   fTargetThetay = src.fTargetThetay;
   fTargetThetaz = src.fTargetThetaz;
   for (int i=0; i < 3; ++i)
      for (int j=0; j < 3; ++j)
         fTargetRmatrix[i][j] = src.fTargetRmatrix[i][j];
   fBeamEnergy = src.fBeamEnergy;
   fBeamErms = src.fBeamErms;
   fBeamEmittance = src.fBeamEmittance;
   fCollimatorSpotrms = src.fCollimatorSpotrms;
   fCollimatorDistance = src.fCollimatorDistance;
   fCollimatorDiameter = src.fCollimatorDiameter;
   fQ2theta2 = src.fQ2theta2;
   fQ2weight = src.fQ2weight;
}

CobremsGenerator &CobremsGenerator::operator=(const CobremsGenerator &src)
{
   // assignment operator

   fTargetCrystal = src.fTargetCrystal;
   fTargetThickness = src.fTargetThickness;
   fTargetThetax = src.fTargetThetax;
   fTargetThetay = src.fTargetThetay;
   fTargetThetaz = src.fTargetThetaz;
   for (int i=0; i < 3; ++i)
      for (int j=0; j < 3; ++j)
         fTargetRmatrix[i][j] = src.fTargetRmatrix[i][j];
   fBeamEnergy = src.fBeamEnergy;
   fBeamErms = src.fBeamErms;
   fBeamEmittance = src.fBeamEmittance;
   fCollimatorSpotrms = src.fCollimatorSpotrms;
   fCollimatorDistance = src.fCollimatorDistance;
   fCollimatorDiameter = src.fCollimatorDiameter;
   fQ2theta2 = src.fQ2theta2;
   fQ2weight = src.fQ2weight;
   return *this;
}

CobremsGenerator::~CobremsGenerator() { }

double CobremsGenerator::CoherentEnhancement(double x)
{
   // Returns ratio of total bremsstrahlung yield over incoherent yield
   // for photon energy k = x*fBeamEnergy

   double yc = Rate_dNcdx(x);
   double yi = Rate_dNidx(x);
   return (yi + yc) / (yi + 1e-99);
}

double CobremsGenerator::Rate_dNtdx(double x)
{
   // Returns total bremsstrahlung probability density differential in
   // x (scaled photon energy) at photon energy k = x*fBeamEnergy.

   return Rate_dNcdx(x) + Rate_dNidx(x);
}

double CobremsGenerator::Rate_dNtdx(double x, 
                                    double distance_m, double diameter_m)
{
   // Returns total bremsstrahlung probability density differential in x
   // (scaled photon energy) at photon energy k = x*fBeamEnergy with
   // user-specified variations in the collimator distance and diameter,
   // for plotting. Special case: if diameter_m < 0 then interpret its
   // absolute value as the collimator radius in characteristic units m/E.

   double dist = fCollimatorDistance;
   double diam = fCollimatorDiameter;
   fCollimatorDistance = (distance_m > 0)? distance_m : fCollimatorDistance;
   fCollimatorDiameter = (diameter_m > 0)? diameter_m : (diameter_m < 0)?
                         -2 * distance_m * diameter_m * me / fBeamEnergy :
                         fCollimatorDiameter;
   double rate = Rate_dNtdx(x);
   fCollimatorDistance = dist;
   fCollimatorDiameter = diam;
   return rate;
}

double CobremsGenerator::Rate_dNtdk(double k_GeV)
{
   // Returns total bremsstrahlung probability density differential
   // in photon energy k (GeV).

   return Rate_dNtdx(k_GeV / fBeamEnergy) / fBeamEnergy;
}

double CobremsGenerator::Rate_dNcdx(double x)
{
   // Returns the coherent bremsstrahlung probability density differential
   // in x (scaled photon energy) at photon energy k = x*fBeamEnergy.

   return 2 * dpi * Rate_dNcdxdp(x, dpi/4);
}

double CobremsGenerator::Rate_dNcdx(double x, 
                                    double distance_m, double diameter_m)
{
   // Returns the coherent bremsstrahlung probability density differential
   // in x (scaled photon energy) at photon energy k = x*fBeamEnergy with
   // user-specified variations in the collimator distance and diameter.
   // Special case: if diameter_m < 0 then interpret its absolute
   // value as the collimator radius in characteristic units m/E.

   double dist = fCollimatorDistance;
   double diam = fCollimatorDiameter;
   fCollimatorDistance = (distance_m > 0)? distance_m : fCollimatorDistance;
   fCollimatorDiameter = (diameter_m > 0)? diameter_m : (diameter_m < 0)?
                         -2 * distance_m * diameter_m * me / fBeamEnergy :
                         fCollimatorDiameter;
   double rate = 2 * dpi * Rate_dNcdxdp(x, dpi/4);
   fCollimatorDistance = dist;
   fCollimatorDiameter = diam;
   return rate;
}

double CobremsGenerator::Rate_dNcdxdp(double x, double phi)
{
   // Returns the coherent bremsstrahlung probabililty density differential
   // in x (scaled photon energy) and phi (azimuthal emission angle) for
   // fixed photon energy k = x*fBeamEnergy and phi. If fPolarizedFlag is
   // false (0, default) then the total yield is returned, otherwise it is
   // only the polarized fraction. If fCollimatedFlag is false (0) then
   // the total yield is returned, otherwise only the part that passes the
   // collimator is counted (default).

   double Z = fTargetCrystal.Z;
   double a = fTargetCrystal.lattice_constant;
   double sigma0 = 16 * dpi * fTargetThickness * Z*Z * pow(alpha, 3) *
                   fBeamEnergy * hbarc/(a*a) * pow(hbarc / (a * me), 4);

   fQ2theta2.clear();
   fQ2weight.clear();
   double qzmin = 99;
   int hmin, kmin, lmin;
   double sum = 0;
   // can restrict to h=0 for cpu speedup, if crystal alignment is "reasonable"
   for (int h = -4; h <= 4; ++h) {
      for (int k = -10; k <= 10; ++k) {
         for (int l = -10; l <= 10; ++l) {
            if (h/2 * 2 == h) {
               if (k/2 * 2 != k || l/2 * 2 != l ||
                   (h + k + l)/4 * 4 != h + k + l)
               {
                  continue;
               }
            }
            else if (k/2 * 2 == k || l/2 * 2 == l) {
              continue;
            }
            double ReS = 0;
            double ImS = 0;
            for (int i=0; i < fTargetCrystal.nsites; ++i) {
              double qdota = 2 * dpi * (h * fTargetCrystal.ucell_site[i].x +
                                        k * fTargetCrystal.ucell_site[i].y +
                                        l * fTargetCrystal.ucell_site[i].z);
              ReS += cos(qdota);
              ImS += sin(qdota);
            }
            double S2 = ReS*ReS + ImS*ImS;
            if (S2 < 1e-4)
               continue;
            double qnorm = hbarc * 2 * dpi / a;
            double q[3];
            q[0] = qnorm * (fTargetRmatrix[0][0] * h +
                            fTargetRmatrix[0][1] * k +
                            fTargetRmatrix[0][2] * l);
            q[1] = qnorm * (fTargetRmatrix[1][0] * h +
                            fTargetRmatrix[1][1] * k +
                            fTargetRmatrix[1][2] * l);
            q[2] = qnorm * (fTargetRmatrix[2][0] * h +
                            fTargetRmatrix[2][1] * k +
                            fTargetRmatrix[2][2] * l);
            double q2 = q[0]*q[0] + q[1]*q[1] + q[2]*q[2];
            double qT2 = q[0]*q[0] + q[1]*q[1];
            double xmax = 2 * fBeamEnergy * q[2];
            xmax /= xmax + me*me;
            if (x > xmax || xmax > 1) {
               continue;
            }

#if COBREMS_GENERATOR_VERBOSITY > 2
            else {
               std::cout << h << "," << k << "," << l << ","
                         << S2 << "," << q2 << "," << xmax
                         << std::endl;
            }
#endif

            if (q[2] < qzmin) {
               qzmin = q[2];
               hmin = h;
               kmin = k;
               lmin = l;
            }
            double theta2 = (1 - x) * xmax / (x * (1 - xmax)) - 1;
            double betaFF2 = pow(fTargetCrystal.betaFF, 2);
            double FF = 1 / (1 + q2 * betaFF2);
            sum += sigma0 * qT2 * S2 * pow(FF * betaFF2, 2) *
                   exp(-q2 * fTargetCrystal.Debye_Waller_const) *
                   ((1 - x) / pow(x * (1 + theta2), 2)) *
                   ((1 + pow(1 - x, 2)) - 8 * (theta2 / pow(1 + theta2, 2) * 
                                              (1 - x) * pow(cos(phi), 2))) *
                   ((fCollimatedFlag)? Acceptance(theta2) : 1) *
                   ((fPolarizedFlag)? Polarization(x, theta2) : 1);
            fQ2theta2.push_back(theta2);
            fQ2weight.push_back(sum);
         }
      }
   }

#if COBREMS_GENERATOR_VERBOSITY > 1
   if (qzmin < 99) {
      std::cout << hmin << "," << kmin << "," << lmin
                << " is the best plane at x=" << x
                << std::endl;
   }
#endif

   return sum;
}

double CobremsGenerator::Rate_dNidx(double x)
{
   // Returns the incoherent bremsstrahlung probabililty density differential
   // in x (scaled photon energy) at fixed photon energy k = x*fBeamEnergy.
 
   if (x > 1)
      return 0;

   // Numerical integration in d(theta**2) over [0,inf]
   // is mapped onto u=1/(1+theta^2) as (1/u^2) d(u) over [0,1]
   int niter = 50;
   double dNidx = 0;
   double du = 1. / niter;
   for (int iter = 0; iter < niter; ++iter) {
      double u = (iter + 0.5) / niter;
      double theta2 = (1 - u) / u;
      dNidx += Rate_dNidxdt2(x, theta2) * du/(u*u);
   }
   return dNidx;
}

double CobremsGenerator::Rate_dNBidx(double x)
{
   // In the following paper, a closed form is given for the integral that
   // is being performed analytically by dNidx.  I include this second form
   // here in case some time it might be useful as a cross check.
   //
   // "Coherent bremsstrahlung in crystals as a tool for producing high
   // energy photon beams to be used in photoproduction experiments at
   // CERN SPS", Nucl. Instr. Meth. 204 (1983) pp.299-310. 
   //
   // Note: in this paper they have swapped subscripts for coherent and
   // incoherent intensities.  This is not very helpful to the reader!
   //
   // The result is some 15% lower radiation rate than the result of dNidx.
   // I take the latter to be more detailed (because it gives a more
   // realistic behaviour at the endpoint and agrees better with the PDG
   // radiation length for carbon).  Most of this deficiency is remedied
   // by simply replacing Z**2 in the cross section with Z*(Z+zeta) as
   // recommended by Kaune et.al., and followed by the PDG in their fit
   // to radiation lengths.
   //
   //                            WARNING
   // dNidx and dNBidx give the incoherent radiation rate for crystalline
   // radiators.  If you take the incoherent radiation formulae here and
   // integrate them you will NOT obtain the radiation length for amorphous
   // radiators; it will be overestimated by some 15%.  The reason is that
   // the part of the integral in q-space that is covered by the discrete
   // sum has been subtracted to avoid double-counting with the coherent
   // part.  If you were to spin the crystal fast enough, the coherent
   // spectrum should average out to yield the remaining 15% with a 
   // spectral shape resembling the Bethe-Heitler result.

   double Z = fTargetCrystal.Z;
   double betaFF = fTargetCrystal.betaFF;
   double a = fTargetCrystal.lattice_constant;
   double AoverB2 = fTargetCrystal.Debye_Waller_const / (betaFF * betaFF);
   double Tfact = -(1 + AoverB2) * exp(AoverB2) *
                  MyExpInt( AoverB2);
   double psiC1 = 2 * (2 * log(betaFF * me) + Tfact + 2);
   double psiC2 = psiC1 - 2/3.;
   double zeta = log(1440 * pow(Z, -2/3.)) / log(183 * pow(Z, -1/3.));
   double dNBidx = fTargetCrystal.nsites * fTargetThickness *
                   Z * (Z + zeta) * pow(alpha, 3) * 
                   pow(hbarc / (a*me), 2) / (a * x) *
                   (psiC1 * (1 + pow(1 - x, 2)) - psiC2 * (1 - x) * 2/3.);
   return dNBidx;
}

double CobremsGenerator::Rate_dNidxdt2(double x, double theta2)
{
   // Returns the incoherent bremsstrahlung probabililty density differential
   // in x (scaled photon energy) and theta^2 at fixed photon energy 
   // k = x*fBeamEnergy and production angle theta. Argument theta2 is equal
   // to theta^2 expressed in units of (me/fBeamEnergy)^2. If internal flag
   // fCollimatedFlag is false (0) then the total yield is returned,
   // otherwise only the part that passes the collimator is counted (default).
 
   double delta = 1.02;
   double Z = fTargetCrystal.Z;
   double betaFF = fTargetCrystal.betaFF;
   double a = fTargetCrystal.lattice_constant;
   double zeta = log(1440 * pow(Z, -2/3.)) / log(183 * pow(Z, -1/3.));
   double MSchiff = 1 / (pow((me * x) / (2 * fBeamEnergy * (1 - x)), 2) +
                         1 / pow(betaFF * me * (1 + theta2), 2));
   double dNidxdt2 = 2 * fTargetCrystal.nsites * fTargetThickness * Z *
                     (Z + zeta) * pow(alpha, 3) * pow(hbarc/(a*me), 2) / (a*x) *
                     ( ((1 + pow(1 - x, 2)) - 4 * theta2 * (1 - x) / 
                        pow(1 + theta2, 2)) /
                       pow(1 + theta2, 2) *
                       (log(MSchiff) - 2 * delta * Z / (Z + zeta)) +
                        16 * theta2 * (1 - x) / pow(1 + theta2, 4) -
                        pow(2 - x, 2) / pow(1 + theta2, 2) ) *
                   ((fCollimatedFlag)? Acceptance(theta2) : 1);

_DBG_<<"fTargetCrystal.nsites="<<fTargetCrystal.nsites<<" fTargetThickness="<<fTargetThickness<<" Z="<<Z<<endl;
_DBG_<<"line2="<<(Z + zeta) * pow(alpha, 3) * pow(hbarc/(a*me), 2) / (a*x)<<endl;
_DBG_<<"line3="<< (1 + pow(1 - x, 2)) - 4 * theta2 * (1 - x)<<endl;
_DBG_<<"line4="<<pow(1 + theta2, 2)<<endl;
_DBG_<<"line6="<<(log(MSchiff) - 2 * delta * Z / (Z + zeta))<<endl;
_DBG_<<"line7="<<16 * theta2 * (1 - x) / pow(1 + theta2, 4)<<endl;
_DBG_<<"line8="<<pow(2 - x, 2) / pow(1 + theta2, 2)<<endl;
_DBG_<<"line9="<<((fCollimatedFlag)? Acceptance(theta2) : 1)<<endl;

   return dNidxdt2;
}

double CobremsGenerator::Rate_para(double x, double theta2, double phi)
{
   // Returns the relative rate of in-plane polarized flux from coherent
   // bremsstrahlung at production angles theta and phi and photon energy
   // k = x*fBeamEnergy. The units are arbitrary, but the same as Rate_ortho
   // (see below). The argument theta2 is the production polar angle theta^2
   // expressed in units of (me/fBeamEnergy)^2.

   return 0.5 * pow((2 - x) * (1 + theta2), 2) -
          8 * theta2 * (1 - x) * pow(cos(phi), 2) -
          8 * pow(theta2, 2) * (1 - x) * pow(cos(phi), 2) * pow(sin(phi), 2);
}

double CobremsGenerator::Rate_ortho(double x, double theta2, double phi)
{
   // Returns the relative rate of out-of-plane polarized flux from coherent
   // bremsstrahlung at production angles theta and phi and photon energy k
   // = x*fBeamEnergy. The units are arbitrary, but the same as Rate_para
   // (see above). The argument theta2 is the production polar angle theta^2
   // expressed in units of (me/fBeamEnergy)^2.

   return 0.5 * pow(x * (1 + theta2), 2) +
          8 * pow(theta2, 2) * (1 - x) * pow(cos(phi), 2) * pow(sin(phi), 2);
}

double CobremsGenerator::Polarization(double x, double theta2)
{
   // Returns the degree of linear polarization in a coherent bremsstrahlung
   // beam at photon energy k = x*fBeamEnergy and production angle theta.
   // The formula evaluated below is the azimuthal average of the ratio
   //       (Rate_para - Rate_ortho) / (Rate_para + Rate_ortho)
   // The argument theta2 is the production polar angle theta^2 expressed
   // in units of (me/fBeamEnergy)^2.

   return 2 * (1 - x) / (pow(1 + theta2, 2) * (pow(1 - x, 2) + 1) - 
                         4 * theta2 * (1 - x));
}

double CobremsGenerator::Acceptance(double theta2, double phi,
                                    double xshift_m, double yshift_m)
{
   // Returns the acceptance of the collimator for photons emitted at
   // polar angle theta and azimuthal angle phi at the radiator. Both
   // beam emittance and multiple-scattering in the target contribute
   // to smearing of the angular acceptance at the the collimator edge.
   // The argument theta2 is the production polar angle theta^2
   // expressed in units of (me/fBeamEnergy)^2. Misalignment of the
   // collimator with the beam axis is taken into account by the 
   // arguments xshift,yshift.

   double theta = sqrt(theta2) * (me/fBeamEnergy);
   double xc = fCollimatorDistance * tan(theta) * cos(phi) + xshift_m;
   double yc = fCollimatorDistance * tan(theta) * sin(phi) + yshift_m;
   double thetaprime = atan2(sqrt(xc*xc + yc*yc), fCollimatorDistance);
   return Acceptance(pow(thetaprime * fBeamEnergy/me, 2));
}

double CobremsGenerator::Acceptance(double theta2)
{
   // Returns the acceptance of the collimator for photons emitted at
   // polar angle theta and azimuthal angle phi at the radiator. Both
   // beam emittance and multiple-scattering in the target contribute
   // to smearing of the angular acceptance at the the collimator edge.
   // The argument theta2 is the production polar angle theta^2
   // expressed in units of (me/fBeamEnergy)^2.

   double acceptance = 0;
   double niter = 50;
   double theta = sqrt(theta2);
   double thetaC = fCollimatorDiameter / (2 * fCollimatorDistance) *
                                              fBeamEnergy / me;
   double var0 = pow((fCollimatorSpotrms / fCollimatorDistance) *
                                              fBeamEnergy / me, 2);
   double varMS = Sigma2MS(fTargetThickness) * pow(fBeamEnergy / me, 2);
   if (theta < thetaC) {
      double u1 = thetaC - theta;
      if (u1*u1 / (var0 + varMS) > 20) {
         return 1;
      }
      for (int iter = 0; iter < niter; ++iter) {
         double u = u1 * (iter + 0.5) / niter;
         double u2 = u * u;
         double du2 = 2 * u * u1 / niter;
         double pu;
         if (varMS / var0 > 1e-4) {
            pu = (MyExpInt( u2 / (2 * (var0 + varMS))) -
                  MyExpInt( u2 / (2 * var0))) / (2 * varMS);
         }
         else {
            pu = exp(-u2 / (2 * var0)) / (2 * var0);
         }
_DBG_<<"pu="<<pu<<" du2="<<du2<<endl;
         acceptance += pu * du2;
      }
   }
   double u0 = fabs(theta - thetaC);
   double u1 = fabs(theta + thetaC);
   for (int iter = 0; iter < niter; ++iter) {
      double u = u0 + (u1 - u0) * (iter + 0.5) / niter;
      double u2 = u * u;
      double du2 = 2 * u * (u1 - u0) / niter;
      double pu;
      if (varMS / var0 > 1e-4) {
         pu = (MyExpInt( u2 / (2 * (var0 + varMS))) -
               MyExpInt( u2 / (2 * var0))) / (2 * varMS);
      }
      else {
         pu = exp(-u2 / (2 * var0)) / (2 * var0);
      }
      acceptance += pu * du2/dpi *
                    atan2(sqrt((theta2 - pow(thetaC - u, 2)) *
                                         (pow(thetaC + u, 2) - theta2)),
                          theta2 - pow(thetaC, 2) + u2);
   }
   return acceptance;
}
   
void CobremsGenerator::RotateTarget(double thetax, 
                                    double thetay,
                                    double thetaz)
{
   // Apply a sequence of rotations to the target crystal as
   //   Rmatrix(out) = Rx(thx) Ry(thy) Rz(thz) Rmatrix(in)
   // with rotations understood in the passive sense.

   if (thetaz != 0) {
      double sint = sin(thetaz);
      double cost = cos(thetaz);
      for (int i=0; i < 3; ++i) {
         double x = fTargetRmatrix[0][i];
         double y = fTargetRmatrix[1][i];
         fTargetRmatrix[0][i] = cost * x + sint * y;
         fTargetRmatrix[1][i] = cost * y - sint * x;
      }
   }
   if (thetay != 0) {
      double sint = -sin(thetay);
      double cost = cos(thetay);
      for (int i=0; i < 3; ++i) {
         double x = fTargetRmatrix[0][i];
         double z = fTargetRmatrix[2][i];
         fTargetRmatrix[0][i] = cost * x + sint * z;
         fTargetRmatrix[2][i] = cost * z - sint * x;
      }
   }
   if (thetax != 0) {
      double sint = sin(thetax);
      double cost = cos(thetax);
      for (int i=0; i < 3; ++i) {
         double y = fTargetRmatrix[1][i];
         double z = fTargetRmatrix[2][i];
         fTargetRmatrix[1][i] = cost * y + sint * z;
         fTargetRmatrix[2][i] = cost * z - sint * y;
       }
   }
}

double CobremsGenerator::Sigma2MS(double thickness_m)
{
   // Returns the mean-square multiple-scattering angle of the
   // electron beam inside the radiator crystal target, in radians.
   // This method wraps one of the concrete implementations, see below.
   // Some formulas, although valid for a reasonable range of target
   // thickness, can go negative for extremely small target thicknesses.
   // Here I protect against these unusual cases by taking the absolute value.

   return fabs(Sigma2MS_Geant(thickness_m));
}

double CobremsGenerator::Sigma2MS_Kaune(double thickness_m)
{
   // Multiple scattering formula of Kaune et.al. 
   // with a correction factor from a multiple-scattering calculation
   // taking into account the atomic and nuclear form factors for carbon.
   //
   // Note by RTJ, Oct. 13, 2008:
   // I think this formula overestimates multiple scattering in thin targets
   // like these diamond radiators, because it scales simply like sqrt(t).
   // Although the leading behavior is sqrt(t/radlen), it should increase
   // faster than that because of the 1/theta^2 tail of the Rutherford
   // distribution that makes the central gaussian region swell with increasing
   // number of scattering events.  For comparison, I include below the PDG
   // formula (sigma2MS_PDG), the Moliere formula used in the Geant3 simulation
   // of gaussian multiple scattering (sigma2MS_Geant), and a Moliere fit for
   // thin targets taken from reference Phys.Rev. vol.3 no.2, (1958), p.647
   // (sigma2MS_Hanson).  The latter two separate the gaussian part from the
   // tails in different ways, but both agree that the central part is much
   // more narrow than the formulation by Kaune et.al. below.

   double carboncor = 4.2 / 4.6;
   double Z = fTargetCrystal.Z;
   double a = fTargetCrystal.lattice_constant;
   return 8 * dpi * fTargetCrystal.nsites * pow(alpha * Z, 2) *
          thickness_m * pow(hbarc / (fBeamEnergy * a), 2) / a *
          log(183 * pow(Z, -1/3.)) * 
          carboncor;
}

double CobremsGenerator::Sigma2MS_PDG(double thickness_m)
{
   // Evaluates the PDG formula for multiple scattering of the beam electron
   // inside the target crystal, with beta=1, charge=1. This formula is said
   // to be within 11% for t > 1e-3 rad.len.

   double t = thickness_m / fTargetCrystal.radiation_length;
   return pow(13.6e-3 / fBeamEnergy, 2) * t * pow(1 + 0.038 * log(t), 2);
}

double CobremsGenerator::Sigma2MS_Geant(double thickness_m)
{
   // Returns the Geant3 formula for the rms multiple-scattering angle
   // This formula is based on the theory of Moliere scattering.  It contains
   // a cutoff parameter F that is used for the fractional integral of the
   // scattering probability distribution that is included in computing the
   // rms.  This is needed because the complete distribution of scattering
   // angles connects smoothly from a central gaussian (small-angle
   // multiple-scattering regime) to a 1/theta^2 tail (large-angle Rutherford
   // scattering regime) through the so-called plural scattering region.

   double rBohr = 0.52917721e-10; // m
   double F = 0.98; // probability cutoff in definition of sigma2MS
   double Z = fTargetCrystal.Z;
   double chi2cc = pow(0.39612e-2, 2) * Z * (Z + 1) *
                   fTargetCrystal.density / 12; // GeV^2/m
   double chi2c = chi2cc * thickness_m / pow(fBeamEnergy, 2);
   double chi2alpha = 1.13 * pow(hbarc / (fBeamEnergy * rBohr * 0.885), 2) *
                      pow(Z, 2/3.) * (1 + 3.34 * pow(alpha * Z, 2));
   double omega0 = chi2c / (1.167 * chi2alpha); // mean number of scatters
   double gnu = omega0 / (2 * (1 - F));
   return chi2c / (1 + pow(F, 2)) * ((1 + gnu) / gnu * log(1 + gnu) -1);
}

double CobremsGenerator::Sigma2MS_Hanson(double thickness_m)
{
   // Formulation of the rms projected angle attributed to Hanson et.al.
   // in reference Phys.Rev. vol.3 no.2, (1958), p.647.  This is just Moliere
   // theory used to give the 1/e angular width of the scattering distribution.
   // In the paper, though, they compare it with experiment for a variety of
   // metal foils down to 1e-4 rad.len. in thickness, and show excellent
   // agreement with the gaussian approximation out to 4 sigma or so.  I
   // like this paper because of the excellent agreement between the theory
   // and experimental data.

   double Z = fTargetCrystal.Z;
   double ttingcm2 = thickness_m * 100 * fTargetCrystal.density;
   double EinMeV = fBeamEnergy * 1000;
   double theta2max = 0.157 * Z * (Z + 1) / fTargetCrystal.A *
                      ttingcm2 / pow(EinMeV, 2);
   double theta2screen = theta2max * fTargetCrystal.A * 
                         (1 + 3.35 * pow(Z * alpha, 2)) /
                         (7800 * (Z + 1) * pow(Z, 1/3.) * ttingcm2);
   double BminuslogB = log(theta2max / theta2screen) - 0.154;
   double Blast = 1;
   double B;
   for (int i=0; i < 999; ++i) {
      B = BminuslogB + log(Blast);
      if (B < 1.2) {
         B = 1.21;
         break;
      }
      else if (fabs(B - Blast) > 1e-6) {
         Blast = B;
      }
      else {
         break;
      }
   }
   return theta2max * (B - 1.2) / 2;
}

#ifdef BOOST_PYTHON_WRAPPING

void CobremsGenerator::pyApplyBeamCrystalConvolution(int nbins, pyobject xarr,
                                                                pyobject yarr)
{
   using boost::python::extract;
   typedef boost::python::tuple pytuple;
   pytuple xtuple = extract<pytuple>(xarr.attr("buffer_info")());
   pytuple ytuple = extract<pytuple>(yarr.attr("buffer_info")());
   double *xbuf = reinterpret_cast<double*>((int)extract<int>(xtuple[0]));
   double *ybuf = reinterpret_cast<double*>((int)extract<int>(ytuple[0]));
   applyBeamCrystalConvolution(nbins, xbuf, ybuf);
}

double (CobremsGenerator::*Rate_dNtdx_1)(double) = &CobremsGenerator::Rate_dNtdx;
double (CobremsGenerator::*Rate_dNtdx_3)(double, double, double) = &CobremsGenerator::Rate_dNtdx;
double (CobremsGenerator::*Rate_dNcdx_1)(double) = &CobremsGenerator::Rate_dNcdx;
double (CobremsGenerator::*Rate_dNcdx_3)(double, double, double) = &CobremsGenerator::Rate_dNcdx;
double (CobremsGenerator::*Acceptance_1)(double) = &CobremsGenerator::Acceptance;
double (CobremsGenerator::*Acceptance_4)(double, double, double, double) = &CobremsGenerator::Acceptance;

BOOST_PYTHON_MODULE(libcobrems)
{
   using boost::python::class_;
   using boost::python::enum_;
   using boost::python::def;

   class_<CobremsGenerator, CobremsGenerator*>
         ("CobremsGenerator", 
          "coherent bremsstrahlung spectrum and polarization calculator, "
          "with methods for generating random Monte Carlo samples",
          boost::python::init<double, double>())
      .def("setBeamEnergy", &CobremsGenerator::setBeamEnergy)
      .def("setBeamErms", &CobremsGenerator::setBeamErms)
      .def("setBeamEmittance", &CobremsGenerator::setBeamEmittance)
      .def("setCollimatorSpotrms", &CobremsGenerator::setCollimatorSpotrms)
      .def("setCollimatorDistance", &CobremsGenerator::setCollimatorDistance)
      .def("setCollimatorDiameter", &CobremsGenerator::setCollimatorDiameter)
      .def("setTargetThickness", &CobremsGenerator::setTargetThickness)
      .def("setTargetCrystal", &CobremsGenerator::setTargetCrystal)
      .def("resetTargetOrientation", &CobremsGenerator::resetTargetOrientation)
      .def("setCoherentEdge", &CobremsGenerator::setCoherentEdge)
      .def("setTargetThetax", &CobremsGenerator::setTargetThetax)
      .def("setTargetThetay", &CobremsGenerator::setTargetThetay)
      .def("setTargetThetaz", &CobremsGenerator::setTargetThetaz)
      .def("RotateTarget", &CobremsGenerator::RotateTarget)
      .def("getBeamEnergy", &CobremsGenerator::getBeamEnergy)
      .def("getBeamErms", &CobremsGenerator::getBeamErms)
      .def("getBeamEmittance", &CobremsGenerator::getBeamEmittance)
      .def("getCollimatorSpotrms", &CobremsGenerator::getCollimatorSpotrms)
      .def("getCollimatorDistance", &CobremsGenerator::getCollimatorDistance)
      .def("getCollimatorDiameter", &CobremsGenerator::getCollimatorDiameter)
      .def("getTargetThickness", &CobremsGenerator::getTargetThickness)
      .def("getTargetCrystal", &CobremsGenerator::getTargetCrystal)
      .def("getTargetCrystalNsites", &CobremsGenerator::getTargetCrystalNsites)
      .def("getTargetCrystalAtomicNumber", &CobremsGenerator::getTargetCrystalAtomicNumber)
      .def("getTargetCrystalAtomicWeight", &CobremsGenerator::getTargetCrystalAtomicWeight)
      .def("getTargetCrystalDensity", &CobremsGenerator::getTargetCrystalDensity)
      .def("getTargetCrystalLatticeConstant", &CobremsGenerator::getTargetCrystalLatticeConstant)
      .def("getTargetCrystalRadiationLength", &CobremsGenerator::getTargetCrystalRadiationLength)
      .def("getTargetCrystalDebyeWallerConst", &CobremsGenerator::getTargetCrystalDebyeWallerConst)
      .def("getTargetCrystalMosaicSpread", &CobremsGenerator::getTargetCrystalMosaicSpread)
      .def("getTargetCrystalBetaFF", &CobremsGenerator::getTargetCrystalBetaFF)
      .def("getTargetThetax", &CobremsGenerator::getTargetThetax)
      .def("getTargetThetay", &CobremsGenerator::getTargetThetay)
      .def("getTargetThetaz", &CobremsGenerator::getTargetThetaz)
      .def("getTargetRadiationLength_PDG", &CobremsGenerator::getTargetRadiationLength_PDG)
      .def("getTargetRadiationLength_Schiff", &CobremsGenerator::getTargetRadiationLength_Schiff)
      .def("getTargetDebyeWallerConstant", &CobremsGenerator::getTargetDebyeWallerConstant)
      .def("getCollimatedFlag", &CobremsGenerator::getCollimatedFlag)
      .def("setCollimatedFlag", &CobremsGenerator::setCollimatedFlag)
      .def("getPolarizedFlag", &CobremsGenerator::getPolarizedFlag)
      .def("setPolarizedFlag", &CobremsGenerator::setPolarizedFlag)
      .def("applyBeamCrystalConvolution", &CobremsGenerator::pyApplyBeamCrystalConvolution)
      .def("printBeamlineInfo", &CobremsGenerator::printBeamlineInfo)
      .def("printTargetCrystalInfo", &CobremsGenerator::printTargetCrystalInfo)
      .def("CoherentEnhancement", &CobremsGenerator::CoherentEnhancement)
      .def("Rate_dNtdx", Rate_dNtdx_1)
      .def("Rate_dNtdx", Rate_dNtdx_3)
      .def("Rate_dNtdk", &CobremsGenerator::Rate_dNtdk)
      .def("Rate_dNcdx", Rate_dNcdx_1)
      .def("Rate_dNcdx", Rate_dNcdx_3)
      .def("Rate_dNcdxdp", &CobremsGenerator::Rate_dNcdxdp)
      .def("Rate_dNidx", &CobremsGenerator::Rate_dNidx)
      .def("Rate_dNBidx", &CobremsGenerator::Rate_dNBidx)
      .def("Rate_dNidxdt2", &CobremsGenerator::Rate_dNidxdt2)
      .def("Rate_para", &CobremsGenerator::Rate_para)
      .def("Rate_ortho", &CobremsGenerator::Rate_ortho)
      .def("Polarization", &CobremsGenerator::Polarization)
      .def("Acceptance", Acceptance_1)
      .def("Acceptance", Acceptance_4)
      .def("Sigma2MS", &CobremsGenerator::Sigma2MS)
      .def("Sigma2MS_Kaune", &CobremsGenerator::Sigma2MS_Kaune)
      .def("Sigma2MS_PDG", &CobremsGenerator::Sigma2MS_PDG)
      .def("Sigma2MS_Geant", &CobremsGenerator::Sigma2MS_Geant)
      .def("Sigma2MS_Hanson", &CobremsGenerator::Sigma2MS_Hanson)
      .def_readonly("dpi", &CobremsGenerator::dpi)
      .def_readonly("me", &CobremsGenerator::me)
      .def_readonly("alpha", &CobremsGenerator::alpha)
      .def_readonly("hbarc", &CobremsGenerator::hbarc)
   ;
}

#endif
