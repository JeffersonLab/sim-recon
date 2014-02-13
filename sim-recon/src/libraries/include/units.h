#if !defined(UNITS_H)
#define UNITS_H
// -*- C++ -*-
//
// Standard unit definitions.
// 
// Units can be made clear by writing, 1*k_cm or 1*k_m for example.

#include <math.h>

const double k_in_to_cm = 2.54;

const double k_m  = 100;
const double k_cm = 1.0;  // lengths are in centimeters
const double k_mm = 0.1;
const double k_um = 1.0E-4; // microns
const double k_in = k_in_to_cm * k_cm;
const double k_mil = k_in * 1.0E-3;
const double k_to_m = 1/k_m;
const double k_to_cm = 1/k_cm;
const double k_to_mm = 1/k_mm;
const double k_to_um = 1/k_um;
const double k_to_in = 1/k_in;
const double k_to_mil = 1/k_mil;

const double k_kg = 1000;
const double k_g  = 1.0;
const double k_u  = 1.6605402E-24;
const double k_to_kg = 1/k_kg;
const double k_to_g = 1/k_g;
const double k_to_u = 1/k_u;

const double k_kg_per_m3 = 0.001;
const double k_g_per_cm3 = 1;  // density is in g/cm^3
const double k_to_kg_per_m3 = 1/k_kg_per_m3;
const double k_to_g_per_cm3 = 1/k_g_per_cm3;

const double k_Pa  = 1; // pressure is in pascals
const double k_atm = 101325;
const double k_to_Pa = 1/k_Pa;
const double k_to_atm = 1/k_atm;

const double k_GeV = 1; // energy is in GeV
const double k_MeV = 1.0E-3;
const double k_keV = 1.0E-6;
const double k_eV = 1.0E-9;
const double k_to_GeV = 1/k_GeV;
const double k_to_MeV = 1/k_MeV;
const double k_to_keV = 1/k_keV;
const double k_to_eV = 1/k_eV;
const double k_to_J = 1.60217733E-10;
const double k_J = 1/k_to_J;

const double k_kGauss = 0.1;
const double k_Gauss = 1.0E-4;
const double k_Tesla = 1; // magnetic field in Tesla
const double k_to_kGauss = 1/k_kGauss;
const double k_to_Gauss = 1/k_Gauss;
const double k_to_Tesla = 1/k_Tesla;

const double k_radians = 1;
const double k_degrees = M_PI / 180;
const double k_to_radians = 1/k_radians;
const double k_to_degrees = 1/k_degrees;

const double k_psec = 1.0E-3;
const double k_nsec = 1; // time in nanoseconds
const double k_usec = 1.0E3;
const double k_msec = 1.0E6;
const double k_sec  = 1.0E9;
const double k_to_psec = 1/k_psec;
const double k_to_nsec = 1/k_nsec;
const double k_to_usec = 1/k_usec;
const double k_to_msec = 1/k_msec;
const double k_to_sec  = 1/k_sec;

#endif /* UNITS_H */
