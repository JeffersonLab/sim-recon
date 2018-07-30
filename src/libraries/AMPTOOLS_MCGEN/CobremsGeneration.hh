//
// CobremsGeneration class header
//
// author: richard.t.jones at uconn.edu
// version: july 27, 2015
//
// notes:
//
// This class computes differential rates and polarization factors
// for coherent bremsstrahlung by an electron beam passing through
// a crystal radiator. A beamline geometry similar to that in Hall D
// at Jefferson Lab is assumed, consisting of a single radiator
// followed by a collimator located some distance away. Rates are
// computed for both the pre-collimated and post-collimated beams.
//
// This code was ported from cobrems.f, written in Fortran 77.
//
// units:
// Any length is in m; energy,momentum,mass in GeV (c=1); angles in
// radians; time in seconds; current in microAmps.

#ifndef CobremsGeneration_h
#define CobremsGeneration_h 1

#include <string>
#include <vector>

#if BOOST_PYTHON_WRAPPING
#include <boost/python.hpp>
#endif

class CobremsGeneration {
 public:
   CobremsGeneration(double Emax_GeV, double Epeak_GeV);
   CobremsGeneration(const CobremsGeneration &src);
   CobremsGeneration &operator=(const CobremsGeneration &src);
   ~CobremsGeneration();

   void setBeamEnergy(double Ebeam_GeV);
   void setBeamErms(double Erms_GeV);
   void setBeamEmittance(double emit_m_r);
   void setCollimatorSpotrms(double spotrms_m);
   void setCollimatorDistance(double distance_m);
   void setCollimatorDiameter(double diameter_m);
   void setTargetThickness(double thickness_m);
   void setTargetCrystal(std::string crystal);
   void setCoherentEdge(double Epeak_GeV);
   void setTargetThetax(double thetax);
   void setTargetThetay(double thetay);
   void setTargetThetaz(double thetaz);
   void setTargetOrientation(double thetax, double thetay, double thetaz);
   void setPhotonEnergyMin(double Emin_GeV);
   void RotateTarget(double thetax, double thetay, double thetaz);
   void setCollimatedFlag(bool flag);
   void setPolarizedFlag(bool flag);

   double getBeamEnergy() {
      return fBeamEnergy; // (GeV)
   }
   double getBeamErms() {
      return fBeamErms; // (GeV)
   }
   double getBeamEmittance() {
      return fBeamEmittance; // (m rad)
   }
   double getCollimatorSpotrms() {
      return fCollimatorSpotrms; // (m)
   }
   double getCollimatorDistance() {
      return fCollimatorDistance; // (m)
   }
   double getCollimatorDiameter() {
      return fCollimatorDiameter; // (m)
   }
   double getTargetThickness() {
      return fTargetThickness; // (m)
   }
   std::string getTargetCrystal() {
      return fTargetCrystal.name;
   }
   int getTargetCrystalNsites() {
      return fTargetCrystal.nsites;
   }
   double getTargetCrystalAtomicNumber() {
      return fTargetCrystal.Z;
   }
   double getTargetCrystalAtomicWeight() {
      return fTargetCrystal.A;  // (amu)
   }
   double getTargetCrystalDensity() {
      return fTargetCrystal.density; // (g/cm^3)
   }
   double getTargetCrystalLatticeConstant() {
      return fTargetCrystal.lattice_constant; // (m)
   }
   double getTargetCrystalRadiationLength() {
      return fTargetCrystal.radiation_length; // (m)
   }
   double getTargetCrystalDebyeWallerConst() {
      return fTargetCrystal.Debye_Waller_const; // (1/GeV^2)
   }
   double getTargetCrystalMosaicSpread() {
      return fTargetCrystal.mosaic_spread; // (rad)
   }
   double getTargetCrystalBetaFF() {
      return fTargetCrystal.betaFF; // (1/GeV^2)
   }
   double getTargetThetax() {
      return fTargetThetax; // (rad)
   }
   double getTargetThetay() {
      return fTargetThetay; // (rad)
   }
   double getTargetThetaz() {
      return fTargetThetaz; // (rad)
   }
   double getPhotonEnergyMin() {
      return fPhotonEnergyMin; // (GeV)
   }
   bool getCollimatedFlag() {
      return fCollimatedFlag;
   }
   bool getPolarizedFlag() {
      return fPolarizedFlag;
   }

   double getTargetRadiationLength_PDG();
   double getTargetRadiationLength_Schiff();
   double getTargetDebyeWallerConstant(double DebyeT_K, double T_K);
   void applyBeamCrystalConvolution(int nbins, double *xvalues, 
                                               double *yvalues);
#if BOOST_PYTHON_WRAPPING
   typedef boost::python::object pyobject;
   void pyApplyBeamCrystalConvolution(int nbins, pyobject xarr, pyobject yarr);
#endif
   void printBeamlineInfo();
   void printTargetCrystalInfo();
   double CoherentEnhancement(double x);
   double Rate_dNtdx(double x);
   double Rate_dNtdx(double x, double distance_m, double diameter_m);
   double Rate_dNtdk(double k_GeV);
   double Rate_dNcdx(double x);
   double Rate_dNcdx(double x, double distance_m, double diameter_m);
   double Rate_dNcdxdp(double x, double phi);
   double Rate_dNidx(double x);
   double Rate_dNBidx(double x);
   double Rate_dNidxdt2(double x, double theta2);
   double Rate_para(double x, double theta2, double phi);
   double Rate_ortho(double x, double theta2, double phi);
   double Polarization(double x, double theta2);
   double Polarization(double x, double theta2, double phi);
   double AbremsPolarization(double x, double theta2, double phi);
   double Acceptance(double theta2, double phi, 
                     double xshift_m, double yshift_m);
   double Acceptance(double theta2);
   double Sigma2MS(double thickness_m);
   double Sigma2MS_Kaune(double thickness_m);
   double Sigma2MS_PDG(double thickness_m);
   double Sigma2MS_Geant(double thickness_m);
   double Sigma2MS_Hanson(double thickness_m);

   // some math and physical constants
   static const double dpi;
   static const double me;
   static const double alpha;
   static const double hbarc;

   // statistical record from last sum over reciprocal lattice
   std::vector<double> fQ2theta2;
   std::vector<double> fQ2weight;

 private:
   void resetTargetOrientation();
   void updateTargetOrientation();

   // description of the radiator crystal lattice, here configured for diamond
   // but may be customized to describe any regular crystal
   struct lattice_vector {
      double x;
      double y;
      double z;
      lattice_vector()
       : x(0), y(0), z(0) {}
      lattice_vector(double ux, double uy, double uz)
       : x(ux), y(uy), z(uz) {}
      lattice_vector(const lattice_vector &src)
       : x(src.x), y(src.y), z(src.z) {}
      lattice_vector &operator=(const lattice_vector &src) {
         x = src.x;
         y = src.y;
         z = src.z;
         return *this;
      }
   };
   struct crystal_parameters_t {
      std::string name;
      int nsites;
      double Z;
      double A;                        // amu
      double density;                  // g/cm^3
      double lattice_constant;         // m
      double radiation_length;         // m
      double Debye_Waller_const;       // 1/GeV^2
      double mosaic_spread;            // rms radians
      double betaFF;                   // 1/GeV^2
      std::vector<lattice_vector> ucell_site;
      lattice_vector primaryHKL;
   } fTargetCrystal;
   double fTargetThickness;

   // orientation of the radiator with respect to the beam axis
   double fTargetThetax;		// the "small" angle
   double fTargetThetay;		// the "large" angle
   double fTargetThetaz;
   double fTargetRmatrix[3][3];

   // description of the beam at the radiator
   double fBeamEnergy;                 // GeV
   double fBeamErms;                   // GeV
   double fBeamEmittance;              // m radians
   double fCollimatorSpotrms;          // m
   double fCollimatorDistance;         // m
   double fCollimatorDiameter;         // m

   // flags to select kind of flux to be computed
   bool fCollimatedFlag;
   bool fPolarizedFlag;

   // parameters controlling Monte Carlo generation of photons
   double fPhotonEnergyMin;            // GeV
};

inline void CobremsGeneration::setBeamEmittance(double emit_m_r) {
   fBeamEmittance = emit_m_r;
}

inline void CobremsGeneration::setBeamEnergy(double Ebeam_GeV) {
   fBeamEnergy = Ebeam_GeV;
}

inline void CobremsGeneration::setBeamErms(double Erms_GeV) {
   fBeamErms = Erms_GeV;
}

inline void CobremsGeneration::setCollimatorSpotrms(double spotrms_m) {
   fCollimatorSpotrms = spotrms_m;
}

inline void CobremsGeneration::setCollimatorDistance(double distance_m) {
   fCollimatorDistance = distance_m;
}

inline void CobremsGeneration::setCollimatorDiameter(double diameter_m) {
   fCollimatorDiameter = diameter_m;
}

inline void CobremsGeneration::setTargetThickness(double thickness_m) {
   fTargetThickness = thickness_m;
}

inline void CobremsGeneration::setTargetThetax(double thetax) {
   fTargetThetax = thetax;
   updateTargetOrientation();
}

inline void CobremsGeneration::setTargetThetay(double thetay) {
   fTargetThetay = thetay;
   updateTargetOrientation();
}

inline void CobremsGeneration::setTargetThetaz(double thetaz) {
   fTargetThetaz = thetaz;
   updateTargetOrientation();
}

inline void CobremsGeneration::setTargetOrientation(double thetax,
                                                   double thetay,
                                                   double thetaz) {
   fTargetThetax = thetax;
   fTargetThetay = thetay;
   fTargetThetaz = thetaz;
   updateTargetOrientation();
}

inline void CobremsGeneration::setPhotonEnergyMin(double Emin_GeV) {
   fPhotonEnergyMin = Emin_GeV;
}

inline void CobremsGeneration::setCollimatedFlag(bool flag) {
   fCollimatedFlag = flag;
}

inline void CobremsGeneration::setPolarizedFlag(bool flag) {
   fPolarizedFlag = flag;
}

inline void CobremsGeneration::resetTargetOrientation() {
   fTargetRmatrix[0][0] = 1;
   fTargetRmatrix[0][1] = 0;
   fTargetRmatrix[0][2] = 0;
   fTargetRmatrix[1][0] = 0;
   fTargetRmatrix[1][1] = 1;
   fTargetRmatrix[1][2] = 0;
   fTargetRmatrix[2][0] = 0;
   fTargetRmatrix[2][1] = 0;
   fTargetRmatrix[2][2] = 1;
}

#endif
