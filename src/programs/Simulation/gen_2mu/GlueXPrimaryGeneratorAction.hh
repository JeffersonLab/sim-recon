//
// WARNING!!!!
//
// This file was taken from the HDGeant4 project and has been
// modified for use in this program. It may be out of sync with
// that and has definitely had some of its original funtionality
// disabled. You can compare this to the HDGeant4 version by looking
// for the HDGeant4 project on github.com.
//

//
// GlueXPrimaryGeneratorAction class header
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state.
// Separate object instances are created for each worker thread,
// but virtually all of its functions need to be serialized, so
// it maintains its own interlocks for this purpose. Resources
// are created once when the first object is instantiated, and
// destroyed once when the last object is destroyed.

#ifndef _GLUEXPRIMARYGENERATORACTION_H_
#define _GLUEXPRIMARYGENERATORACTION_H_

//#include "G4Threading.hh"
//#include "G4AutoLock.hh"

#include "CobremsGenerator.hh"
//#include "G4VUserPrimaryGeneratorAction.hh"
//#include "G4ParticleDefinition.hh"
//#include "GlueXParticleGun.hh"
//#include "G4SystemOfUnits.hh"
//
//#include "globals.hh"
//
//#include <HDDM/hddm_s.hpp>
//
//#include <fstream>
//
//class G4Event;

#include <TVector3.h>

//class GlueXPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
class GlueXPrimaryGeneratorAction
{
 public:
   
//   enum source_type_t {
//      SOURCE_TYPE_NONE,
//      SOURCE_TYPE_PARTICLE_GUN,
//      SOURCE_TYPE_COBREMS_GEN,
//      SOURCE_TYPE_HDDM
//   };

   GlueXPrimaryGeneratorAction();
//   GlueXPrimaryGeneratorAction(const GlueXPrimaryGeneratorAction &src);
//   GlueXPrimaryGeneratorAction &operator=(const GlueXPrimaryGeneratorAction &src);
   ~GlueXPrimaryGeneratorAction();
   
//   virtual void GeneratePrimaries(G4Event* anEvent);
//   void GeneratePrimariesHDDM(G4Event* anEvent);
//   void GeneratePrimariesParticleGun(G4Event* anEvent);
//   void GeneratePrimariesCobrems(G4Event* anEvent);
//   void GenerateBeamPhoton(G4Event* anEvent, double t0);
   void GenerateBeamPhoton(TVector3 &pgamma, TVector3 &pol);

//   int ConvertGeant3ToPdg(int Geant3number) const;
 
 private:
//   static int instanceCount;
//   static source_type_t fSourceType;
//   static std::ifstream *fHDDMinfile;
//   static hddm_s::istream *fHDDMistream;
   static CobremsGenerator *fCobremsGenerator;
//   static G4ParticleTable *fParticleTable;
//   static GlueXParticleGun *fParticleGun;

 public:
//   struct single_particle_gun_t {
//      int geantType;
//      int pdgType;
//      G4ParticleDefinition *partDef;
//      G4ThreeVector pos;
//      double mom;
//      double theta;
//      double phi;
//      double deltaR;
//      double deltaZ;
//      double deltaMom;
//      double deltaTheta;
//      double deltaPhi;
//   };

 private:
//   static single_particle_gun_t fGunParticle;

   static double fBeamBucketPeriod;
   static double fBeamBackgroundRate;
   static double fBeamBackgroundGateStart;
   static double fBeamBackgroundGateStop;
   static double fL1triggerTimeSigma;
   static double fBeamStartZ;

//   static int fEventCount;

   // The following parameters describe the dimensions of the target
   // that are used when generating the primary interaction vertex for
   // events from an external generator. An external event generator
   // knows nothing about the simulation geometry, so it makes sense
   // that this should be modeled in the simulation. They only apply
   // to the HDDM input source.  They are initialized to default values
   // for the GlueX liquid hydrogen target in the constructor, but
   // can be accessed/changed by the getter/setter methods below.
   static double fTargetCenterZ;
//   static double fTargetLength;
//   static double fBeamDiameter;

 public:
//   void setTargetCenterZ(double Z_cm) {
//      fTargetCenterZ = Z_cm * cm;
//   }
//   void setTargetLength(double L_cm) {
//      fTargetLength = L_cm * cm;
//   }
//   void setBeamDiameter(double D_cm) {
//      fBeamDiameter = D_cm * cm;
//   }
//   double getTargetCenterZ() {
//      return fTargetCenterZ / cm;
//   }
//   double getTargetLength() {
//      return fTargetLength / cm;
//   }
//   double getBeamDiameter() {
//      return fBeamDiameter / cm;
//   }
//
//   double getBeamBucketPeriod(int runno=0);
//
//   int getEventCount() {
//      return fEventCount;
//   }
//
//   void setBeamBucketPeriod(double period_ns) {
//      fBeamBucketPeriod = period_ns * ns;
//   }
//   void setL1triggerTimeSigma(double sigma_ns) {
//      fL1triggerTimeSigma = sigma_ns;
//   }
//   void setBeamStartZ(double Z_cm) {
//      fBeamStartZ = Z_cm * cm;
//   }
//   double getL1triggerTimeSigma() {
//      return fL1triggerTimeSigma;
//   }
//   double getBeamStartZcm() {
//      return fBeamStartZ / cm;
//   }

   // The following tables contain PDFs for importance-sampling the
   // kinematic variables in coherent bremsstrahlung beam generation.
 
   struct ImportanceSampler {
      std::vector<double> randvar;
      std::vector<double> density;
      std::vector<double> integral;
      double Psum;
      double Pcut;
      double Pmax;
      int Nfailed;
      int Npassed;

      ImportanceSampler()
       : Psum(1.0), Pcut(1), Pmax(0), Nfailed(0), Npassed(1) {}
//       : Psum(0), Pcut(1), Pmax(0), Nfailed(0), Npassed(0) {}
   };

   static ImportanceSampler fCoherentPDFx; 
   static ImportanceSampler fIncoherentPDFlogx;
   static ImportanceSampler fIncoherentPDFy;

 private:

   static double fIncoherentPDFtheta02;

   void prepareCobremsImportanceSamplingPDFs();

 private:
//   static G4Mutex fMutex;
};

#endif
