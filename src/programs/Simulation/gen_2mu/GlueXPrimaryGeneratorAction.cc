//
// class implementation for GlueXPrimaryGeneratorAction
//
// author: richard.t.jones at uconn.edu
// version: may 12, 2012

#include <stdlib.h>

#include <iostream>
#include <map>
#include <cmath>
using namespace std;

#include <TRandom2.h>
#include <TVector3.h>

#define s      (1.0)
#define ns     (1.0E-9)
#define m      (1.0)
#define cm     (1.0E-2)
#define GeV    (1.0)
#define radian (1.0)

#define twopi 6.28318530717958623
#define pi    3.14159265358979312
#define electron_mass_c2 0.000511

#ifndef _DBG_
#define _DBG_ cout<<__FILE__<<":"<<__LINE__<<" "
#define _DBG__ cout<<__FILE__<<":"<<__LINE__<<endl
#endif

static TRandom2 RAND(1);

#include "GlueXPrimaryGeneratorAction.hh"
//#include "GlueXUserEventInformation.hh"
//#include "GlueXUserOptions.hh"
//
//#include "G4Event.hh"
//#include "G4ParticleGun.hh"
//#include "G4ParticleTable.hh"
//#include "G4ParticleDefinition.hh"
//#include "G4SystemOfUnits.hh"
//#include "G4PhysicalConstants.hh"
//#include "Randomize.hh"
//
//#include <JANA/jerror.h>
//#include <JANA/JApplication.h>
//#include <JANA/JCalibration.h>
//
//typedef GlueXPrimaryGeneratorAction::source_type_t source_type_t;
//typedef GlueXPrimaryGeneratorAction::single_particle_gun_t particle_gun_t;
//typedef GlueXPrimaryGeneratorAction::ImportanceSampler ImportanceSampler;
//
//int GlueXPrimaryGeneratorAction::instanceCount = 0;
//source_type_t GlueXPrimaryGeneratorAction::fSourceType = SOURCE_TYPE_NONE;
//
//std::ifstream *GlueXPrimaryGeneratorAction::fHDDMinfile;
//hddm_s::istream *GlueXPrimaryGeneratorAction::fHDDMistream;
CobremsGenerator *GlueXPrimaryGeneratorAction::fCobremsGenerator;
//G4ParticleTable *GlueXPrimaryGeneratorAction::fParticleTable;
//GlueXParticleGun *GlueXPrimaryGeneratorAction::fParticleGun;
//particle_gun_t GlueXPrimaryGeneratorAction::fGunParticle;

//double GlueXPrimaryGeneratorAction::fBeamBucketPeriod = 0;
//double GlueXPrimaryGeneratorAction::fBeamBackgroundRate = 0;
//double GlueXPrimaryGeneratorAction::fBeamBackgroundGateStart = 0;
//double GlueXPrimaryGeneratorAction::fBeamBackgroundGateStop = 0;
//double GlueXPrimaryGeneratorAction::fL1triggerTimeSigma = 10 * ns;
//double GlueXPrimaryGeneratorAction::fBeamStartZ = -24 * m;
double GlueXPrimaryGeneratorAction::fTargetCenterZ = 1 * cm;
//double GlueXPrimaryGeneratorAction::fTargetLength = 29.9746 * cm;
//double GlueXPrimaryGeneratorAction::fBeamDiameter = 0.5 * cm;

//int GlueXPrimaryGeneratorAction::fEventCount = 0;

GlueXPrimaryGeneratorAction::ImportanceSampler GlueXPrimaryGeneratorAction::fCoherentPDFx; 
GlueXPrimaryGeneratorAction::ImportanceSampler GlueXPrimaryGeneratorAction::fIncoherentPDFlogx;
GlueXPrimaryGeneratorAction::ImportanceSampler GlueXPrimaryGeneratorAction::fIncoherentPDFy;
double GlueXPrimaryGeneratorAction::fIncoherentPDFtheta02;

//G4Mutex GlueXPrimaryGeneratorAction::fMutex = G4MUTEX_INITIALIZER;

//--------------------------------------------
// GlueXPrimaryGeneratorAction (constructor)
//--------------------------------------------

GlueXPrimaryGeneratorAction::GlueXPrimaryGeneratorAction()
{
//   G4AutoLock barrier(&fMutex);
//   ++instanceCount;
//
//   // Initializaton is driven by the control.in file, which
//   // gets read and parsed only once, by the first constructor.
//
//   if (fSourceType != SOURCE_TYPE_NONE) {
//      return;
//   }
//
//   fParticleGun = new GlueXParticleGun();
//   fParticleTable = G4ParticleTable::GetParticleTable();
//   
//   GlueXUserOptions *user_opts = GlueXUserOptions::GetInstance();
//   if (user_opts == 0) {
//      cerr << "Error in GlueXPrimaryGeneratorAction constructor - "
//             << "GlueXUserOptions::GetInstance() returns null, "
//             << "cannot continue." << endl;
//      exit(-1);
//   }
//
//   fHDDMinfile = 0;
//   fHDDMistream = 0;
//   fCobremsGenerator = 0;
//   std::map<int,std::string> infile;
   std::map<int,double> beampars;
   std::map<int,double> kinepars;

//   // Three event source options are supported:
//   // 1) external generator, hddm input stream source
//   // 2) internal coherent bremsstrahlung beam generator
//   // 3) internal particle gun generator
// 
//   if (user_opts->Find("INFILE", infile) ||
//       user_opts->Find("INFI", infile))
//   {
//      fHDDMinfile = new std::ifstream(infile[1].c_str());
//      if (!fHDDMinfile->is_open()) {
//         cerr << "GlueXPrimaryGeneratorAction error: "
//                << "Unable to open HDDM input file: " << infile[1]
//                << endl;
//         exit(-1);
//      }
//      fHDDMistream = new hddm_s::istream(*fHDDMinfile);
//      cout << "Opened input file: " << infile[1] << endl;
//      fSourceType = SOURCE_TYPE_HDDM;
//   }
//
//   else if (user_opts->Find("BEAM", beampars))
//   {
      double beamE0 = beampars[1];
      double beamEpeak = beampars[2];
      double beamEmin = (beampars[3] > 0)? beampars[3] : 0.120;
      double radColDist = (beampars[4] > 0)? beampars[4] : 76.;
      double colDiam = (beampars[5] > 0)? beampars[5] : 0.0034;
      double beamEmit = (beampars[6] > 0)? beampars[6] : 2.5e-9;
      double radThick = (beampars[7] > 0)? beampars[7] : 20e-6;

beamE0=12.0;
beamEpeak = 9.0;
beamEmin = 1.0;

      if (beamE0 == 0 || beamEpeak == 0) {
         cerr << "GlueXPrimaryGeneratorAction error: "
                << "BEAM card specified in control.in but required values "
                << "Emax and/or Epeak are missing, cannot continue."
                << endl;
         exit(-1);
      }

      fCobremsGenerator = new CobremsGenerator(beamE0, beamEpeak);
      fCobremsGenerator->setPhotonEnergyMin(beamEmin);
      fCobremsGenerator->setCollimatorDistance(radColDist);
      fCobremsGenerator->setCollimatorDiameter(colDiam);
      fCobremsGenerator->setBeamEmittance(beamEmit);
      fCobremsGenerator->setTargetThickness(radThick);
      prepareCobremsImportanceSamplingPDFs();

//      std::map<int, double> bgratepars;
//      std::map<int, double> bggatepars;
//      if (user_opts->Find("BGRATE", bgratepars) &&
//          user_opts->Find("BGGATE", bggatepars))
//      {
//         fBeamBackgroundRate = bgratepars[1] * 1/s;
//         fBeamBackgroundGateStart = bgratepars[1] * ns;
//         fBeamBackgroundGateStop = bgratepars[2] * ns;
//         if (fBeamBackgroundRate > 0 &&
//             fBeamBackgroundGateStart >= fBeamBackgroundGateStop)
//         {
//            cerr << "GlueXPrimaryGeneratorAction error: "
//                   << "BGRATE is non-zero, but the time window specified "
//                   << "in BGGATE is invalid."
//                   << endl;
//            exit(-1);
//         }
//      }
//      fSourceType = SOURCE_TYPE_COBREMS_GEN;
//   }
//
//   else if (user_opts->Find("KINE", kinepars))
//   {
//      if (kinepars[1] == 1000) {
//         fGunParticle.geantType = 0;
//         fGunParticle.pdgType = 999999;
//         fGunParticle.partDef = fParticleTable->FindParticle("geantino");
//      }
//      else if (kinepars[1] == 1001) {
//         fGunParticle.geantType = 0;
//         fGunParticle.pdgType = 999999;
//         fGunParticle.partDef = fParticleTable->FindParticle("chargedgeantino");
//      }
//      else {
//         if (kinepars[1] > 100)
//            fGunParticle.geantType = kinepars[1] - 100;
//         else
//            fGunParticle.geantType = kinepars[1];
//         fGunParticle.pdgType = ConvertGeant3ToPdg(fGunParticle.geantType);
//         fGunParticle.partDef = fParticleTable->FindParticle(fGunParticle.pdgType);
//      }
//      if (fGunParticle.partDef == 0) {   
//         cerr << "GlueXPrimaryGeneratorAction constructor error - "
//                << "Unknown GEANT particle type: " << kinepars[1] 
//                << " was specified in the control.in file." << endl;
//         exit(-1);
//      }
//      fParticleGun->SetParticleDefinition(fGunParticle.partDef);
//
//      double x(0), y(0), z(65 * cm);
//      std::map<int,double> scappars;
//      if (user_opts->Find("SCAP", scappars)) {
//         x = scappars[1] * cm;
//         y = scappars[2] * cm;
//         z = scappars[3] * cm;
//      }
//      fGunParticle.pos.set(x,y,z);
//      std::map<int,double> tgtwidthpars;
//      if (user_opts->Find("tgtwidth", tgtwidthpars)) {
//         fGunParticle.deltaR = tgtwidthpars[1] * cm;
//         fGunParticle.deltaZ = tgtwidthpars[2] * cm;
//      }
//      else {
//         fGunParticle.deltaR = 0;
//         fGunParticle.deltaZ = 0;
//      }
//
//      fGunParticle.mom = kinepars[2] * GeV;
//      if (kinepars[1] > 100) {
//         fGunParticle.theta = kinepars[3] * degree;
//         fGunParticle.phi = kinepars[4] * degree;
//         fGunParticle.deltaMom = kinepars[5];
//         fGunParticle.deltaTheta = kinepars[6];
//         fGunParticle.deltaPhi = kinepars[7];
//      }
//      else {
//         fGunParticle.deltaMom = 0;
//         fGunParticle.theta = 90 * degree;
//         fGunParticle.deltaTheta = 180 * degree;
//         fGunParticle.phi = 0;
//         fGunParticle.deltaPhi = 360 * degree;
//      }
//      fSourceType = SOURCE_TYPE_PARTICLE_GUN;
//   }
//
//   std::map<int,double> trefsigma;
//   if (user_opts->Find("trefsigma", trefsigma)) {
//      fL1triggerTimeSigma = trefsigma[1] * ns;
//   }
//   else {
//      fL1triggerTimeSigma = 10 * ns;
//   }
}

//GlueXPrimaryGeneratorAction::GlueXPrimaryGeneratorAction(const
//                             GlueXPrimaryGeneratorAction &src)
// : G4VUserPrimaryGeneratorAction(src)
//{
//   G4AutoLock barrier(&fMutex);
//   ++instanceCount;
//}
//
//GlueXPrimaryGeneratorAction &GlueXPrimaryGeneratorAction::operator=(const
//                             GlueXPrimaryGeneratorAction &src)
//{
//   *(G4VUserPrimaryGeneratorAction*)this = src;
//   return *this;
//}

//--------------------------------------------
// ~GlueXPrimaryGeneratorAction (destructor)
//--------------------------------------------

GlueXPrimaryGeneratorAction::~GlueXPrimaryGeneratorAction()
{
//   G4AutoLock barrier(&fMutex);
//   if (--instanceCount == 0) {
//      if (fHDDMistream)
//         delete fHDDMistream;
//      if (fHDDMinfile)
//         delete fHDDMinfile;
      if (fCobremsGenerator)
         delete fCobremsGenerator;
//      delete fParticleGun;
//   }
}

void GlueXPrimaryGeneratorAction::prepareCobremsImportanceSamplingPDFs()
{

	cout << "Preparing Cobrems Importance Sampling PDFs ..." << endl;

   // Construct lookup tables representing the PDFs used for
   // importance-sampling the coherent bremsstrahlung kinematics.

   const int Ndim = 500;
   double Emin = fCobremsGenerator->getPhotonEnergyMin() * GeV;
   double Emax = fCobremsGenerator->getBeamEnergy() * GeV;
   double sum;

   // Compute approximate PDF for dNc/dx
   double xmin = Emin / Emax;
   double dx = (1 - xmin) / Ndim;
   double xarr[Ndim + 1], yarr[Ndim + 1];
   for (int i=0; i <= Ndim; ++i) {
      xarr[i] = xmin + i * dx;
      yarr[i] = twopi * fCobremsGenerator->Rate_dNcdxdp(xarr[i], pi/2);
   }
   fCobremsGenerator->applyBeamCrystalConvolution(Ndim + 1, xarr, yarr);
   sum = 0;
   for (int i=0; i <= Ndim; ++i) {
      sum += (i > 0)? (yarr[i] + yarr[i - 1]) / 2 : 0;
      fCoherentPDFx.randvar.push_back(xarr[i]);
      fCoherentPDFx.density.push_back(yarr[i]);
      fCoherentPDFx.integral.push_back(sum);
   }
   for (int i=0; i <= Ndim; ++i) {
      fCoherentPDFx.density[i] /= sum * dx;
      fCoherentPDFx.integral[i] /= sum;
   }

   // Compute approximate PDF for dNi/dx
   double logxmin = log(xmin);
   double dlogx = -logxmin / Ndim;
   sum = 0;
   for (int i=0; i <= Ndim; ++i) {
      double logx = logxmin + i * dlogx;
      double x = exp(logx);
      double dNidx = fCobremsGenerator->Rate_dNidxdt2(x, 0);
      double dNidlogx = dNidx * x;
      fIncoherentPDFlogx.randvar.push_back(logx);
      fIncoherentPDFlogx.density.push_back(dNidlogx);
      fIncoherentPDFlogx.integral.push_back(sum);
      sum += (i < Ndim)? dNidlogx : 0;
   }
   for (int i=0; i <= Ndim; ++i) {
      fIncoherentPDFlogx.density[i] /= sum * dlogx;
      fIncoherentPDFlogx.integral[i] /= sum;
   }
 
   // Compute approximate PDF for dNi/dy
   fIncoherentPDFtheta02 = 1.8;
   double ymin = 1e-3;
   double dy = (1 - ymin) / Ndim;
   sum = 0;
   for (int i=0; i <= Ndim; ++i) {
      double y = ymin + i * dy;
      double theta2 = fIncoherentPDFtheta02 * (1 / y - 1);
      double dNidxdt2 = fCobremsGenerator->Rate_dNidxdt2(0.5, theta2);
      fIncoherentPDFy.randvar.push_back(y);
      fIncoherentPDFy.density.push_back(dNidxdt2);
      fIncoherentPDFy.integral.push_back(sum);
      sum += (i < Ndim)? dNidxdt2 : 0;
   }
   for (int i=0; i <= Ndim; ++i) {
      fIncoherentPDFy.density[i] /= sum * dy;
      fIncoherentPDFy.integral[i] /= sum;
   }

   // These cutoffs should be set empirically, as low as possible
   // for good efficiency, but not too low so as to avoid excessive
   // warnings about Pcut violations.
   fCoherentPDFx.Pcut = .0001;
   fIncoherentPDFlogx.Pcut = .02;

   cout << "Completed Cobrems Importance Sampling PDFs." << endl;
}
//
////--------------------------------------------
//// GeneratePrimaries
////--------------------------------------------
//
//void GlueXPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
//{
//   G4AutoLock barrier(&fMutex);
//
//   switch(fSourceType){
//      case SOURCE_TYPE_HDDM:
//         GeneratePrimariesHDDM(anEvent);
//         break;
//      case SOURCE_TYPE_COBREMS_GEN:
//         GeneratePrimariesCobrems(anEvent);
//         break;
//      case SOURCE_TYPE_PARTICLE_GUN:
//         GeneratePrimariesParticleGun(anEvent);
//         break;
//      default:
//         cout << "No event source selected, cannot continue!" << endl;
//         exit(-1);
//   }
//}   
//
////--------------------------------------------
//// GeneratePrimariesParticleGun
////--------------------------------------------
//
//void GlueXPrimaryGeneratorAction::GeneratePrimariesParticleGun(G4Event* anEvent)
//{
//   // Unbelievably, GEANT4's G4ParticleGun class insists on printing
//   // a message whenever the momentum or energy is changed, unless
//   // the other is 0. Here, we reset the particle gun energy using 
//   // our own derived class. (Sheesh!!)
//   fParticleGun->Reset();
//
//   // place and smear the particle gun origin
//   G4ThreeVector pos(fGunParticle.pos);
//   if (fGunParticle.deltaR > 0) {
//      double dx, dy;
//      while (true) {
//         double rx = RAND.Rndm() - 0.5;
//         double ry = RAND.Rndm() - 0.5;
//         if (rx*rx + ry*ry <= 0.25) {
//            dx = rx * 2 * fGunParticle.deltaR;
//            dy = ry * 2 * fGunParticle.deltaR;
//            break;
//         }
//      }
//      pos += G4ThreeVector(dx, dy, 0);
//   }
//   if (fGunParticle.deltaZ > 0) {
//      double dz = (RAND.Rndm() - 0.5) * fGunParticle.deltaZ;
//      pos += G4ThreeVector(0, 0, dz);
//   }
//   fParticleGun->SetParticlePosition(pos);
//
//   // Assign and optionally smear the particle momentum
//   double p = fGunParticle.mom;
//   double thetap = fGunParticle.theta;
//   double phip = fGunParticle.phi;
//   if (fGunParticle.deltaMom > 0)
//      p += (RAND.Rndm() - 0.5) * fGunParticle.deltaMom;
//   if (fGunParticle.deltaTheta > 0)
//      thetap += (RAND.Rndm() - 0.5) * fGunParticle.deltaTheta;
//   if (fGunParticle.deltaPhi > 0)
//      phip += (RAND.Rndm() - 0.5) * fGunParticle.deltaPhi;
//   G4ThreeVector mom(p * sin(thetap) * cos(phip),
//                     p * sin(thetap) * sin(phip),
//                     p * cos(thetap));
//   fParticleGun->SetParticleMomentum(mom);
//
//   // Set the event number and fire the gun
//   anEvent->SetEventID(++fEventCount);
//   fParticleGun->GeneratePrimaryVertex(anEvent);
//
//   // Store generated particle info so it can be written to output file
//   pos *= 1 / cm; // convert to cm
//   mom *= 1 / GeV; // convert to GeV
//   int type = fGunParticle.geantType;
//   anEvent->SetUserInformation(new GlueXUserEventInformation(type, pos, mom));
//}
//
////--------------------------------------------
//// GeneratePrimariesHDDM
////--------------------------------------------
//
//void GlueXPrimaryGeneratorAction::GeneratePrimariesHDDM(G4Event* anEvent)
//{
//   if (! fHDDMinfile->good()) {
//      anEvent->SetEventAborted();
//      return;
//   }
//
//   hddm_s::HDDM *hddmevent = new hddm_s::HDDM;
//   try {
//      *fHDDMistream >> *hddmevent;
//   }
//   catch(std::exception e) {
//      cout << e.what() << endl;
//      anEvent->SetEventAborted();
//      return;
//   }
//
//   // Store generated event info so it can be written to output file
//   ++fEventCount;
//   anEvent->SetUserInformation(new GlueXUserEventInformation(hddmevent));
//
//   // Unpack generated event and prepare initial state for simulation
//   int Nprimaries = 0;
//   hddm_s::VertexList vertices = hddmevent->getVertices();
//   if (vertices.size() == 0) {
//      cout << "No vertices in input HDDM event!" << endl;
//      anEvent->SetEventAborted();
//      return;
//   }
//   hddm_s::VertexList::iterator it_vertex;
//   for (it_vertex = vertices.begin();
//        it_vertex != vertices.end(); ++it_vertex)
//   {
//      anEvent->SetEventID(it_vertex->getEventNo());
//      hddm_s::Origin &origin = it_vertex->getOrigin();
//      double x = origin.getVx() * cm;
//      double y = origin.getVy() * cm;
//      double z = origin.getVz() * cm;
//      if (x == 0 && y == 0 && z == 0) {
//         while (true) {
//            x = RAND.Rndm() - 0.5;
//            y = RAND.Rndm() - 0.5;
//            if (x*x + y*y <= 0.25) {
//               x *= fBeamDiameter;
//               y *= fBeamDiameter;
//            }
//         }
//         z = fTargetCenterZ + (RAND.Rndm() - 0.5) * fTargetLength;
//      }
//      G4ThreeVector pos(x, y, z);
//
//      // The primary interaction vertex time is referenced to a clock
//      // whose t=0 is synchronized to the crossing of a beam bunch
//      // through the target midplane. This beam bunch may not contain
//      // the beam particle whose interaction generated the vertex,
//      // but it represents best-guess based on the arrival time of
//      // the L1 trigger signal. The spread in the L1 relative to the
//      // interacting bunch time is parameterized as a Gaussian.
//
//      extern int run_number;
//      if (fBeamBucketPeriod == 0)
//         getBeamBucketPeriod(run_number);
//         // getBeamBucketPeriod(it_vertex->getRunNo());
//
//      double t0, t0rf;
//      double lightSpeed = 2.99792e8 * m/s;
//      t0 = (origin.getT() * ns) + fL1triggerTimeSigma * RAND.Gaus::shoot();
//      t0rf = fBeamBucketPeriod * int(t0 / fBeamBucketPeriod + 0.5);
//      t0 = t0rf + (z - fTargetCenterZ) / lightSpeed;
//      G4PrimaryVertex* vertex = new G4PrimaryVertex(pos, t0);
//
//      hddm_s::ProductList &products = it_vertex->getProducts();
//      hddm_s::ProductList::iterator it_product;
//      for (it_product = products.begin();
//           it_product != products.end(); ++it_product)
//      {
//         // ignore intermediaries in the MC record
//         if (it_product->getType() <= 0)
//           continue;
//
//         int g3type = it_product->getType();
//         int pdgtype = it_product->getPdgtype();
//         G4ParticleDefinition *part;
//         if (pdgtype > 0 && pdgtype < 999999) {
//            part = fParticleTable->FindParticle(pdgtype);
//         }
//         else if (g3type > 0) {
//            pdgtype = ConvertGeant3ToPdg(g3type);
//            part = fParticleTable->FindParticle(pdgtype);
//#if FORCE_PARTICLE_TYPE_CHARGED_GEANTINO
//            part = fParticleTable->FindParticle("chargedgeantino");
//#endif
//         }
//         else {
//            cerr << "Unknown particle found in input MC record, "
//                   << "geant3 type " << g3type 
//                   << ", PDG type " << pdgtype
//                   << ", failing over to geantino!"
//                   << endl;
//            part = fParticleTable->FindParticle("geantino");
//         }
//         hddm_s::Momentum &momentum = it_product->getMomentum();
//         double px = momentum.getPx() * GeV;
//         double py = momentum.getPy() * GeV;
//         double pz = momentum.getPz() * GeV;
//         double Etot = momentum.getE() * GeV;
//         vertex->SetPrimary(new G4PrimaryParticle(part, px, py, pz, Etot));
//         ++Nprimaries;
//      }
//      anEvent->AddPrimaryVertex(vertex);
//   }
//   
//   if (Nprimaries == 0) {
//      cerr << "Number of primaries in event is zero!!" << endl;
//      anEvent->SetEventAborted();
//   }
//
//   // Superimpose any request background minimum-bias beam interactions
//
//   if (fBeamBackgroundRate > 0) {
//      double t = fBeamBackgroundGateStart;
//      while (true) {
//         t += -log(RAND.Rndm()) / fBeamBackgroundRate;
//         if (t > fBeamBackgroundGateStop)
//            break;
//         GenerateBeamPhoton(anEvent, t);
//      }
//   }
//}
//
//void GlueXPrimaryGeneratorAction::GeneratePrimariesCobrems(G4Event* anEvent)
//{
//   GenerateBeamPhoton(anEvent, 0);
//   ++fEventCount;
//}

void GlueXPrimaryGeneratorAction::GenerateBeamPhoton(TVector3 &pgamma, TVector3 &pol)
{
   // Generates a single beam photon according to the coherent bremsstrahlung
   // model defined by class CobremsGenerator.  The photon begins its lifetime
   // just upstream of the primary collimator (WARNING: position is hard-wired
   // in the code below) and is tracked by the simulation from there forward.
   // Its time t0 should identify its beam bucket, ie. the time the photon
   // would reach the midplane of the target. To enable beam motion spreading,
   // define the beam box size below.

#define BEAM_PHOTON_START_Z (-24 * m)
// #define BEAM_BOX_SIZE (5 * mm)

   // The algorithm below generates coherent bremsstrahlung photons using a
   // importance-sampling technique. This algorithm requires that we prepare
   // an approximate probability density function for the generated photons.
   // The function is not in general equal to the true physical PDF, which
   // varies from event to event depending on the direction of the incident
   // electron, and also the exact angle of the crystal planes at the point
   // of scattering which moves because of the mosaic spread of the crystal.
   // The important thing is that the approximate PDF be reasonably close to
   // the average over all beam particles and the crystal mosaic, and that
   // deviations from event to event are sufficiently small that rejection
   // sampling can be used to take them into account with high efficiency.
   //
   // The kinematics of bremsstrahlung are described by three independent
   // variables (x, theta, phi) where x is the photon energy in units of the
   // incident electron energy, and theta,phi are the polar,azimuthal angles
   // of the photon in a lab frame tilted so that the incident electron comes
   // in along the z axis. Polar angle theta is represented by dimensionless
   // variable y = theta0^2 / (theta^2 + theta0^2) where contant theta0 is
   // chosen to optimize the uniformity of the PDF in y. On each event,
   // a new random tuple (x, phi, y) is generated on the interval x:[0,1],
   // phi:[0,2pi], y:[0,1] using a split-and-recombine strategy. One side 
   // of the split covers the coherent process and the other side covers the
   // incoherent one.
   //
   //  1) coherent process - the PDF here is continuous in x,phi according
   //     the dNc/(dx dphi), and the dependence on y is a sequence of delta 
   //     functions corresponding to the different planes that contribute to
   //     the scattering at the given value of x. Here we take advantage of
   //     the fact that the marginal distribution dNc/dx is proportional to
   //     dNc/(dx dphi) at phi=pi/4. This allows us to decompose the generation
   //     into two stages, first generating x from dNc/dx and then generating
   //     phi from dNc/(dx dphi) at fixed x. The x generation step is performed
   //     using importance sampling based on the average PDF stored in table
   //     fCoherentPDF, followed by rejection sampling based on the value of
   //     dNc/(dx dphi) computed for the particular kinematics of each event.
   //     The y value is obtained by sampling the weighted list of q2 values
   //     that contributed the to q-sum in the calculation of dNc/(dx dphi).
   //
   //  2) incoherent process - the PDF here is continuous in x,phi,y but it
   //     is uniform in phi, so it is effectively a 2D distribution. Here we
   //     take advantage of the fact that x and y are independent variables
   //     to a good approximation, which allows us to generate x using
   //     importance sampling from an approximation to dNi/(dx dtheta^2) at
   //     theta=0 and y ~ uniform [0,1], then employ rejection sampling based
   //     on the exact PDF dNi/(dx dtheta2) to get a true sample.
   //
   // Recombination after the split is very simple. First we integrate over
   // phi in both cases to obtain values dNc/dx and dNi/(dx dy). It turns
   // out that in spite of the fact that the y-dependence is discrete in the
   // coherent case and continuous in the incoherent case, the sum over the
   // probabilities for all values of y in dNc/dx is always normalized to 1
   // independently for all values of x. Hence we can treat y as a psuedo
   // coordinate y' ~ Unif[0,1] and form a 2D PDF dNc/(dx dy') which is
   // numerically equal to dNc/dx, do the rejection sampling in combination
   // with that applied to dNi/(dx dy) and then replace the fake variable y'
   // with the true y that was sampled as described above.

   double phiMosaic = twopi * RAND.Rndm();
   double rhoMosaic = sqrt(-2 * log(RAND.Rndm()));
   rhoMosaic *= fCobremsGenerator->getTargetCrystalMosaicSpread() * m*radian;
   double thxMosaic = rhoMosaic * cos(phiMosaic);
   double thyMosaic = rhoMosaic * sin(phiMosaic);

   double xemittance = fCobremsGenerator->getBeamEmittance() * m;
   double yemittance = xemittance / 2.5; // nominal, should be checked
   double xspotsize = fCobremsGenerator->getCollimatorSpotrms() * m;
   double yspotsize = xspotsize; // nominal, should be checked
   double thxBeam = (xemittance / xspotsize) * sqrt(-2 * log(RAND.Rndm()));
   double thyBeam = (yemittance / yspotsize) * sqrt(-2 * log(RAND.Rndm()));

   double raddz = fCobremsGenerator->getTargetThickness() * m;
   double varMS = fCobremsGenerator->Sigma2MS(raddz * RAND.Rndm());
   double thxMS = sqrt(-2 * varMS * log(RAND.Rndm()));
   double thyMS = sqrt(-2 * varMS * log(RAND.Rndm()));

   double targetThetax = fCobremsGenerator->getTargetThetax() * radian;
   double targetThetay = fCobremsGenerator->getTargetThetay() * radian;
   double targetThetaz = fCobremsGenerator->getTargetThetaz() * radian;
   double thetax = thxBeam + thxMS - targetThetax - thxMosaic;
   double thetay = thyBeam + thyMS - targetThetay - thyMosaic;
   double thetaz = -targetThetaz;
   fCobremsGenerator->resetTargetOrientation();
   fCobremsGenerator->RotateTarget(0, pi/2, 0);    // point (1,0,0) along beam
   fCobremsGenerator->RotateTarget(0, 0, pi/4);    // point (0,1,1) vertically
   fCobremsGenerator->RotateTarget(thetax, thetay, thetaz);

   // Generate with importance sampling
   double x, phi, theta2=0;
   double polarization = 0;
   double Scoherent = fCoherentPDFx.Npassed / 
                      (fCoherentPDFx.Psum / fCoherentPDFx.Npassed);
   double Sincoherent = fIncoherentPDFlogx.Npassed /
                        (fIncoherentPDFlogx.Psum / fIncoherentPDFlogx.Npassed);
   if (Scoherent < Sincoherent) {
      while (true) {                             // try coherent generation
         double dNcdxPDF=0.0;
         double u = RAND.Rndm();
         for (unsigned int i=1; i < fCoherentPDFx.randvar.size(); ++i) {
            if (u <= fCoherentPDFx.integral[i]) {
               double x0 = fCoherentPDFx.randvar[i - 1];
               double x1 = fCoherentPDFx.randvar[i];
               double f0 = fCoherentPDFx.density[i - 1];
               double f1 = fCoherentPDFx.density[i];
               double u0 = fCoherentPDFx.integral[i - 1];
               double u1 = fCoherentPDFx.integral[i];
               x = (x0 * (u1 - u) + x1 * (u - u0)) / (u1 - u0);
               dNcdxPDF = (f0 * (u1 - u) + f1 * (u - u0)) / (u1 - u0);
               break;
            }
         }
         double dNcdx = twopi * fCobremsGenerator->Rate_dNcdxdp(x, pi / 4);
         double Pfactor = dNcdx / dNcdxPDF;
         if (Pfactor > fCoherentPDFx.Pmax)
            fCoherentPDFx.Pmax = Pfactor;
         if (Pfactor > fCoherentPDFx.Pcut) {
            cout << "Warning in GenerateBeamPhoton - Pfactor " << Pfactor
                   << " exceeds fCoherentPDFx.Pcut = " << fCoherentPDFx.Pcut
                   << ", please increase." << endl;
         }
         if (RAND.Rndm() * fCoherentPDFx.Pcut > Pfactor) {
            ++fCoherentPDFx.Nfailed;
            continue;
         }
         fCoherentPDFx.Psum += Pfactor;
         ++fCoherentPDFx.Npassed;

         double fmax = dNcdx / pi;
         while (true) {
            phi = twopi * RAND.Rndm();
            double f = fCobremsGenerator->Rate_dNcdxdp(x, phi);
            if (RAND.Rndm() * fmax < f)
               break;
         }
         double uq = RAND.Rndm();
         for (unsigned int i=0; i < fCobremsGenerator->fQ2theta2.size(); ++i) {
            if (uq <= fCobremsGenerator->fQ2weight[i]) {
               theta2 = fCobremsGenerator->fQ2theta2[i];
               break;
            }
         }
         polarization = fCobremsGenerator->Polarization(x, theta2);
		 break;
      }
   }
   else {
      while (true) {                           // try incoherent generation
         double dNidxdyPDF;
         double u = RAND.Rndm();
         for (unsigned int i=1; i < fIncoherentPDFlogx.randvar.size(); ++i) {
            if (u <= fIncoherentPDFlogx.integral[i]) {
               double logx0 = fIncoherentPDFlogx.randvar[i - 1];
               double logx1 = fIncoherentPDFlogx.randvar[i];
               double f0 = fIncoherentPDFlogx.density[i - 1];
               double f1 = fIncoherentPDFlogx.density[i];
               double u0 = fIncoherentPDFlogx.integral[i - 1];
               double u1 = fIncoherentPDFlogx.integral[i];
               double logx = (logx0 * (u1 - u) + logx1 * (u - u0)) / (u1 - u0);
               dNidxdyPDF = (f0 * (u1 - u) + f1 * (u - u0)) / (u1 - u0);
               x = exp(logx);
               break;
            }
         }
         double y=0.0;
         double uy = RAND.Rndm();
         for (unsigned int i=1; i < fIncoherentPDFy.randvar.size(); ++i) {
            if (uy <= fIncoherentPDFy.integral[i]) {
               double y0 = fIncoherentPDFy.randvar[i - 1];
               double y1 = fIncoherentPDFy.randvar[i];
               double f0 = fIncoherentPDFy.density[i - 1];
               double f1 = fIncoherentPDFy.density[i];
               double u0 = fIncoherentPDFy.integral[i - 1];
               double u1 = fIncoherentPDFy.integral[i];
               y = (y0 * (u1 - uy) + y1 * (uy - u0)) / (u1 - u0);
               dNidxdyPDF *= (f0 * (u1 - uy) + f1 * (uy - u0)) / (u1 - u0);
               break;
            }
         }
         theta2 = fIncoherentPDFtheta02 * (1 / (y + 1e-99) - 1);
         double dNidxdy = fCobremsGenerator->Rate_dNidxdt2(x, theta2) *
                          fIncoherentPDFtheta02 / (y*y + 1e-99);
         double Pfactor = dNidxdy / dNidxdyPDF;
         if (Pfactor > fIncoherentPDFlogx.Pmax)
            fIncoherentPDFlogx.Pmax = Pfactor;
         if (Pfactor > fIncoherentPDFlogx.Pcut) {
            cout << "Warning in GenerateBeamPhoton - Pfactor " << Pfactor
                   << " exceeds fIncoherentPDFlogx.Pcut = " 
                   << fIncoherentPDFlogx.Pcut << ", please increase."
                   << endl;
         }
         if (RAND.Rndm() * fIncoherentPDFlogx.Pcut > Pfactor) {
            ++fIncoherentPDFlogx.Nfailed;
            continue;
         }
         ++fIncoherentPDFlogx.Psum += Pfactor;
         ++fIncoherentPDFlogx.Npassed;

         phi = twopi * RAND.Rndm();
         polarization = 0;
		 break;
      }
   }

   // Define the particle kinematics and polarization in lab coordinates
//   G4ParticleDefinition *part = fParticleTable->FindParticle("gamma");
   double Emax = fCobremsGenerator->getBeamEnergy() * GeV;
   double Erms = fCobremsGenerator->getBeamErms() * GeV;
   double Ebeam = Emax + Erms * RAND.Gaus();
   double theta = sqrt(theta2) * electron_mass_c2 / Emax;
   double alphax = thxBeam + thxMS + theta * cos(phi);
   double alphay = thyBeam + thyMS + theta * sin(phi);
   double pabs = Ebeam * x;
   double px = pabs * alphax;
   double py = pabs * alphay;
   double pz = sqrt(pabs*pabs - px*px - py*py);
//   double colphi = twopi * RAND.Rndm();
//   double vspotrms = fCobremsGenerator->getCollimatorSpotrms() * m;
//   double colrho = vspotrms * sqrt(-2 * log(RAND.Rndm()));
//   double colDist = fCobremsGenerator->getCollimatorDistance() * m;
//   double radx = colrho * cos(colphi) - colDist * thxBeam;
//   double rady = colrho * sin(colphi) - colDist * thyBeam;
//   double colx = radx + colDist * alphax;
//   double coly = rady + colDist * alphay;
//#if defined BEAM_BOX_SIZE
//   colx += BEAM_BOX_SIZE * (RAND.Rndm() - 0.5);
//   coly += BEAM_BOX_SIZE * (RAND.Rndm() - 0.5);
//#endif
//   TVector3 vtx(colx, coly, fBeamStartZ);
//   TVector3 pol(0, polarization, -polarization * py / pz);
	pgamma.SetXYZ(px, py, pz);
	pol.SetXYZ(0, polarization, -polarization * py / pz);

//   // Generate a new primary for the beam photon
//   double lightSpeed = 2.99792e8 * m/s;
//   double t0rf = fBeamBucketPeriod * int(t0 / fBeamBucketPeriod + 0.5);
//   t0rf += (fBeamStartZ - fTargetCenterZ) / lightSpeed;
//   G4PrimaryVertex* vertex = new G4PrimaryVertex(vtx, t0rf);
//   G4PrimaryParticle* photon = new G4PrimaryParticle(part, px, py, pz);
//   photon->SetPolarization(pol);
//   vertex->SetPrimary(photon);
//   anEvent->AddPrimaryVertex(vertex);
      
   // call hitTagger(vertex,vertex,plab,plab,0.,1,0,0)
}
//
//// Convert particle types from Geant3 types to PDG scheme
//
//int GlueXPrimaryGeneratorAction::ConvertGeant3ToPdg(int Geant3number) const
//{
//   // This method was imported from ROOT source file TDatabasePDG.cc
//
//   switch(Geant3number) {
//
//      case 1   : return 22;       // photon
//      case 25  : return -2112;    // anti-neutron
//      case 2   : return -11;      // e+
//      case 26  : return -3122;    // anti-Lambda
//      case 3   : return 11;       // e-
//      case 27  : return -3222;    // Sigma-
//      case 4   : return 12;       // e-neutrino (NB: flavour undefined by Geant)
//      case 28  : return -3212;    // Sigma0
//      case 5   : return -13;      // mu+
//      case 29  : return -3112;    // Sigma+ (PB)*/
//      case 6   : return 13;       // mu-
//      case 30  : return -3322;    // Xi0
//      case 7   : return 111;      // pi0
//      case 31  : return -3312;    // Xi+
//      case 8   : return 211;      // pi+
//      case 32  : return -3334;    // Omega+ (PB)
//      case 9   : return -211;     // pi-
//      case 33  : return -15;      // tau+
//      case 10  : return 130;      // K long
//      case 34  : return 15;       // tau-
//      case 11  : return 321;      // K+
//      case 35  : return 411;      // D+
//      case 12  : return -321;     // K-
//      case 36  : return -411;     // D-
//      case 13  : return 2112;     // n
//      case 37  : return 421;      // D0
//      case 14  : return 2212;     // p
//      case 38  : return -421;     // D0
//      case 15  : return -2212;    // anti-proton
//      case 39  : return 431;      // Ds+
//      case 16  : return 310;      // K short
//      case 40  : return -431;     // anti Ds-
//      case 17  : return 221;      // eta
//      case 41  : return 4122;     // Lamba_c+
//      case 18  : return 3122;     // Lambda
//      case 42  : return 24;       // W+
//      case 19  : return 3222;     // Sigma+
//      case 43  : return -24;      // W-
//      case 20  : return 3212;     // Sigma0
//      case 44  : return 23;       // Z
//      case 21  : return 3112;     // Sigma-
//      case 45  : return 0;        // deuteron
//      case 22  : return 3322;     // Xi0
//      case 46  : return 0;        // triton
//      case 23  : return 3312;     // Xi-
//      case 47  : return 0;        // alpha
//      case 24  : return 3334;     // Omega- (PB)
//      case 48  : return 0;        // G nu ? PDG ID 0 is undefined
//
//      default  : return 0;
//
//   }
//}
//
//double GlueXPrimaryGeneratorAction::getBeamBucketPeriod(int runno)
//{
//   // Look up the beam bucket period for this run in ccdb
//   // unless the user has already set the value by hand.
//
//   if (runno > 0) {
//      jana::JCalibration *jcalib = japp->GetJCalibration(runno);
//      std::map<std::string, double> result;
//      std::string map_key("/PHOTON_BEAM/RF/rf_period");
//      if (jcalib->Get(map_key, result)) {
//         cerr << "Error in GeneratePrimariesHDDM - "
//                << "error fetching " << map_key << " from ccdb, "
//                << "cannot continue." << endl;
//         exit(-1);
//      }
//      else if (result.find("rf_period") != result.end()) {
//         fBeamBucketPeriod = result["rf_period"] * ns;
//      }
//      else {
//         cerr << "Error in GeneratePrimariesHDDM - "
//                << "error finding value for " << map_key
//                << " in ccdb, cannot continue." << endl;
//         exit(-1);
//      }
//   }
//   return fBeamBucketPeriod;
//}
