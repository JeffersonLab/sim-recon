// $Id$
//
// Created Oct 10, 2013  David Lawrence

#include "hddm_select_events.h"
#include "particleType.h"

using namespace std;

extern TRandom2 *rndm;

#define SQR(X) ((X)*(X))

//-----------------------------------------------------------------
// selectEvent_s
//-----------------------------------------------------------------
bool selectEvent_s(int select_type, hddm_s::HDDM &record, 
                   int nevents, bool debug)
{
  if (record.getPhysicsEvents().size() == 0) {
    std::cout << "PhysicsEvent not given" << std::endl;
    return -1;
  }

  // Select Lambda -> p pi- events
  if (select_type == 1) {
    bool foundLambda = false;
    int  idLambda    = -999;
    int  idProton    = -999;
    bool foundProton = false;
    if (debug)
    std::cout << "---------------------------------"
              << "---------------------------------"
              << std::endl;

    hddm_s::VertexList vs = record.getVertices();
    hddm_s::VertexList::iterator viter;
    for (viter = vs.begin(); viter != vs.end(); ++viter) {
       hddm_s::OriginList orig = viter->getOrigins();
       hddm_s::ProductList ps = viter->getProducts();
       if (ps.size() > 0 && orig.size() > 0) {
          hddm_s::ProductList::iterator piter;
          for (piter = ps.begin(); piter != ps.end(); ++piter) {
             int type     = (int)piter->getType();
             int id       = piter->getId();
             int parentid = piter->getParentid();
             double E  = piter->getMomentum().getE();
             double px = piter->getMomentum().getPx();
             double py = piter->getMomentum().getPy();
             double pz = piter->getMomentum().getPz();
             double mass = sqrt(E*E - (px*px + py*py + pz*pz));
             if (! finite(mass))
                mass = 0.0;
             if (debug) {
                std::cout << "\t\t\t--- new particle ---" << std::endl;
                std::cout << "\t\t\ttype = " << type << " " 
                          << ParticleType((Particle_t)type) 
                          << " parentid = " << parentid << std::endl;
                std::cout << "\t\t\t(" << px << "," << py << "," << pz << ","
                          << E << "), mass = " << mass << std::endl;
                std::cout << "\t\t\t(" << orig().getVx() << "," 
                          << orig().getVy() << "," << orig().getVz() << ","
                          << orig().getT() << ")" << std::endl;
             }

             // Find Lambda
             if (type == Lambda) {
                foundLambda = true;
                idLambda = id;
             }

             // Find decay products of Lambda
             if (foundLambda) {
                // if parent was Lambda
                if (parentid == idLambda) {
                   // find what type it is
                   if (type == Proton) {
                      foundProton = true;
                      idProton = id;
                   }
                }
             }

             // If we found a proton decaying from Lambda,
             // save this event
             if (foundProton) {
                return true;
             }
          } // end of loop over products
       } // end of having products and origin
    } // end of loop over vertices
  } // end of select_type 1
  //___________________________________________________________________________

  // select Lambda -> p pi-, eta -> gamma gamma events
  else if (select_type == 2) {
    bool foundLambda = false;
    int  idLambda    = -999;
    int  idProton    = -999;
    bool foundProton = false;
    bool foundEta = false;
    int  idEta    = -999;
    int  idPhoton[2];
    int  nPhoton = 0;
    if (debug)
       std::cout << "------------------------------"
                    "------------------------------------"
                 << std::endl;
      
    hddm_s::VertexList vs = record.getVertices();
    hddm_s::VertexList::iterator viter;
    for (viter = vs.begin(); viter != vs.end(); ++viter) {
       hddm_s::OriginList orig = viter->getOrigins();
       hddm_s::ProductList ps = viter->getProducts();
       if (ps.size() > 0 && orig.size() > 0) {
          hddm_s::ProductList::iterator piter;
          for (piter = ps.begin(); piter != ps.end(); ++piter) {
             int type     = (int)piter->getType();
             int id       = piter->getId();
             int parentid = piter->getParentid();
             double E  = piter->getMomentum().getE();
             double px = piter->getMomentum().getPx();
             double py = piter->getMomentum().getPy();
             double pz = piter->getMomentum().getPz();
             double mass = sqrt(E*E - (px*px + py*py + pz*pz));
             if (! finite(mass))
                mass = 0.0;
             if (debug) {
                std::cout << "\t\t\t--- new particle ---" << std::endl;
                std::cout << "\t\t\ttype = " << type << " " 
                          << ParticleType((Particle_t)type) 
                          << " parentid = " << parentid << std::endl;
                std::cout << "\t\t\t(" << px << "," << py << "," << pz << ","
                          << E << "), mass = " << mass << std::endl;
                std::cout << "\t\t\t(" << orig().getVx() << "," 
                          << orig().getVy() << "," << orig().getVz() << ","
                          << orig().getT() << ")" << std::endl;
             }

             // Find Lambda
             if (type == Lambda) {
               foundLambda = true;
               idLambda = id;
             }

             // Find decay products of Lambda
             if (foundLambda) {
                // if parent was Lambda
                if (parentid == idLambda) {
                   // find what type it is
                   if (type == Proton) {
                      foundProton = true;
                      idProton = id;
                   }
                }
             }

             // Find eta
             if (type == Eta) {
                foundEta = true;
                idEta = id;
             }

             // Find decay products of Eta
             if (foundEta) {
                // if parent was Eta
                if (parentid == idEta) {
                   // find what type it is
                   if (type == Gamma) {
                      idPhoton[nPhoton] = id;
                      nPhoton++;
                   }
                }
             }
          }
       }

       // If we found a proton decaying from Lambda,
       // and we found 2 photons from eta decay,
       // save this event
       if (foundProton && nPhoton == 2) {
          return true;
       }
    } // end of loop over vertices
  }
  //___________________________________________________________________________

  // Select Lambda -> p pi- events, decay length is less than VERTEXDIFFMAX
  else if (select_type == 3) {
    const double VERTEXDIFFMAX = 2.0;
    bool foundLambda = false;
    int  idLambda    = -999;
    int  idProton    = -999;
    bool foundProton = false;
    double vertexLambda[3];
    double vertexProton[3];
    for (int v=0; v<3; v++) {
       vertexLambda[v] = -999;
       vertexProton[v] = 999;
    }
    double vertexdiff = 999;
    if (debug)
       std::cout << "----------------------------------"
                    "--------------------------------" << std::endl;
 
    hddm_s::VertexList vs = record.getVertices();
    hddm_s::VertexList::iterator viter;
    for (viter = vs.begin(); viter != vs.end(); ++viter) {
       hddm_s::OriginList orig = viter->getOrigins();
       hddm_s::ProductList ps = viter->getProducts();
       if (ps.size() > 0 && orig.size() > 0) {
          hddm_s::ProductList::iterator piter;
          for (piter = ps.begin(); piter != ps.end(); ++piter) {
             int type     = (int)piter->getType();
             int id       = piter->getId();
             int parentid = piter->getParentid();
             double E  = piter->getMomentum().getE();
             double px = piter->getMomentum().getPx();
             double py = piter->getMomentum().getPy();
             double pz = piter->getMomentum().getPz();
             double mass = sqrt(E*E - (px*px + py*py + pz*pz));
             if (! finite(mass))
                mass = 0.0;
             if (debug) {
                std::cout << "\t\t\t--- new particle ---" << std::endl;
                std::cout << "\t\t\ttype = " << type << " " 
                          << ParticleType((Particle_t)type) 
                          << " parentid = " << parentid << std::endl;
                std::cout << "\t\t\t(" << px << "," << py << "," << pz << ","
                          << E << "), mass = " << mass << std::endl;
                std::cout << "\t\t\t(" << orig().getVx() << "," 
                          << orig().getVy() << "," << orig().getVz() << ","
                          << orig().getT() << ")" << std::endl;
             }

             // Find Lambda
             if (type == Lambda) {
                foundLambda = true;
                idLambda = id;
                vertexLambda[0] = orig().getVx();
                vertexLambda[1] = orig().getVy();
                vertexLambda[2] = orig().getVz();
             }

             // Find decay products of Lambda
             if (foundLambda) {
                // if parent was Lambda
                if (parentid == idLambda) {
                  // find what type it is
                  if (type == Proton) {
                    foundProton = true;
                    idProton = id;
                    vertexProton[0] = orig().getVx();
                    vertexProton[1] = orig().getVy();
                    vertexProton[2] = orig().getVz();
                    vertexdiff = sqrt(SQR(vertexLambda[0] - vertexProton[0]) +
                                      SQR(vertexLambda[1] - vertexProton[1]) +
                                      SQR(vertexLambda[2] - vertexProton[2]));
                  }
                }
             }

             // If we found a proton decaying from Lambda,
             // and the vertex difference is less than VERTEXDIFFMAX,
             // save this event
             if (foundProton && vertexdiff < VERTEXDIFFMAX) {
               return true;
             }

          } // end of loop over products
       } // end of having products and origin
    } // end of loop over vertices
  } // end of select_type 3
  //___________________________________________________________________________

  // Select Lambda -> p pi- events, carve out weak decay asymmetry
  else if (select_type == 4) {
    bool foundKPlus = false;
    bool foundLambda = false;
    int  idLambda    = -999;
    int  idProton    = -999;
    bool foundProton = false;
    TLorentzVector p4photon_init;
    TLorentzVector p4proton_init(0,0,0,ParticleMass(Proton));
    TLorentzVector p4kp;
    TLorentzVector p4Lambda;
    TLorentzVector p4proton;
    TVector3 beamdir(0,0,1);
    TVector3 normaldir;
    double phi_n = 0;
    if (debug) 
       std::cout << "---------------------------------"
                    "---------------------------------" << std::endl;

    hddm_s::ReactionList res = record.getReactions();
    hddm_s::ReactionList::iterator riter;
    for (riter = res.begin(); riter != res.end(); ++riter) {
       // Get initial photon
       p4photon_init.SetXYZT(riter->getBeam().getMomentum().getPx(),
                             riter->getBeam().getMomentum().getPy(),
                             riter->getBeam().getMomentum().getPz(),
                             riter->getBeam().getMomentum().getE());
       if (debug)
          std::cout << "photon.        : (" << setw(5) << p4photon_init.X()
                    << ", " << setw(5) << p4photon_init.Y() << ", " 
                    << p4photon_init.Z() << ", " 
                    << setw(5) << p4photon_init.E() << ")" << std::endl;

       hddm_s::VertexList vs = res().getVertices();
       hddm_s::VertexList::iterator viter;
       for (viter = vs.begin(); viter != vs.end(); ++viter) {
          hddm_s::OriginList orig = viter->getOrigins();
          hddm_s::ProductList ps = viter->getProducts();
          if (ps.size() > 0 && orig.size() > 0) {
             hddm_s::ProductList::iterator piter;
             for (piter = ps.begin(); piter != ps.end(); ++piter) {
                int type     = (int)piter->getType();
                int id       = piter->getId();
                int parentid = piter->getParentid();
                double E  = piter->getMomentum().getE();
                double px = piter->getMomentum().getPx();
                double py = piter->getMomentum().getPy();
                double pz = piter->getMomentum().getPz();
                double mass = sqrt(E*E - (px*px + py*py + pz*pz));
                if (! finite(mass))
                   mass = 0.0;
                if (debug) {
                std::cout << "\t\t\t--- new particle ---" << std::endl;
                   std::cout << "\t\t\ttype = " << type << " " 
                             << ParticleType((Particle_t)type) 
                             << " parentid = " << parentid << std::endl;
                   std::cout << "\t\t\t(" << px << "," << py << "," << pz << ","
                             << E << "), mass = " << mass << std::endl;
                   std::cout << "\t\t\t(" << orig().getVx() << "," 
                             << orig().getVy() << "," << orig().getVz() << ","
                             << orig().getT() << ")" << std::endl;
                }

                // Find K+, set up polarization axis
                if (type == KPlus) {
                   foundKPlus = true;
                   p4kp.SetXYZT(px,py,pz,E);
                   normaldir = p4photon_init.Vect().Cross(p4kp.Vect()).Unit();
                   phi_n = normaldir.Phi();

                   if (debug) {
                      std::cout << "K+ dir.    : (" << setw(5) << p4kp.X() 
                                << ", " << setw(5) << p4kp.Y() << ", " 
                                << p4kp.Z() << ")" << std::endl;
                      std::cout << "normal dir.: (" << setw(5)
                                << normaldir.X() << ", " << setw(5) 
                                << normaldir.Y() << ", " 
                                << normaldir.Z() << ")" << std::endl;
                   }
                }

                // Find Lambda
                if (type == Lambda) {
                   foundLambda = true;
                   idLambda = id;
                   p4Lambda.SetXYZT(px,py,pz,E);
                }

                // Find decay products of Lambda
                if (foundLambda) {
                   // if parent was Lambda
                   if (parentid == idLambda) {
                      // find what type it is
                      if (type == Proton) {
                         foundProton = true;
                         idProton = id;
                         p4proton.SetXYZT(px,py,pz,E);
                         if (debug)
                            std::cout << "proton : (" 
                                      << setw(5) << p4proton.X() << ", "
                                      << setw(5) << p4proton.Y() << ", "
                                      << p4proton.Z() << ")" << std::endl;
                      }
                   }
                }
             } // end of loop over products
          } // end of having products and origin
       } // end of loop over vertices
    } // end of loop over reaction

    // If we found a proton decaying from Lambda,
    // save this event
    if (foundKPlus && foundProton) {
       // Get proton angular dist
       // First boost to c.m. frame
       p4Lambda.Boost(-(p4photon_init + p4proton_init).BoostVector());
       p4proton.Boost(-(p4photon_init + p4proton_init).BoostVector());

       // 2. Rotate around z by pi/2 - phi_n
       // This takes K+ into xz plane, and normal into y-dir.
       p4Lambda.RotateZ(TMath::Pi()/2. - phi_n);
       p4proton.RotateZ(TMath::Pi()/2. - phi_n);

       // 3. Rotate around x by pi/2
       // This takes K+ into xy plane, and normal into z-dir.
       p4Lambda.RotateX(TMath::Pi()/2);
       p4proton.RotateX(TMath::Pi()/2);

       // 4. Boost to Lambda rest frame.
       p4proton.Boost(-p4Lambda.BoostVector());
       // Now get decay angle
       double costheta_proton = p4proton.CosTheta();

       // Throw random number and select events
       double rand =rndm->Rndm();
       if (debug) {
          std::cout << "random number: " << rand << std::endl;
          std::cout << "compare to   : " 
                    << (1. + alpha * costheta_proton) / (1. + alpha) 
                    << std::endl;
       }
       if (rand < (1. + alpha * costheta_proton) / (1. + alpha)) {
          std::cout << " --- " << costheta_proton << std::endl;
          return true;
       }
    } // found K+ and proton
  } // end of select_type 4
  //___________________________________________________________________________

  // Select p pi+ pi- events
  else if (select_type == 5) {
    bool foundProton = false;
    bool foundPip    = false;
    bool foundPim    = false;
    int  idProton    = -999;
    int  idPip       = -999;
    int  idPim       = -999;
    TLorentzVector p4photon_init;
    TLorentzVector p4proton_init(0,0,0,ParticleMass(Proton));
    TLorentzVector p4proton;
    TLorentzVector p4pip;
    TLorentzVector p4pim;
    TLorentzVector p4diff;
    if (debug)
       std::cout << "----------------------------------"
                    "--------------------------------" << std::endl;

    hddm_s::ReactionList res = record.getReactions();
    hddm_s::ReactionList::iterator riter;
    for (riter = res.begin(); riter != res.end(); ++riter) {
       // Get initial photon
       p4photon_init.SetXYZT(riter->getBeam().getMomentum().getPx(),
                             riter->getBeam().getMomentum().getPy(),
                             riter->getBeam().getMomentum().getPz(),
                             riter->getBeam().getMomentum().getE());
       if (debug)
          std::cout << "photon.        : (" << setw(5) 
                    << p4photon_init.X() << ", " << setw(5) 
                    << p4photon_init.Y() << ", " << p4photon_init.Z() 
                    << ", " << setw(5) << p4photon_init.E() << ")"
                    << std::endl;

       hddm_s::VertexList vs = res().getVertices();
       hddm_s::VertexList::iterator viter;
       for (viter = vs.begin(); viter != vs.end(); ++viter) {
          hddm_s::OriginList orig = viter->getOrigins();
          hddm_s::ProductList ps = viter->getProducts();
          if (ps.size() > 0 && orig.size() > 0) {
             hddm_s::ProductList::iterator piter;
             for (piter = ps.begin(); piter != ps.end(); ++piter) {
                int type     = (int)piter->getType();
                int id       = piter->getId();
                int parentid = piter->getParentid();
                double E  = piter->getMomentum().getE();
                double px = piter->getMomentum().getPx();
                double py = piter->getMomentum().getPy();
                double pz = piter->getMomentum().getPz();
                double mass = sqrt(E*E - (px*px + py*py + pz*pz));
                if (! finite(mass))
                   mass = 0.0;
                if (debug) {
                std::cout << "\t\t\t--- new particle ---" << std::endl;
                   std::cout << "\t\t\ttype = " << type << " " 
                             << ParticleType((Particle_t)type) 
                             << " parentid = " << parentid << std::endl;
                   std::cout << "\t\t\t(" << px << "," << py << "," << pz << ","
                             << E << "), mass = " << mass << std::endl;
                   std::cout << "\t\t\t(" << orig().getVx() << "," 
                             << orig().getVy() << "," << orig().getVz() << ","
                             << orig().getT() << ")" << std::endl;
                }

                // Find proton
                if (type == Proton) {
                   foundProton = true;
                   idProton = id;
                   p4proton.SetXYZT(px,py,pz,E);
                }

                // Find pi+
                if (type == PiPlus) {
                   foundPip = true;
                   idPip = id;
                   p4pip.SetXYZT(px,py,pz,E);
                }

                // Find pi-
                if (type == PiMinus) {
                   foundPim = true;
                   idPim = id;
                   p4pim.SetXYZT(px,py,pz,E);
                }
             } // end of loop over products
          } // end of having products and origin
       } // end of loop over vertices
    } // end of loop over reaction

    // If we found a proton decaying from Lambda,
    // save this event
    if (foundProton && foundPip && foundPim) {

       // make sure that total 4-mom of p, pi+, pi- are close to
       // initial state
       p4diff = p4photon_init + p4proton_init - p4proton - p4pip - p4pim;

       // Make sure that remaining |p| < 1 MeV, M < 0.1 MeV
       if (p4diff.P() < 0.001 && p4diff.M() < 0.0001)
          return true;
    } // found proton, pi+, pi-

  } // end of select_type 5
  //___________________________________________________________________________

  // Select K+ K+ Xi- -> K+ K+ p pi- pi- events
  else if (select_type == 6) {
    Int_t nKp             = 0;
    Int_t nXiMinus        = 0;
    Int_t nLambda         = 0;
    Int_t nProton         = 0;
    Int_t nPimFromLambda  = 0;
    Int_t nPimFromXiMinus = 0;
    // need intermediate IDs
    int  idXiMinus        = -999;
    int  idLambda         = -999;
    TLorentzVector p4photon_init;
    TLorentzVector p4proton_init(0,0,0,ParticleMass(Proton));
    TLorentzVector p4kp1;
    TLorentzVector p4kp2;
    TLorentzVector p4XiMinus;
    TLorentzVector p4Lambda;
    TLorentzVector p4proton;
    TLorentzVector p4pimFromXiMinus;
    TLorentzVector p4pimFromLambda;
    TLorentzVector p4diff;
    TLorentzVector p4total;
    if (debug) 
       std::cout << "-------------------------------"
                    "-----------------------------------" << std::endl;

    hddm_s::ReactionList res = record.getReactions();
    hddm_s::ReactionList::iterator riter;
    for (riter = res.begin(); riter != res.end(); ++riter) {
       // Get initial photon
       p4photon_init.SetXYZT(riter->getBeam().getMomentum().getPx(),
                             riter->getBeam().getMomentum().getPy(),
                             riter->getBeam().getMomentum().getPz(),
                             riter->getBeam().getMomentum().getE());
       if (debug)
          std::cout << "photon.        : (" << setw(5) 
                    << p4photon_init.X() << ", " << setw(5) 
                    << p4photon_init.Y() << ", " << p4photon_init.Z() 
                    << ", " << setw(5) << p4photon_init.E() << ")"
                    << std::endl;

       hddm_s::VertexList vs = res().getVertices();
       hddm_s::VertexList::iterator viter;
       for (viter = vs.begin(); viter != vs.end(); ++viter) {
          hddm_s::OriginList orig = viter->getOrigins();
          hddm_s::ProductList ps = viter->getProducts();
          if (ps.size() > 0 && orig.size() > 0) {
             hddm_s::ProductList::iterator piter;
             for (piter = ps.begin(); piter != ps.end(); ++piter) {
                int type     = (int)piter->getType();
                int id       = piter->getId();
                int parentid = piter->getParentid();
                double E  = piter->getMomentum().getE();
                double px = piter->getMomentum().getPx();
                double py = piter->getMomentum().getPy();
                double pz = piter->getMomentum().getPz();
                double mass = sqrt(E*E - (px*px + py*py + pz*pz));
                if (! finite(mass))
                   mass = 0.0;
                if (debug) {
                std::cout << "\t\t\t--- new particle ---" << std::endl;
                   std::cout << "\t\t\ttype = " << type << " " 
                             << ParticleType((Particle_t)type) 
                             << " parentid = " << parentid << std::endl;
                   std::cout << "\t\t\t(" << px << "," << py << "," << pz << ","
                             << E << "), mass = " << mass << std::endl;
                   std::cout << "\t\t\t(" << orig().getVx() << "," 
                             << orig().getVy() << "," << orig().getVz() << ","
                             << orig().getT() << ")" << std::endl;
                }

                // Find first K+
                if (nKp == 0 && type == KPlus) {
                   nKp++;
                   p4kp1.SetXYZT(px,py,pz,E);
                   if (debug)
                      std::cout << "K+ 1  : (" << p4kp1.X() << ", " 
                                << p4kp1.Y() << ", " << p4kp1.Z() << ", "
                                << p4kp1.M() << ")" << std::endl;

                   // Need this break since without it we would count the
                   // the first K+ as the second K+ too
                   continue;
                }

                // Find second K+
                if (nKp == 1 && type == KPlus) {
                   nKp++;
                   p4kp2.SetXYZT(px,py,pz,E);
                   if (debug)
                      std::cout << "K+ 2  : (" << p4kp2.X() << ", " 
                                << p4kp2.Y() << ", " << p4kp2.Z() << ", "
                                << p4kp2.M() << ")" << std::endl;
                }

                // Find Xi-
                if (type == XiMinus) {
                   nXiMinus++;
                   idXiMinus = id;
                   p4XiMinus.SetXYZT(px,py,pz,E);
                }

                // Find pi- from Xi-
                if (type == PiMinus && parentid == idXiMinus) {
                   nPimFromXiMinus++;
                   p4pimFromXiMinus.SetXYZT(px,py,pz,E);
                }

                // Find Lambda from Xi-
                if (type == Lambda && parentid == idXiMinus) {
                   nLambda++;
                   idLambda = id;
                   p4Lambda.SetXYZT(px,py,pz,E);
                }

                // Find proton from Lambda
                if (type == Proton && parentid == idLambda) {
                   nProton++;
                   p4proton.SetXYZT(px,py,pz,E);
                }

                // Find pi- from Lambda
                if (type == PiMinus && parentid == idLambda) {
                   nPimFromLambda++;
                   p4pimFromLambda.SetXYZT(px,py,pz,E);
                }
             } // end of loop over products
          } // end of having products and origin
       } // end of loop over vertices
    } // end of loop over reaction

    // If we found a proton decaying from Lambda,
    // save this event
    if (nKp == 2 && nXiMinus == 1 && nLambda == 1 && nProton == 1 &&
        nPimFromXiMinus == 1 && nPimFromLambda == 1)
    {
       std::cout << "found all particles" << std::endl;

       // make sure that total 4-mom of p, pi+, pi- are close to
       // initial state
       p4diff  = p4photon_init + p4proton_init - p4kp1 - p4kp2 - 
                 p4proton - p4pimFromXiMinus - p4pimFromLambda;
       p4total = p4kp1 + p4kp2 + p4proton + p4pimFromXiMinus + p4pimFromLambda;

       if (debug) {
          std::cout << "photon: (" << p4photon_init.X() << ", " 
                    << p4photon_init.Y() << ", " << p4photon_init.Z() 
                    << ", " << p4photon_init.M() << ")" << std::endl;
          std::cout << "proton: (" << p4proton_init.X() << ", " 
                    << p4proton_init.Y() << ", " << p4proton_init.Z() 
                    << ", " << p4proton_init.M() << ")" << std::endl;
          std::cout << "Xi-   : (" << p4XiMinus.X() << ", " 
                    << p4XiMinus.Y() << ", " << p4XiMinus.Z() 
                    << ", " << p4XiMinus.M() << ")" << std::endl;
          std::cout << "Lambda: (" << p4Lambda.X() << ", " 
                    << p4Lambda.Y() << ", " << p4Lambda.Z() << ", " 
                    << p4Lambda.M() << ")" << std::endl;
          std::cout << "K+ 1  : (" << p4kp1.X() << ", " << p4kp1.Y() 
                    << ", " << p4kp1.Z() << ", " << p4kp1.M() << ")" 
                    << std::endl;
          std::cout << "K+ 2  : (" << p4kp2.X() << ", " << p4kp2.Y() 
                    << ", " << p4kp2.Z() << ", " << p4kp2.M() << ")" 
                    << std::endl;
          std::cout << "proton: (" << p4proton.X() << ", " << p4proton.Y() 
                    << ", " << p4proton.Z() << ", " << p4proton.M() << ")"
                    << std::endl;
          std::cout << "pi- 1 : (" << p4pimFromXiMinus.X() << ", " 
                    << p4pimFromXiMinus.Y() << ", " << p4pimFromXiMinus.Z()
                    << ", " << p4pimFromXiMinus.M() << ")" << std::endl;
          std::cout << "pi- 2 : (" << p4pimFromLambda.X() << ", " 
                    << p4pimFromLambda.Y() << ", " << p4pimFromLambda.Z() 
                    << ", " << p4pimFromLambda.M() << ")" << std::endl;
       }

       // Make sure that remaining |p| < 1 MeV, M < 0.1 MeV
       if (p4diff.P() < 0.001 && p4diff.M() < 0.0001) {
          std::cout << "returning true" << std::endl;
          return true;
       }
       else {
          std::cout << "Added 4-mom of final state particles does not match"
                       " generated by more than 0.1 MeV!!!!" << std::endl;
          std::cout << "diff in |p| = " << p4diff.P()* 1000. << " MeV/c" 
                    << std::endl;
          std::cout << "diff  = (" << p4diff.X() << ", " << p4diff.Y() 
                    << ", " << p4diff.Z() << ", " << p4diff.E() << ")" 
                    << std::endl;
          std::cout << "total = (" << p4total.X() << ", " << p4total.Y()
                    << ", " << p4total.Z() << ", " << p4total.E() << ")"
                    << std::endl;
          abort();
       }
    } // found proton, pi+, pi-
  } // end of select_type 6
  //___________________________________________________________________________

  // Select p pi+ pi- pi0 events
  else if (select_type == 7) {
    bool foundProton = false;
    bool foundPip    = false;
    bool foundPim    = false;
    bool foundPi0    = false;
    int  idProton    = -999;
    int  idPip       = -999;
    int  idPim       = -999;
    int  idPi0       = -999;
    TLorentzVector p4photon_init;
    TLorentzVector p4proton_init(0,0,0,ParticleMass(Proton));
    TLorentzVector p4proton;
    TLorentzVector p4pip;
    TLorentzVector p4pim;
    TLorentzVector p4pi0;
    TLorentzVector p4diff;
    if (debug) 
       std::cout << "---------------------------------"
                    "---------------------------------" << std::endl;

    hddm_s::ReactionList res = record.getReactions();
    hddm_s::ReactionList::iterator riter;
    for (riter = res.begin(); riter != res.end(); ++riter) {
       // Get initial photon
       p4photon_init.SetXYZT(riter->getBeam().getMomentum().getPx(),
                             riter->getBeam().getMomentum().getPy(),
                             riter->getBeam().getMomentum().getPz(),
                             riter->getBeam().getMomentum().getE());
       if (debug)
          std::cout << "photon.        : (" << setw(5) 
                    << p4photon_init.X() << ", " << setw(5) 
                    << p4photon_init.Y() << ", " << p4photon_init.Z() 
                    << ", " << setw(5) << p4photon_init.E() << ")"
                    << std::endl;

       hddm_s::VertexList vs = res().getVertices();
       hddm_s::VertexList::iterator viter;
       for (viter = vs.begin(); viter != vs.end(); ++viter) {
          hddm_s::OriginList orig = viter->getOrigins();
          hddm_s::ProductList ps = viter->getProducts();
          if (ps.size() > 0 && orig.size() > 0) {
             hddm_s::ProductList::iterator piter;
             for (piter = ps.begin(); piter != ps.end(); ++piter) {
                int type     = (int)piter->getType();
                int id       = piter->getId();
                int parentid = piter->getParentid();
                double E  = piter->getMomentum().getE();
                double px = piter->getMomentum().getPx();
                double py = piter->getMomentum().getPy();
                double pz = piter->getMomentum().getPz();
                double mass = sqrt(E*E - (px*px + py*py + pz*pz));
                if (! finite(mass))
                   mass = 0.0;
                if (debug) {
                std::cout << "\t\t\t--- new particle ---" << std::endl;
                   std::cout << "\t\t\ttype = " << type << " " 
                             << ParticleType((Particle_t)type) 
                             << " parentid = " << parentid << std::endl;
                   std::cout << "\t\t\t(" << px << "," << py << "," << pz << ","
                             << E << "), mass = " << mass << std::endl;
                   std::cout << "\t\t\t(" << orig().getVx() << "," 
                             << orig().getVy() << "," << orig().getVz() << ","
                             << orig().getT() << ")" << std::endl;
                }

                // Find proton
                if (type == Proton) {
                   foundProton = true;
                   idProton = id;
                   p4proton.SetXYZT(px,py,pz,E);
                }

                // Find pi+
                if (type == PiPlus) {
                   foundPip = true;
                   idPip = id;
                   p4pip.SetXYZT(px,py,pz,E);
                }

                // Find pi-
                if (type == PiMinus) {
                   foundPim = true;
                   idPim = id;
                   p4pim.SetXYZT(px,py,pz,E);
                }

                // Find pi0
                if (type == Pi0) {
                   foundPi0 = true;
                   idPi0 = id;
                   p4pi0.SetXYZT(px,py,pz,E);
                }

             } // end of loop over products
          } // end of having products and origin
       } // end of loop over vertices
    } // end of loop over reaction

    // If we found a proton decaying from Lambda,
    // save this event
    if (foundProton && foundPip && foundPim && foundPi0) {

      // make sure that total 4-mom of p, pi+, pi- are close to
      // initial state
       p4diff = p4photon_init + p4proton_init - p4proton - 
                p4pip - p4pim - p4pi0;

      // Make sure that remaining |p| < 1 MeV, M < 0.1 MeV
      if (p4diff.P() < 0.001 && p4diff.M() < 0.0001)
         return true;
    } // found proton, pi+, pi-
  } // end of select_type 7
  //__________________________________________________________________________________________________

  // default is to return false
  return false;
}
