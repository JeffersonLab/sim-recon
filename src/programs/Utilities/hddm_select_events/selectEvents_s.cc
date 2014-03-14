// $Id$
//
// Created Oct 10, 2013  David Lawrence

#include "hddm_select_events.h"
#include "particleType.h"

using namespace std;

extern TRandom2 *rndm;

//-----------------------------------------------------------------
// selectEvent_s
//-----------------------------------------------------------------
bool selectEvent_s(int select_type, s_HDDM_t* hddm_s, int nevents, bool debug){

  // Select Lambda -> p pi- events
  if(select_type==1){

    s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
    if(!PE){
      cout << "PhysicsEvent not given" << endl;
      return -1;
    }
    
    for(unsigned int i=0; i<PE->mult; i++){
      // ------------ Reactions --------------
      s_Reactions_t *reactions=PE->in[i].reactions;
      if(!reactions)continue;

      bool foundLambda = false;
      int  idLambda    = -999;
      int  idProton    = -999;
      bool foundProton = false;
      if(debug) cout << "------------------------------------------------------------------" << endl;
      
      for(unsigned int j=0; j<reactions->mult; j++){
	s_Vertices_t *vertices = reactions->in[j].vertices;
	if(vertices){
	  for(unsigned int k=0; k<vertices->mult; k++){
	    s_Origin_t *origin = vertices->in[k].origin;
	    s_Products_t *products = vertices->in[k].products;
	    if(products && origin){
	      for(unsigned int m=0;m<products->mult;m++){

		if(m==0){
		  if(debug) cout << setw(6) << nevents << " i = " << i << " (PhysicsEvent) j = " << j << " (reactions) k = " << k << " (vertices)" << endl;
		}
		s_Product_t *product = &products->in[m];
		
		int type     = (int)product->type;
		int id       = product->id;
		int parentid = product->parentid;
		double E  = product->momentum->E;
		double px = product->momentum->px;
		double py = product->momentum->py;
		double pz = product->momentum->pz;
		double mass = sqrt(E*E - (px*px + py*py + pz*pz));
		if(!finite(mass))mass = 0.0;
		if(debug){
		  cout << "\t\t\t--- m = " << m << " ---" << endl;
		  cout << "\t\t\ttype = " << type << " " << ParticleType((Particle_t)type) << " parentid = " << parentid << endl;
		  cout << "\t\t\t(" << px << "," << py << "," << pz << "," << E << "), mass = " << mass << endl;
		  cout << "\t\t\t(" << origin->vx << "," << origin->vy << "," << origin->vz << "," << origin->t << ")" << endl;
		}

		// Find Lambda
		if(type == Lambda){
		  foundLambda = true;
		  idLambda = id;
		}

		// Find decay products of Lambda
		if(foundLambda){
		  // if parent was Lambda
		  if(parentid == idLambda){
		    // find what type it is
		    if(type == Proton){
		      foundProton = true;
		      idProton = id;
		    }
		  }
		}

		// If we found a proton decaying from Lambda,
		// save this event
		if(foundProton){
		  return true;
		}

	      } // end of loop over products
	    } // end of having products and origin
	  } // end of loop over vertices
	} // end of having vertices
      } // end of loop over reaction
    } // end of loop over PhysicsEvent

  } // end of select_type 1
  //__________________________________________________________________________________________________

  // select Lambda -> p pi-, eta -> gamma gamma events
  else if(select_type==2){
    
    s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
    if(!PE){
      cout << "PhysicsEvent not given" << endl;
      return -1;
    }
    
    for(unsigned int i=0; i<PE->mult; i++){
      // ------------ Reactions --------------
      s_Reactions_t *reactions=PE->in[i].reactions;
      if(!reactions)continue;

      bool foundLambda = false;
      int  idLambda    = -999;
      int  idProton    = -999;
      bool foundProton = false;

      bool foundEta = false;
      int  idEta    = -999;
      int  idPhoton[2];
      int  nPhoton = 0;
      if(debug) cout << "------------------------------------------------------------------" << endl;
      
      for(unsigned int j=0; j<reactions->mult; j++){
	s_Vertices_t *vertices = reactions->in[j].vertices;
	if(vertices){
	  for(unsigned int k=0; k<vertices->mult; k++){
	    s_Origin_t *origin = vertices->in[k].origin;
	    s_Products_t *products = vertices->in[k].products;
	    if(products && origin){
	      for(unsigned int m=0;m<products->mult;m++){

		if(m==0){
		  if(debug) cout << setw(6) << nevents << " i = " << i << " (PhysicsEvent) j = " << j << " (reactions) k = " << k << " (vertices)" << endl;
		}
		s_Product_t *product = &products->in[m];
		
		int type     = (int)product->type;
		int id       = product->id;
		int parentid = product->parentid;
		double E  = product->momentum->E;
		double px = product->momentum->px;
		double py = product->momentum->py;
		double pz = product->momentum->pz;
		double mass = sqrt(E*E - (px*px + py*py + pz*pz));
		if(!finite(mass))mass = 0.0;
		if(debug){
		  cout << "\t\t\t--- m = " << m << " ---" << endl;
		  cout << "\t\t\ttype = " << type << " " << ParticleType((Particle_t)type) << " parentid = " << parentid << endl;
		  cout << "\t\t\t(" << px << "," << py << "," << pz << "," << E << "), mass = " << sqrt(E*E - px*px-py*py-pz*pz) << endl;
		  cout << "\t\t\t(" << origin->vx << "," << origin->vy << "," << origin->vz << "," << origin->t << ")" << endl;
		}

		// Find Lambda
		if(type == Lambda){
		  foundLambda = true;
		  idLambda = id;
		}

		// Find decay products of Lambda
		if(foundLambda){
		  // if parent was Lambda
		  if(parentid == idLambda){
		    // find what type it is
		    if(type == Proton){
		      foundProton = true;
		      idProton = id;
		    }
		  }
		}

		// Find eta
		if(type == Eta){
		  foundEta = true;
		  idEta = id;
		}

		// Find decay products of Eta
		if(foundEta){
		  // if parent was Eta
		  if(parentid == idEta){
		    // find what type it is
		    if(type == Gamma){
		      idPhoton[nPhoton] = id;
		      nPhoton++;
		    }
		  }
		}

		// If we found a proton decaying from Lambda,
		// and we found 2 photons from eta decay,
		// save this event
		if(foundProton && nPhoton==2){
		  return true;
		}

	      } // end of loop over products
	    } // end of having products and origin
	  } // end of loop over vertices
	} // end of having vertices
      } // end of loop over reaction
    } // end of loop over PhysicsEvent

  }
  //__________________________________________________________________________________________________

  // Select Lambda -> p pi- events, decay length is less than VERTEXDIFFMAX
  else if(select_type==3){

    const double VERTEXDIFFMAX = 2.0;

    s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
    if(!PE){
      cout << "PhysicsEvent not given" << endl;
      return -1;
    }
    
    for(unsigned int i=0; i<PE->mult; i++){
      // ------------ Reactions --------------
      s_Reactions_t *reactions=PE->in[i].reactions;
      if(!reactions)continue;

      bool foundLambda = false;
      int  idLambda    = -999;
      int  idProton    = -999;
      bool foundProton = false;
      double vertexLambda[3];
      double vertexProton[3];
      for(int v=0;v<3;v++){
	vertexLambda[v] = -999;
	vertexProton[v] = 999;
      }
      double vertexdiff = 999;
      if(debug) cout << "------------------------------------------------------------------" << endl;
      
      for(unsigned int j=0; j<reactions->mult; j++){
	s_Vertices_t *vertices = reactions->in[j].vertices;
	if(vertices){
	  for(unsigned int k=0; k<vertices->mult; k++){
	    s_Origin_t *origin = vertices->in[k].origin;
	    s_Products_t *products = vertices->in[k].products;
	    if(products && origin){
	      for(unsigned int m=0;m<products->mult;m++){

		if(m==0){
		  if(debug) cout << setw(6) << nevents << " i = " << i << " (PhysicsEvent) j = " << j << " (reactions) k = " << k << " (vertices)" << endl;
		}
		s_Product_t *product = &products->in[m];
		
		int type     = (int)product->type;
		int id       = product->id;
		int parentid = product->parentid;
		double E  = product->momentum->E;
		double px = product->momentum->px;
		double py = product->momentum->py;
		double pz = product->momentum->pz;
		double mass = sqrt(E*E - (px*px + py*py + pz*pz));
		if(!finite(mass))mass = 0.0;
		if(debug){
		  cout << "\t\t\t--- m = " << m << " ---" << endl;
		  cout << "\t\t\ttype = " << type << " " << ParticleType((Particle_t)type) << " parentid = " << parentid << endl;
		  cout << "\t\t\t(" << px << "," << py << "," << pz << "," << E << "), mass = " << sqrt(E*E - px*px-py*py-pz*pz) << endl;
		  cout << "\t\t\t(" << origin->vx << "," << origin->vy << "," << origin->vz << "," << origin->t << ")" << endl;
		}

		// Find Lambda
		if(type == Lambda){
		  foundLambda = true;
		  idLambda = id;
		  vertexLambda[0] = origin->vx;
		  vertexLambda[1] = origin->vy;
		  vertexLambda[2] = origin->vz;
		}

		// Find decay products of Lambda
		if(foundLambda){
		  // if parent was Lambda
		  if(parentid == idLambda){
		    // find what type it is
		    if(type == Proton){
		      foundProton = true;
		      idProton = id;
		      vertexProton[0] = origin->vx;
		      vertexProton[1] = origin->vy;
		      vertexProton[2] = origin->vz;
		      vertexdiff = sqrt(pow(vertexLambda[0] - vertexProton[0],2.) + pow(vertexLambda[1] - vertexProton[1],2.) + pow(vertexLambda[2] - vertexProton[2],2.));
		    }
		  }
		}

		// If we found a proton decaying from Lambda,
		// and the vertex difference is less than VERTEXDIFFMAX,
		// save this event
		if(foundProton && vertexdiff < VERTEXDIFFMAX){
		  return true;
		}

	      } // end of loop over products
	    } // end of having products and origin
	  } // end of loop over vertices
	} // end of having vertices
      } // end of loop over reaction
    } // end of loop over PhysicsEvent

  } // end of select_type 3
  //__________________________________________________________________________________________________

  // Select Lambda -> p pi- events, carve out weak decay asymmetry
  else if(select_type==4){

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

    s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
    if(!PE){
      cout << "PhysicsEvent not given" << endl;
      return -1;
    }
    
    for(unsigned int i=0; i<PE->mult; i++){
      // ------------ Reactions --------------
      s_Reactions_t *reactions=PE->in[i].reactions;
      if(!reactions)continue;

      if(debug) cout << "------------------------------------------------------------------" << endl;

      for(unsigned int j=0; j<reactions->mult; j++){

	// Get initial photon
	p4photon_init.SetXYZT(reactions->in[j].beam->momentum->px,reactions->in[j].beam->momentum->py,reactions->in[j].beam->momentum->pz,reactions->in[j].beam->momentum->E);
	if(debug){
	  cout << "photon.        : (" << setw(5) << p4photon_init.X() << ", " << setw(5) << p4photon_init.Y() << ", " << p4photon_init.Z() << ", " << setw(5) << p4photon_init.E() << ")" << endl;
	}
      

	s_Vertices_t *vertices = reactions->in[j].vertices;
	if(vertices){
	  for(unsigned int k=0; k<vertices->mult; k++){
	    s_Origin_t *origin = vertices->in[k].origin;
	    s_Products_t *products = vertices->in[k].products;
	    if(products && origin){
	      for(unsigned int m=0;m<products->mult;m++){

		if(m==0){
		  if(debug) cout << setw(6) << nevents << " i = " << i << " (PhysicsEvent) j = " << j << " (reactions) k = " << k << " (vertices)" << endl;
		}
		s_Product_t *product = &products->in[m];
		
		int type     = (int)product->type;
		int id       = product->id;
		int parentid = product->parentid;
		double E  = product->momentum->E;
		double px = product->momentum->px;
		double py = product->momentum->py;
		double pz = product->momentum->pz;
		double mass = sqrt(E*E - (px*px + py*py + pz*pz));
		if(!finite(mass))mass = 0.0;
		if(debug){
		  cout << "\t\t\t--- m = " << m << " ---" << endl;
		  cout << "\t\t\ttype = " << type << " " << ParticleType((Particle_t)type) << " parentid = " << parentid << endl;
		  cout << "\t\t\t(" << px << "," << py << "," << pz << "," << E << "), mass = " << sqrt(E*E - px*px-py*py-pz*pz) << endl;
		  cout << "\t\t\t(" << origin->vx << "," << origin->vy << "," << origin->vz << "," << origin->t << ")" << endl;
		}

		// Find K+, set up polarization axis
		if(type == KPlus){
		  foundKPlus = true;
		  p4kp.SetXYZT(px,py,pz,E);
		  normaldir = p4photon_init.Vect().Cross(p4kp.Vect()).Unit();
		  phi_n = normaldir.Phi();

		  if(debug){
		    cout << "K+ dir.    : (" << setw(5) << p4kp.X() << ", " << setw(5) << p4kp.Y() << ", " << p4kp.Z() << ")" << endl;
		    cout << "normal dir.: (" << setw(5) << normaldir.X() << ", " << setw(5) << normaldir.Y() << ", " << normaldir.Z() << ")" << endl;
		  }

		}

		// Find Lambda
		if(type == Lambda){
		  foundLambda = true;
		  idLambda = id;
		  p4Lambda.SetXYZT(px,py,pz,E);
		}

		// Find decay products of Lambda
		if(foundLambda){
		  // if parent was Lambda
		  if(parentid == idLambda){
		    // find what type it is
		    if(type == Proton){
		      foundProton = true;
		      idProton = id;
		      p4proton.SetXYZT(px,py,pz,E);
		      if(debug) cout << "proton : (" << setw(5) << p4proton.X() << ", " << setw(5) << p4proton.Y() << ", " << p4proton.Z() << ")" << endl;
		    }
		  }
		}
	      } // end of loop over products
	    } // end of having products and origin
	  } // end of loop over vertices
	} // end of having vertices
      } // end of loop over reaction
    } // end of loop over PhysicsEvent

    // If we found a proton decaying from Lambda,
    // save this event
    if(foundKPlus && foundProton){
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
      if(debug){
	cout << "random number: " << rand << endl;
	cout << "compare to   : " << (1. + alpha * costheta_proton) / (1. + alpha) << endl;
      }
      if(rand < (1. + alpha * costheta_proton) / (1. + alpha)){
	cout << " --- " << costheta_proton << endl;
	return true;
      }
    } // found K+ and proton

  } // end of select_type 4

  // Select p pi+ pi- events
  else if(select_type==5){

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

    s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
    if(!PE){
      cout << "PhysicsEvent not given" << endl;
      return -1;
    }
    
    for(unsigned int i=0; i<PE->mult; i++){
      // ------------ Reactions --------------
      s_Reactions_t *reactions=PE->in[i].reactions;
      if(!reactions)continue;

      if(debug) cout << "------------------------------------------------------------------" << endl;

      for(unsigned int j=0; j<reactions->mult; j++){

	// Get initial photon
	p4photon_init.SetXYZT(reactions->in[j].beam->momentum->px,reactions->in[j].beam->momentum->py,reactions->in[j].beam->momentum->pz,reactions->in[j].beam->momentum->E);
	if(debug){
	  cout << "photon.        : (" << setw(5) << p4photon_init.X() << ", " << setw(5) << p4photon_init.Y() << ", " << p4photon_init.Z() << ", " << setw(5) << p4photon_init.E() << ")" << endl;
	}

	s_Vertices_t *vertices = reactions->in[j].vertices;
	if(vertices){
	  for(unsigned int k=0; k<vertices->mult; k++){
	    s_Origin_t *origin = vertices->in[k].origin;
	    s_Products_t *products = vertices->in[k].products;
	    if(products && origin){
	      for(unsigned int m=0;m<products->mult;m++){

		if(m==0){
		  if(debug) cout << setw(6) << nevents << " i = " << i << " (PhysicsEvent) j = " << j << " (reactions) k = " << k << " (vertices)" << endl;
		}
		s_Product_t *product = &products->in[m];
		
		int type     = (int)product->type;
		int id       = product->id;
		int parentid = product->parentid;
		double E  = product->momentum->E;
		double px = product->momentum->px;
		double py = product->momentum->py;
		double pz = product->momentum->pz;
		double mass = sqrt(E*E - (px*px + py*py + pz*pz));
		if(!finite(mass))mass = 0.0;
		if(debug){
		  cout << "\t\t\t--- m = " << m << " ---" << endl;
		  cout << "\t\t\ttype = " << type << " " << ParticleType((Particle_t)type) << " parentid = " << parentid << endl;
		  cout << "\t\t\t(" << px << "," << py << "," << pz << "," << E << "), mass = " << sqrt(E*E - px*px-py*py-pz*pz) << endl;
		  cout << "\t\t\t(" << origin->vx << "," << origin->vy << "," << origin->vz << "," << origin->t << ")" << endl;
		}

		// Find proton
		if(type == Proton){
		  foundProton = true;
		  idProton = id;
		  p4proton.SetXYZT(px,py,pz,E);

		}

		// Find pi+
		if(type == PiPlus){
		  foundPip = true;
		  idPip = id;
		  p4pip.SetXYZT(px,py,pz,E);
		}

		// Find pi-
		if(type == PiMinus){
		  foundPim = true;
		  idPim = id;
		  p4pim.SetXYZT(px,py,pz,E);
		}
	      } // end of loop over products
	    } // end of having products and origin
	  } // end of loop over vertices
	} // end of having vertices
      } // end of loop over reaction
    } // end of loop over PhysicsEvent

    // If we found a proton decaying from Lambda,
    // save this event
    if(foundProton && foundPip && foundPim){

      // make sure that total 4-mom of p, pi+, pi- are close to
      // initial state
      p4diff = p4photon_init + p4proton_init - p4proton - p4pip - p4pim;

      // Make sure that remaining |p| < 1 MeV, M < 0.1 MeV
      if(p4diff.P() < 0.001 && p4diff.M() < 0.0001){
	return true;
      }
    } // found proton, pi+, pi-

  } // end of select_type 5

  // Select K+ K+ Xi- -> K+ K+ p pi- pi- events
  else if(select_type==6){

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

    s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
    if(!PE){
      cout << "PhysicsEvent not given" << endl;
      return -1;
    }
    
    for(unsigned int i=0; i<PE->mult; i++){
      // ------------ Reactions --------------
      s_Reactions_t *reactions=PE->in[i].reactions;
      if(!reactions)continue;

      if(debug) cout << "------------------------------------------------------------------" << endl;

      for(unsigned int j=0; j<reactions->mult; j++){

	// Get initial photon
	p4photon_init.SetXYZT(reactions->in[j].beam->momentum->px,reactions->in[j].beam->momentum->py,reactions->in[j].beam->momentum->pz,reactions->in[j].beam->momentum->E);
	if(debug){
	  cout << "photon.        : (" << setw(5) << p4photon_init.X() << ", " << setw(5) << p4photon_init.Y() << ", " << p4photon_init.Z() << ", " << setw(5) << p4photon_init.E() << ")" << endl;
	}

	s_Vertices_t *vertices = reactions->in[j].vertices;
	if(vertices){
	  for(unsigned int k=0; k<vertices->mult; k++){
	    s_Origin_t *origin = vertices->in[k].origin;
	    s_Products_t *products = vertices->in[k].products;
	    if(products && origin){
	      for(unsigned int m=0;m<products->mult;m++){

		if(m==0){
		  if(debug) cout << setw(6) << nevents << " i = " << i << " (PhysicsEvent) j = " << j << " (reactions) k = " << k << " (vertices)" << endl;
		}
		s_Product_t *product = &products->in[m];
		
		int type     = (int)product->type;
		int id       = product->id;
		int parentid = product->parentid;
		double E  = product->momentum->E;
		double px = product->momentum->px;
		double py = product->momentum->py;
		double pz = product->momentum->pz;
		double mass = sqrt(E*E - (px*px + py*py + pz*pz));
		if(!finite(mass))mass = 0.0;
		if(debug){
		  cout << "\t\t\t--- m = " << m << " ---" << endl;
		  cout << "\t\t\t" << ParticleType((Particle_t)type) << " type = " << type << " parentid = " << parentid << endl;
		  cout << "\t\t\t(" << px << "," << py << "," << pz << "," << E << "), mass = " << sqrt(E*E - px*px-py*py-pz*pz) << endl;
		  cout << "\t\t\t(" << origin->vx << "," << origin->vy << "," << origin->vz << "," << origin->t << ")" << endl;
		}

		// Find first K+
		if(nKp==0 && type == KPlus){
		  nKp++;
		  p4kp1.SetXYZT(px,py,pz,E);
		  if(debug)
		    cout << "K+ 1  : (" << p4kp1.X() << ", " << p4kp1.Y() << ", " << p4kp1.Z() << ", " << p4kp1.M() << ")" << endl;

		  // Need this break since without it we would count the
		  // the first K+ as the second K+ too
		  continue;
		}

		// Find second K+
		if(nKp==1 && type == KPlus){
		  nKp++;
		  p4kp2.SetXYZT(px,py,pz,E);
		  if(debug)
		    cout << "K+ 2  : (" << p4kp2.X() << ", " << p4kp2.Y() << ", " << p4kp2.Z() << ", " << p4kp2.M() << ")" << endl;
		}

		// Find Xi-
		if(type == XiMinus){
		  nXiMinus++;
		  idXiMinus = id;
		  p4XiMinus.SetXYZT(px,py,pz,E);
		}

		// Find pi- from Xi-
		if(type == PiMinus && parentid == idXiMinus){
		  nPimFromXiMinus++;
		  p4pimFromXiMinus.SetXYZT(px,py,pz,E);
		}

		// Find Lambda from Xi-
		if(type == Lambda && parentid == idXiMinus){
		  nLambda++;
		  idLambda = id;
		  p4Lambda.SetXYZT(px,py,pz,E);
		}

		// Find proton from Lambda
		if(type == Proton && parentid == idLambda){
		  nProton++;
		  p4proton.SetXYZT(px,py,pz,E);
		}

		// Find pi- from Lambda
		if(type == PiMinus && parentid == idLambda){
		  nPimFromLambda++;
		  p4pimFromLambda.SetXYZT(px,py,pz,E);
		}
	      
	      } // end of loop over products
	    } // end of having products and origin
	  } // end of loop over vertices
	} // end of having vertices
      } // end of loop over reaction
    } // end of loop over PhysicsEvent

    // If we found a proton decaying from Lambda,
    // save this event
    if(nKp==2 && nXiMinus==1 && nLambda==1 && nProton==1 && nPimFromXiMinus==1 && nPimFromLambda==1){
      cout << "found all particles" << endl;

      // make sure that total 4-mom of p, pi+, pi- are close to
      // initial state
      p4diff  = p4photon_init + p4proton_init - p4kp1 - p4kp2 - p4proton - p4pimFromXiMinus - p4pimFromLambda;
      p4total = p4kp1 + p4kp2 + p4proton + p4pimFromXiMinus + p4pimFromLambda;

      if(debug){
	cout << "photon: (" << p4photon_init.X() << ", " << p4photon_init.Y() << ", " << p4photon_init.Z() << ", " << p4photon_init.M() << ")" << endl;
	cout << "proton: (" << p4proton_init.X() << ", " << p4proton_init.Y() << ", " << p4proton_init.Z() << ", " << p4proton_init.M() << ")" << endl;
	cout << "Xi-   : (" << p4XiMinus.X() << ", " << p4XiMinus.Y() << ", " << p4XiMinus.Z() << ", " << p4XiMinus.M() << ")" << endl;
	cout << "Lambda: (" << p4Lambda.X() << ", " << p4Lambda.Y() << ", " << p4Lambda.Z() << ", " << p4Lambda.M() << ")" << endl;
	cout << "K+ 1  : (" << p4kp1.X() << ", " << p4kp1.Y() << ", " << p4kp1.Z() << ", " << p4kp1.M() << ")" << endl;
	cout << "K+ 2  : (" << p4kp2.X() << ", " << p4kp2.Y() << ", " << p4kp2.Z() << ", " << p4kp2.M() << ")" << endl;
	cout << "proton: (" << p4proton.X() << ", " << p4proton.Y() << ", " << p4proton.Z() << ", " << p4proton.M() << ")" << endl;
	cout << "pi- 1 : (" << p4pimFromXiMinus.X() << ", " << p4pimFromXiMinus.Y() << ", " << p4pimFromXiMinus.Z() << ", " << p4pimFromXiMinus.M() << ")" << endl;
	cout << "pi- 2 : (" << p4pimFromLambda.X() << ", " << p4pimFromLambda.Y() << ", " << p4pimFromLambda.Z() << ", " << p4pimFromLambda.M() << ")" << endl;
      }

      // Make sure that remaining |p| < 1 MeV, M < 0.1 MeV
      if(p4diff.P() < 0.001 && p4diff.M() < 0.0001){
	cout << "returning true" << endl;
	return true;
      }
      else{
	cout << "Added 4-mom of final state particles does not match generated by more than 0.1 MeV!!!!" << endl;
	cout << "diff in |p| = " << p4diff.P()* 1000. << " MeV/c" << endl;
	cout << "diff  = (" << p4diff.X() << ", " << p4diff.Y() << ", " << p4diff.Z() << ", " << p4diff.E() << ")" << endl;
	cout << "total = (" << p4total.X() << ", " << p4total.Y() << ", " << p4total.Z() << ", " << p4total.E() << ")" << endl;
	abort();
      }
    } // found proton, pi+, pi-

  } // end of select_type 6

  // Select p pi+ pi- pi0 events
  else if(select_type==7){

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

    s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
    if(!PE){
      cout << "PhysicsEvent not given" << endl;
      return -1;
    }
    
    for(unsigned int i=0; i<PE->mult; i++){
      // ------------ Reactions --------------
      s_Reactions_t *reactions=PE->in[i].reactions;
      if(!reactions)continue;

      if(debug) cout << "------------------------------------------------------------------" << endl;

      for(unsigned int j=0; j<reactions->mult; j++){

	// Get initial photon
	p4photon_init.SetXYZT(reactions->in[j].beam->momentum->px,reactions->in[j].beam->momentum->py,reactions->in[j].beam->momentum->pz,reactions->in[j].beam->momentum->E);
	if(debug){
	  cout << "photon.        : (" << setw(5) << p4photon_init.X() << ", " << setw(5) << p4photon_init.Y() << ", " << p4photon_init.Z() << ", " << setw(5) << p4photon_init.E() << ")" << endl;
	}

	s_Vertices_t *vertices = reactions->in[j].vertices;
	if(vertices){
	  for(unsigned int k=0; k<vertices->mult; k++){
	    s_Origin_t *origin = vertices->in[k].origin;
	    s_Products_t *products = vertices->in[k].products;
	    if(products && origin){
	      for(unsigned int m=0;m<products->mult;m++){

		if(m==0){
		  if(debug) cout << setw(6) << nevents << " i = " << i << " (PhysicsEvent) j = " << j << " (reactions) k = " << k << " (vertices)" << endl;
		}
		s_Product_t *product = &products->in[m];
		
		int type     = (int)product->type;
		int id       = product->id;
		int parentid = product->parentid;
		double E  = product->momentum->E;
		double px = product->momentum->px;
		double py = product->momentum->py;
		double pz = product->momentum->pz;
		double mass = sqrt(E*E - (px*px + py*py + pz*pz));
		if(!finite(mass))mass = 0.0;
		if(debug){
		  cout << "\t\t\t--- m = " << m << " ---" << endl;
		  cout << "\t\t\ttype = " << type << " " << ParticleType((Particle_t)type) << " parentid = " << parentid << endl;
		  cout << "\t\t\t(" << px << "," << py << "," << pz << "," << E << "), mass = " << sqrt(E*E - px*px-py*py-pz*pz) << endl;
		  cout << "\t\t\t(" << origin->vx << "," << origin->vy << "," << origin->vz << "," << origin->t << ")" << endl;
		}

		// Find proton
		if(type == Proton){
		  foundProton = true;
		  idProton = id;
		  p4proton.SetXYZT(px,py,pz,E);

		}

		// Find pi+
		if(type == PiPlus){
		  foundPip = true;
		  idPip = id;
		  p4pip.SetXYZT(px,py,pz,E);
		}

		// Find pi-
		if(type == PiMinus){
		  foundPim = true;
		  idPim = id;
		  p4pim.SetXYZT(px,py,pz,E);
		}

		// Find pi0
		if(type == Pi0){
		  foundPi0 = true;
		  idPi0 = id;
		  p4pi0.SetXYZT(px,py,pz,E);
		}

	      } // end of loop over products
	    } // end of having products and origin
	  } // end of loop over vertices
	} // end of having vertices
      } // end of loop over reaction
    } // end of loop over PhysicsEvent

    // If we found a proton decaying from Lambda,
    // save this event
    if(foundProton && foundPip && foundPim && foundPi0){

      // make sure that total 4-mom of p, pi+, pi- are close to
      // initial state
      p4diff = p4photon_init + p4proton_init - p4proton - p4pip - p4pim - p4pi0;

      // Make sure that remaining |p| < 1 MeV, M < 0.1 MeV
      if(p4diff.P() < 0.001 && p4diff.M() < 0.0001){
	return true;
      }
    } // found proton, pi+, pi-

  } // end of select_type 7

  //__________________________________________________________________________________________________


  // default is to return false
  return false;
}
