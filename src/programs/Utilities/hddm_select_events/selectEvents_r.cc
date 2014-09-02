// $Id$
//
// Created Oct 10, 2013  David Lawrence

#include "hddm_select_events.h"
#include <HDDM/hddm_r.hpp>
#include "particleType.h"

using namespace std;
using std::distance;

//-----------------------------------------------------------------
// selectEvent_s
//-----------------------------------------------------------------
bool selectEvent_r(int select_type, hddm_r::HDDM &record, int nevents, bool debug){

  // Select Lambda -> p pi- events
  if(select_type==1){

    bool foundLambda = false;
    int  idLambda    = -999;
    int  idProton    = -999;
    bool foundProton = false;
    if(debug) cout << "event: " << nevents << "----------------------------------------------------------------" << endl;

    hddm_r::ReconstructedPhysicsEvent &re = record.getReconstructedPhysicsEvent();
    hddm_r::ReactionList reactions = re.getReactions();
    // Loop over ReactionList
    hddm_r::ReactionList::iterator iter;
    for (iter = reactions.begin(); iter != reactions.end(); ++iter) {
      int ireaction = std::distance(reactions.begin(), iter);
      
      hddm_r::VertexList vertices = iter->getVertices();
      hddm_r::VertexList::iterator iter_vertex;
      // Loop over VertexList
      for (iter_vertex = vertices.begin(); iter_vertex != vertices.end(); ++iter_vertex){
	int ivertex = std::distance(vertices.begin(), iter_vertex);

	hddm_r::ProductList products = iter_vertex->getProducts();
	hddm_r::ProductList::iterator iter_product;
	// Loop over ProductList
	for (iter_product = products.begin(); iter_product != products.end(); ++iter_product){
	  int iparticle = std::distance(products.begin(), iter_product);

	  if(debug) cout << "reaction: " << ireaction << " vertex: " << ivertex << " particle: " << iparticle
			 << " id: " << setw(3) << iter_product->getId() << " parent: " << iter_product->getParentId()
			 << " PDG type: " << setw(7) << iter_product->getPdgtype()
			 << " type: " << ParticleType((Particle_t)PDGtoPType((Particle_t)iter_product->getPdgtype()))
			 << " vertex: (" << iter_vertex->getOrigin().getVx() << ", " << iter_vertex->getOrigin().getVy() << ", " << iter_vertex->getOrigin().getVz() << ")"
			 << endl;
	  
	  // Find Lambda
	  if(iter_product->getPdgtype() == PDGtype(Lambda)){
	    foundLambda = true;
	    idLambda = iter_product->getId();
	  }

	  // Find decay products of Lambda
	  if(foundLambda){
	    // if parent was Lambda
	    if(iter_product->getParentId() == idLambda){
	      // find what type it is
	      if(iter_product->getPdgtype() == PDGtype(Proton)){
		foundProton = true;
		idProton = iter_product->getId();
	      }
	    }
	  }

	  // If we found a proton decaying from Lambda,
	  // save this event
	  if(foundProton){
	    return true;
	  }

	} // end of loop over ProductList
      }	// end of loop over VertexList
    } // end of loop over ReactiionList
  }
  //__________________________________________________________________________________________________

  // select Lambda -> p pi-, eta -> gamma gamma events
  else if(select_type==2){

    bool foundLambda = false;
    int  idLambda    = -999;
    int  idProton    = -999;
    bool foundProton = false;

    bool foundEta = false;
    int  idEta    = -999;
    int  idPhoton[2];
    int  nPhoton = 0;

    if(debug) cout << "event: " << nevents << "----------------------------------------------------------------" << endl;

    hddm_r::ReconstructedPhysicsEvent &re = record.getReconstructedPhysicsEvent();
    hddm_r::ReactionList reactions = re.getReactions();
    // Loop over ReactionList
    hddm_r::ReactionList::iterator iter;
    for (iter = reactions.begin(); iter != reactions.end(); ++iter) {
      int ireaction = std::distance(reactions.begin(), iter);
      
      hddm_r::VertexList vertices = iter->getVertices();
      hddm_r::VertexList::iterator iter_vertex;
      // Loop over VertexList
      for (iter_vertex = vertices.begin(); iter_vertex != vertices.end(); ++iter_vertex){
	int ivertex = std::distance(vertices.begin(), iter_vertex);

	hddm_r::ProductList products = iter_vertex->getProducts();
	hddm_r::ProductList::iterator iter_product;
	// Loop over ProductList
	for (iter_product = products.begin(); iter_product != products.end(); ++iter_product){
	  int iparticle = std::distance(products.begin(), iter_product);

	  if(debug) cout << "reaction: " << ireaction << " vertex: " << ivertex << " particle: " << iparticle
			 << " id: " << setw(3) << iter_product->getId() << " parent: " << iter_product->getParentId()
			 << " PDG type: " << setw(7) << iter_product->getPdgtype()
			 << " type: " << ParticleType((Particle_t)PDGtoPType((Particle_t)iter_product->getPdgtype()))
			 << " vertex: (" << iter_vertex->getOrigin().getVx() << ", " << iter_vertex->getOrigin().getVy() << ", " << iter_vertex->getOrigin().getVz() << ")"
			 << endl;
	  
	  // Find Lambda
	  if(iter_product->getPdgtype() == PDGtype(Lambda)){
	    foundLambda = true;
	    idLambda = iter_product->getId();
	  }

	  // Find decay products of Lambda
	  if(foundLambda){
	    // if parent was Lambda
	    if(iter_product->getParentId() == idLambda){
	      // find what type it is
	      if(iter_product->getPdgtype() == PDGtype(Proton)){
		foundProton = true;
		idProton = iter_product->getId();
	      }
	    }
	  }

	  // Find Eta
	  if(iter_product->getPdgtype() == PDGtype(Eta)){
	    foundEta = true;
	    idEta = iter_product->getId();
	  }

	  // Find decay products of Eta
	  if(foundEta){
	    // if parent was Eta
	    if(iter_product->getParentId() == idEta){
	      // find what type it is
	      if(iter_product->getPdgtype() == PDGtype(Gamma)){
		idPhoton[nPhoton] = iter_product->getId();
		nPhoton++;
	      }
	    }
	  }

	  // If we found a proton decaying from Lambda,
	  // save this event
	  if(foundProton && nPhoton==2){
	    return true;
	  }

	} // end of loop over ProductList
      }	// end of loop over VertexList
    } // end of loop over ReactiionList
  }
  //__________________________________________________________________________________________________

  // Select Lambda -> p pi- events, decay length is less than VERTEXDIFFMAX
  else if(select_type==3){
    const double VERTEXDIFFMAX = 10.0;

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

    hddm_r::ReconstructedPhysicsEvent &re = record.getReconstructedPhysicsEvent();
    hddm_r::ReactionList reactions = re.getReactions();
    // Loop over ReactionList
    hddm_r::ReactionList::iterator iter;
    for (iter = reactions.begin(); iter != reactions.end(); ++iter) {
      int ireaction = std::distance(reactions.begin(), iter);
      
      hddm_r::VertexList vertices = iter->getVertices();
      hddm_r::VertexList::iterator iter_vertex;
      // Loop over VertexList
      for (iter_vertex = vertices.begin(); iter_vertex != vertices.end(); ++iter_vertex){
	int ivertex = std::distance(vertices.begin(), iter_vertex);

	hddm_r::ProductList products = iter_vertex->getProducts();
	hddm_r::ProductList::iterator iter_product;
	// Loop over ProductList
	for (iter_product = products.begin(); iter_product != products.end(); ++iter_product){
	  int iparticle = std::distance(products.begin(), iter_product);

	  if(debug) cout << "reaction: " << ireaction << " vertex: " << ivertex << " particle: " << iparticle
			 << " id: " << setw(3) << iter_product->getId() << " parent: " << iter_product->getParentId()
			 << " PDG type: " << setw(7) << iter_product->getPdgtype()
			 << " type: " << ParticleType((Particle_t)PDGtoPType((Particle_t)iter_product->getPdgtype()))
			 << " vertex: (" << iter_vertex->getOrigin().getVx() << ", " << iter_vertex->getOrigin().getVy() << ", " << iter_vertex->getOrigin().getVz() << ")"
			 << endl;
	  
	  // Find Lambda
	  if(iter_product->getPdgtype() == PDGtype(Lambda)){
	    foundLambda = true;
	    idLambda = iter_product->getId();
	    vertexLambda[0] = iter_vertex->getOrigin().getVx();
	    vertexLambda[1] = iter_vertex->getOrigin().getVy();
	    vertexLambda[2] = iter_vertex->getOrigin().getVz();
	  }

	  // Find decay products of Lambda
	  if(foundLambda){
	    // if parent was Lambda
	    if(iter_product->getParentId() == idLambda){
	      // find what type it is
	      if(iter_product->getPdgtype() == PDGtype(Proton)){
		foundProton = true;
		idProton = iter_product->getId();
		vertexProton[0] = iter_vertex->getOrigin().getVx();
		vertexProton[1] = iter_vertex->getOrigin().getVy();
		vertexProton[2] = iter_vertex->getOrigin().getVz();
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

	} // end of loop over ProductList
      }	// end of loop over VertexList
    } // end of loop over ReactiionList
  } // end of select_type 3
  //__________________________________________________________________________________________________

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

    hddm_r::ReconstructedPhysicsEvent &re = record.getReconstructedPhysicsEvent();
    hddm_r::ReactionList reactions = re.getReactions();
    // Loop over ReactionList
    hddm_r::ReactionList::iterator iter;
    for (iter = reactions.begin(); iter != reactions.end(); ++iter) {
      int ireaction = std::distance(reactions.begin(), iter);
      
      hddm_r::VertexList vertices = iter->getVertices();
      hddm_r::VertexList::iterator iter_vertex;
      // Loop over VertexList
      for (iter_vertex = vertices.begin(); iter_vertex != vertices.end(); ++iter_vertex){
	int ivertex = std::distance(vertices.begin(), iter_vertex);

	hddm_r::ProductList products = iter_vertex->getProducts();
	hddm_r::ProductList::iterator iter_product;
	// Loop over ProductList
	for (iter_product = products.begin(); iter_product != products.end(); ++iter_product){
	  int iparticle = std::distance(products.begin(), iter_product);

	  if(debug) cout << "reaction: " << ireaction << " vertex: " << ivertex << " particle: " << iparticle
			 << " id: " << setw(3) << iter_product->getId() << " parent: " << iter_product->getParentId()
			 << " PDG type: " << setw(7) << iter_product->getPdgtype()
			 << " type: " << ParticleType((Particle_t)PDGtoPType((Particle_t)iter_product->getPdgtype()))
			 << " vertex: (" << iter_vertex->getOrigin().getVx() << ", " << iter_vertex->getOrigin().getVy() << ", " << iter_vertex->getOrigin().getVz() << ")"
			 << " Ebeam: " << iter_vertex->getOrigin().getEbeam()
			 << endl;

	  p4photon_init.SetXYZT(0,0,iter_vertex->getOrigin().getEbeam(),iter_vertex->getOrigin().getEbeam());
	  
	  // Find proton
	  if(iter_product->getPdgtype() == PDGtype(Proton)){
	    foundProton = true;
	    idProton = iter_product->getId();
	    p4proton.SetXYZT(iter_product->getMomentum().getPx(),iter_product->getMomentum().getPy(),iter_product->getMomentum().getPz(),iter_product->getMomentum().getE());
	  }

	  // Find pi+
	  if(iter_product->getPdgtype() == PDGtype(PiPlus)){
	    foundPip = true;
	    idPip = iter_product->getId();
	    p4pip.SetXYZT(iter_product->getMomentum().getPx(),iter_product->getMomentum().getPy(),iter_product->getMomentum().getPz(),iter_product->getMomentum().getE());
	  }

	  // Find pi-
	  if(iter_product->getPdgtype() == PDGtype(PiMinus)){
	    foundPim = true;
	    idPim = iter_product->getId();
	    p4pim.SetXYZT(iter_product->getMomentum().getPx(),iter_product->getMomentum().getPy(),iter_product->getMomentum().getPz(),iter_product->getMomentum().getE());
	  }
	} // end of loop over ProductList
      }	// end of loop over VertexList
    } // end of loop over ReactiionList

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
  //__________________________________________________________________________________________________

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

    hddm_r::ReconstructedPhysicsEvent &re = record.getReconstructedPhysicsEvent();
    hddm_r::ReactionList reactions = re.getReactions();
    // Loop over ReactionList
    hddm_r::ReactionList::iterator iter;
    for (iter = reactions.begin(); iter != reactions.end(); ++iter) {
      int ireaction = std::distance(reactions.begin(), iter);
      
      hddm_r::VertexList vertices = iter->getVertices();
      hddm_r::VertexList::iterator iter_vertex;
      // Loop over VertexList
      for (iter_vertex = vertices.begin(); iter_vertex != vertices.end(); ++iter_vertex){
	int ivertex = std::distance(vertices.begin(), iter_vertex);

	hddm_r::ProductList products = iter_vertex->getProducts();
	hddm_r::ProductList::iterator iter_product;
	// Loop over ProductList
	for (iter_product = products.begin(); iter_product != products.end(); ++iter_product){
	  int iparticle = std::distance(products.begin(), iter_product);

	  if(debug) cout << "reaction: " << ireaction << " vertex: " << ivertex << " particle: " << iparticle
			 << " id: " << setw(3) << iter_product->getId() << " parent: " << iter_product->getParentId()
			 << " PDG type: " << setw(7) << iter_product->getPdgtype()
			 << " type: " << ParticleType((Particle_t)PDGtoPType((Particle_t)iter_product->getPdgtype()))
			 << " vertex: (" << iter_vertex->getOrigin().getVx() << ", " << iter_vertex->getOrigin().getVy() << ", " << iter_vertex->getOrigin().getVz() << ")"
			 << " Ebeam: " << iter_vertex->getOrigin().getEbeam()
			 << endl;

	  p4photon_init.SetXYZT(0,0,iter_vertex->getOrigin().getEbeam(),iter_vertex->getOrigin().getEbeam());

	  // Find first K+
	  if(nKp==0 && iter_product->getPdgtype() == PDGtype(KPlus)){
	    nKp++;
	    p4kp1.SetXYZT(iter_product->getMomentum().getPx(),iter_product->getMomentum().getPy(),iter_product->getMomentum().getPz(),iter_product->getMomentum().getE());
	    if(debug)
	      cout << "K+ 1  : (" << p4kp1.X() << ", " << p4kp1.Y() << ", " << p4kp1.Z() << ", " << p4kp1.M() << ")" << endl;

	    // Need this break since without it we would count the
	    // the first K+ as the second K+ too
	    continue;
	  }

	  // Find second K+
	  if(nKp==1 && iter_product->getPdgtype() == PDGtype(KPlus)){
	    nKp++;
	    p4kp1.SetXYZT(iter_product->getMomentum().getPx(),iter_product->getMomentum().getPy(),iter_product->getMomentum().getPz(),iter_product->getMomentum().getE());
	    if(debug)
	      cout << "K+ 2  : (" << p4kp2.X() << ", " << p4kp2.Y() << ", " << p4kp2.Z() << ", " << p4kp2.M() << ")" << endl;
	  }

	  // Find Xi-
	  if(iter_product->getPdgtype() == PDGtype(XiMinus)){
	    nXiMinus++;
	    idXiMinus = iter_product->getId();
	    p4XiMinus.SetXYZT(iter_product->getMomentum().getPx(),iter_product->getMomentum().getPy(),iter_product->getMomentum().getPz(),iter_product->getMomentum().getE());
	  }

	  // Find pi- from Xi-
	  if(iter_product->getPdgtype() == PDGtype(PiMinus) && iter_product->getParentId() == idXiMinus){
	    nPimFromXiMinus++;
	    p4pimFromXiMinus.SetXYZT(iter_product->getMomentum().getPx(),iter_product->getMomentum().getPy(),iter_product->getMomentum().getPz(),iter_product->getMomentum().getE());
	  }

	  // Find Lambda
	  if(iter_product->getPdgtype() == PDGtype(Lambda) && iter_product->getParentId() == idXiMinus){
	    nLambda++;
	    idLambda = iter_product->getId();
	    p4Lambda.SetXYZT(iter_product->getMomentum().getPx(),iter_product->getMomentum().getPy(),iter_product->getMomentum().getPz(),iter_product->getMomentum().getE());
	  }
	  
	  // Find proton
	  if(iter_product->getPdgtype() == PDGtype(Proton) && iter_product->getParentId() == idLambda){
	    nProton++;
	    p4proton.SetXYZT(iter_product->getMomentum().getPx(),iter_product->getMomentum().getPy(),iter_product->getMomentum().getPz(),iter_product->getMomentum().getE());
	  }

	  // Find pi- from Lambda
	  if(iter_product->getPdgtype() == PDGtype(PiMinus) && iter_product->getParentId() == idLambda){
	    nPimFromLambda++;
	    p4pimFromLambda.SetXYZT(iter_product->getMomentum().getPx(),iter_product->getMomentum().getPy(),iter_product->getMomentum().getPz(),iter_product->getMomentum().getE());
	  }
	} // end of loop over ProductList
      }	// end of loop over VertexList
    } // end of loop over ReactiionList

    // If we n a proton decaying from Lambda,
    // save this event
    if(nKp==2 && nXiMinus==1 && nLambda==1 && nProton==1 && nPimFromXiMinus==1 && nPimFromLambda==1){

      // make sure that total 4-mom of p, pi+, pi- are close to
      // initial state
      p4diff = p4photon_init + p4proton_init - p4kp1 - p4kp2 - p4proton - p4pimFromXiMinus - p4pimFromLambda;
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
	return true;
      }else{
	cout << "Added 4-mom of final state particles does not match generated by more than 0.1 MeV!!!!" << endl;
	cout << "diff in |p| = " << p4diff.P()* 1000. << " MeV/c" << endl;
	cout << "diff  = (" << p4diff.X() << ", " << p4diff.Y() << ", " << p4diff.Z() << ", " << p4diff.E() << ")" << endl;
	cout << "total = (" << p4total.X() << ", " << p4total.Y() << ", " << p4total.Z() << ", " << p4total.E() << ")" << endl;
	abort();
      }
    } // found proton, pi+, pi-

  } // end of select_type 6
  //__________________________________________________________________________________________________

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

    hddm_r::ReconstructedPhysicsEvent &re = record.getReconstructedPhysicsEvent();
    hddm_r::ReactionList reactions = re.getReactions();
    // Loop over ReactionList
    hddm_r::ReactionList::iterator iter;
    for (iter = reactions.begin(); iter != reactions.end(); ++iter) {
      int ireaction = std::distance(reactions.begin(), iter);
      
      hddm_r::VertexList vertices = iter->getVertices();
      hddm_r::VertexList::iterator iter_vertex;
      // Loop over VertexList
      for (iter_vertex = vertices.begin(); iter_vertex != vertices.end(); ++iter_vertex){
	int ivertex = std::distance(vertices.begin(), iter_vertex);

	hddm_r::ProductList products = iter_vertex->getProducts();
	hddm_r::ProductList::iterator iter_product;
	// Loop over ProductList
	for (iter_product = products.begin(); iter_product != products.end(); ++iter_product){
	  int iparticle = std::distance(products.begin(), iter_product);

	  if(debug) cout << "reaction: " << ireaction << " vertex: " << ivertex << " particle: " << iparticle
			 << " id: " << setw(3) << iter_product->getId() << " parent: " << iter_product->getParentId()
			 << " PDG type: " << setw(7) << iter_product->getPdgtype()
			 << " type: " << ParticleType((Particle_t)PDGtoPType((Particle_t)iter_product->getPdgtype()))
			 << " vertex: (" << iter_vertex->getOrigin().getVx() << ", " << iter_vertex->getOrigin().getVy() << ", " << iter_vertex->getOrigin().getVz() << ")"
			 << " Ebeam: " << iter_vertex->getOrigin().getEbeam()
			 << endl;

	  p4photon_init.SetXYZT(0,0,iter_vertex->getOrigin().getEbeam(),iter_vertex->getOrigin().getEbeam());
	  
	  // Find proton
	  if(iter_product->getPdgtype() == PDGtype(Proton)){
	    foundProton = true;
	    idProton = iter_product->getId();
	    p4proton.SetXYZT(iter_product->getMomentum().getPx(),iter_product->getMomentum().getPy(),iter_product->getMomentum().getPz(),iter_product->getMomentum().getE());
	  }

	  // Find pi+
	  if(iter_product->getPdgtype() == PDGtype(PiPlus)){
	    foundPip = true;
	    idPip = iter_product->getId();
	    p4pip.SetXYZT(iter_product->getMomentum().getPx(),iter_product->getMomentum().getPy(),iter_product->getMomentum().getPz(),iter_product->getMomentum().getE());
	  }

	  // Find pi-
	  if(iter_product->getPdgtype() == PDGtype(PiMinus)){
	    foundPim = true;
	    idPim = iter_product->getId();
	    p4pim.SetXYZT(iter_product->getMomentum().getPx(),iter_product->getMomentum().getPy(),iter_product->getMomentum().getPz(),iter_product->getMomentum().getE());
	  }

	  // Find pi-
	  if(iter_product->getPdgtype() == PDGtype(Pi0)){
	    foundPi0 = true;
	    idPi0 = iter_product->getId();
	    p4pi0.SetXYZT(iter_product->getMomentum().getPx(),iter_product->getMomentum().getPy(),iter_product->getMomentum().getPz(),iter_product->getMomentum().getE());
	  }

	} // end of loop over ProductList
      }	// end of loop over VertexList
    } // end of loop over ReactiionList

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
