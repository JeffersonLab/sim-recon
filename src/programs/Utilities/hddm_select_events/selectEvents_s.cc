// $Id$
//
// Created Oct 10, 2013  David Lawrence

#include "hddm_select_events.h"
#include <HDDM/hddm_s.h>
#include "particleType.h"

using namespace std;

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

  // default is to return false
  return false;
}
