/*
savenewvertex: particle stoped because it decayed
               save the daughter particles and the
               decay vertex.


*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <hddm_s.h>
#include <geant3.h>

extern s_HDDM_t* thisInputEvent;

/*
  The following are used to communicate the "id" fields from the
  particles read in from the input file of generated events here
  so that new particles created via decays can have the proper
  parent id recorded and get assigned their own id that is not
  already assigned.
  
  9/11/2013
*/
extern int MAX_GENERATED_PARTICLES;
extern int MYID[];
extern int NEXT_MYID;


void SaveNewVertex(int kcase, int Npart, float *gkin, 
		   float vertex[3], float tofg, int *iflgk,
                   int ipart, int itra, int istak) {


  // get pointer to all vertices
  s_Vertices_t* verts = thisInputEvent->physicsEvents->in[0].reactions->in[0].vertices;
  int VertexCount = verts->mult;

  // create additional space for one more vertex
  s_Vertices_t* NewVerts = make_s_Vertices(VertexCount+1);
  int i;
  for (i=0;i<VertexCount;i++){
    NewVerts->in[i] = verts->in[i];
  }
  NewVerts->mult = VertexCount;
  thisInputEvent->physicsEvents->in[0].reactions->in[0].vertices = NewVerts;
  free(verts);
  verts = NewVerts;
  //verts->in[VertexCount] = make_s_Vertex();
  verts->mult++;
  //printf("Number of stored Vertices is now: %d\n",verts->mult);
  
  // copy in the new vertex coordinates
  s_Origin_t* or = make_s_Origin();
  verts->in[VertexCount].origin = or;
  or->vx = vertex[0];
  or->vy = vertex[1];
  or->vz = vertex[2];
  or->t = tofg;

  // copy in the new particles at this vertex
  s_Products_t* ps = make_s_Products(Npart);
  verts->in[VertexCount].products = ps;
  ps->mult = Npart;
  for (i = 0;i<Npart;i++){
  	
	/* Give this particle a unique id in the output file */
	int myid = NEXT_MYID++;
	
	/* Find parent id based on parent's track number */
	int parentid = itra < MAX_GENERATED_PARTICLES ? MYID[itra]:0;

	/* Add this particle to the map of track number to id in case it also is decayed */
	int my_trk_num = iflgk[i];
	if(my_trk_num < MAX_GENERATED_PARTICLES) MYID[my_trk_num] = myid;
    
    ps->in[i].momentum = make_s_Momentum();
    ps->in[i].momentum->px = gkin[i*5+0];
    ps->in[i].momentum->py = gkin[i*5+1];
    ps->in[i].momentum->pz = gkin[i*5+2];
    ps->in[i].momentum->E  = gkin[i*5+3];
    ps->in[i].type = gkin[i*5+4];
    ps->in[i].parentid = parentid;
    ps->in[i].id = myid;
    ps->in[i].mech = kcase;
    ps->in[i].decayVertex = VertexCount;

  }
    
}



/* entry point from fortran */

void savenewvertex_ (int *kcase, int *N, float* gkin, 
		     float* vertex, float* tofg, int* iflgk, 
                     int* ipart, int* itra, int* istak) {

  SaveNewVertex(*kcase, *N, gkin, vertex, *tofg, iflgk, *ipart, *itra, *istak);

}
