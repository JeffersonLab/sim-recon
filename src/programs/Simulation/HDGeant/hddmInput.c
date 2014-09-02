/*
 * hddmInput  - functions to handle Monte Carlo generator input to HDGeant
 *		through the standard hddm i/o mechanism.
 *
 * Interface:
 *	openInput(filename) - open input stream to file <filename>
 *	skipInput(count) - skip next <count> events on open input file
 *      nextInput() - advance to next event on open input stream
 *      loadInput() - push current input event to Geant kine structures
 *      storeInput() - pop current input event from Geant kine structures
 *	closeInput() - close currently open input stream
 *
 * Richard Jones
 * University of Connecticut
 * July 13, 2001
 *
 * Usage Notes:
 * 1)	Most Monte Carlo generators do not care where the vertex is placed
 *	inside the target, and specify only the final-state particles'
 *	momenta.  In this case the vertex position has to be randomized by
 *	the simulation within the beam/target overlap volume.  If the vertex
 *	position from the generator is (0,0,0) then the simulation vertex is
 *	generated uniformly inside the cylinder specified by TARGET_LENGTH,
 *	BEAM_DIAMETER, and TARGET_CENTER defined below.
 * 2)   The start time for the event in HDGeant is defined to be the
 *      instant the beam photon passes through the midplane of the target,
 *      or would have passed through the midplane if it had gotten that far.
 *
 * Revision history:
 *
 * > Aug 17, 2007 - David Lawrence
 * Fill in id, parentid, pdgtype, and mech fields of reactions objects in HDDM.
 * Mostly zeros, but it makes it clear the fields are invalid and allows
 * cleaner printing.
 *
 * > Nov 17, 2006 - Richard Jones
 * 	Added code to load_event that sets the Geant tofg parameter so
 *      that the start time of the event conforms to note (2) above.
 * 
 * > Apr 10, 2006 - David Lawrence
 *	Added comments to explain a little what each of these routines is 
 *	doing.  No functional changes.
 *
 * >  Dec 15, 2004 - Richard Jones
 *	Changed former behaviour of simulation to overwrite the vertex
 *	coordinates from the input record, if the simulation decides to
 *	override the input values.  At present this happens whenever the
 *      input record specifies 0,0,0 for the vertex, but in the future it
 * 	may be decided to let the simulator determine the vertex position
 *	in other cases.  Since it is not part of the simulation proper, the
 *	decision was made to store this information in the reaction tag.
 */

#define TARGET_LENGTH 30
#define BEAM_DIAMETER 0.5
#define TARGET_CENTER 65


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <HDDM/hddm_s.h>
#include <geant3.h>

#include "gid_map.h"

void seteventid_(int *runNo, int *eventNo);

float settofg_(float origin[3], float *time0);

s_iostream_t* thisInputStream = 0;
s_HDDM_t* thisInputEvent = 0;

int extractRunNumber(int *runNo){
  thisInputEvent = read_s_HDDM(thisInputStream);
  return *runNo = thisInputEvent->physicsEvents->in[0].runNo;
}



/*-------------------------
 * openInput
 *-------------------------
 */
int openInput (char* filename)
{
	/* Open HDDM file for reading in "thrown" particle kinematics */
   thisInputStream = open_s_HDDM(filename);
   return (thisInputStream == 0);
}

/*-------------------------
 * skipInput
 *-------------------------
 */
int skipInput (int count)
{
   return count - skip_s_HDDM(thisInputStream,count);
}

/*-------------------------
 * nextInput
 *-------------------------
 */
int nextInput ()
{
	/* Read in the next HDDM event. This only reads it into the
	 * HDDM buffer "thisInputEvent" and does not yet define the
	 * particles to GEANT. See loadInput for that.
	 */
   if (thisInputStream == 0)
   {
      return 9;		/* input stream was never opened */
   }
   else if (thisInputEvent)
   {
      flush_s_HDDM(thisInputEvent, 0);
   }
   thisInputEvent = read_s_HDDM(thisInputStream);
   return (thisInputEvent == 0);
}

/*-------------------------
 * loadInput
 *-------------------------
 */
int loadInput ()
{
	/* Extracts the "thrown" particle 4-vectors and types from the
	 * current HDDM buffer "thisInputEvent" and creates a vertex for
	 * them (gsvert) and defines the GEANT (gskine) for tracking.
	 */
   s_Reactions_t* reacts;
   int reactCount, ir;
   int runNo = thisInputEvent->physicsEvents->in[0].runNo;
   int eventNo = thisInputEvent->physicsEvents->in[0].eventNo;
   seteventid_(&runNo,&eventNo);
   reacts = thisInputEvent ->physicsEvents->in[0].reactions;
   reactCount = reacts->mult;
   for (ir = 0; ir < reactCount; ir++)
   {
      s_Vertices_t* verts;
      int vertCount, iv;
      s_Reaction_t* react = &reacts->in[ir];
      verts = react->vertices;
      vertCount = verts->mult;
      for (iv = 0; iv < vertCount; iv++)
      {
         float v[3];
         float time0;
         int ntbeam = 0;
         int nttarg = 0;
         int nubuf = 0;
         float ubuf;
         int nvtx;
         s_Products_t* prods;
         int prodCount, ip;
         s_Vertex_t* vert = &verts->in[iv];
         v[0] = vert->origin->vx;
         v[1] = vert->origin->vy;
         v[2] = vert->origin->vz;
         if ((v[0] == 0) && (v[1] == 0) && (v[2] == 0))
         {
            v[0] = 1;
            v[1] = 1;
            v[2] = TARGET_CENTER;
            while (v[0]*v[0] + v[1]*v[1] > 0.25)
            {
               int len = 3;
               grndm_(v,&len);
               v[0] -= 0.5;
               v[1] -= 0.5;
               v[2] -= 0.5;
            }
            v[0] *= BEAM_DIAMETER;
            v[1] *= BEAM_DIAMETER;
            v[2] *= TARGET_LENGTH;
            v[2] += TARGET_CENTER;
            vert->origin->vx = v[0];
            vert->origin->vy = v[1];
            vert->origin->vz = v[2];
         }
         time0 = vert->origin->t;
         vert->origin->t = settofg_(v,&time0) * 1e9;
         gsvert_(v, &ntbeam, &nttarg, &ubuf, &nubuf, &nvtx);
         prods = vert->products;
         prodCount = prods->mult;
         for (ip = 0; ip < prodCount; ip++)
         {
            int ntrk;
            float p[3];
            Particle_t kind;
            s_Product_t* prod = &prods->in[ip];
            kind = prod->type;
				
	    /* Don't tell geant to track particles that are intermediary types */
	    if (kind <= 0)
              continue;
				
            p[0] = prod->momentum->px;
            p[1] = prod->momentum->py;
            p[2] = prod->momentum->pz;
            if (prod->decayVertex == 0)
            {
               gskine_(p, &kind, &nvtx, &ubuf, &nubuf, &ntrk);
	       gidSet(ntrk, ip + 1);
            }
         }
      }
   }
   return 0;
}

/*-------------------------
 * storeInput
 *-------------------------
 */
int storeInput (int runNo, int eventNo, int ntracks)
{
	/* This is called by the built-in generators (coherent brem. and
	 * single track) in order to store the "thrown" particle parameters
	 * in the output HDDM file. What this actually does is free the 
	 * input buffer "thisInputEvent" if it exists and creates a new
	 * one. When an external generator is used, the thisInputEvent
	 * buffer is kept unmodified and this routine is never called.
         */
   s_PhysicsEvents_t* pes;
   s_Reactions_t* rs;
   s_Vertices_t* vs;
   s_Origin_t* or;
   s_Products_t* ps;
   int nvtx, ntbeam, nttarg, itra, nubuf;
   float vert[3], plab[3], tofg, ubuf[10];
   Particle_t kind;

   if (thisInputEvent)
   {
      flush_s_HDDM(thisInputEvent, 0);
   }
   thisInputEvent = make_s_HDDM();
   thisInputEvent->physicsEvents = pes = make_s_PhysicsEvents(1);
   pes->in[0].reactions = rs = make_s_Reactions(1);
   pes->mult = 1;
   rs->mult = 1;
   rs->in[0].vertices = vs = make_s_Vertices(99);
   /*printf("Make 99 Vertices!!!!!!!\n");*/
   vs->mult = 0;
   for (itra = 1; itra <= ntracks; itra++)
   {
      char chnpar[99];
      int itrtyp;
      float amass,charge,tlife;
      gfkine_(&itra,vert,plab,&kind,&nvtx,ubuf,&nubuf);
      gfpart_(&kind,chnpar,&itrtyp,&amass,&charge,&tlife,ubuf,&nubuf);
      if (nvtx < 1)
      {
         return 1;
      }
      else
      {
         vs->mult = (nvtx < vs->mult)? vs->mult : nvtx;
      }
      gfvert_(&nvtx,vert,&ntbeam,&nttarg,&tofg,ubuf,&nubuf);
      or = vs->in[nvtx-1].origin;
      ps = vs->in[nvtx-1].products;
      if (or == HDDM_NULL)
      {
         or = make_s_Origin();
         vs->in[nvtx-1].origin = or;
         or->vx = vert[0];
         or->vy = vert[1];
         or->vz = vert[2];
         or->t = tofg * 1e9;
      }
      if (ps == HDDM_NULL)
      {
         ps = make_s_Products(ntracks);
         vs->in[nvtx-1].products = ps;
         ps->mult = 0;
      }
      ps->in[ps->mult].type = kind;
      ps->in[ps->mult].pdgtype = 0;	/* don't bother with the PDG type here */
      ps->in[ps->mult].id = itra;	/* unique value for this particle within the event */
      gidSet(itra, itra);	        /* assume same value for geant id */
      ps->in[ps->mult].parentid = 0;	/* All internally generated particles have no parent */
      ps->in[ps->mult].mech = 0;	/* maybe this should be set to something? */
      ps->in[ps->mult].momentum = make_s_Momentum();
      ps->in[ps->mult].momentum->px = plab[0];
      ps->in[ps->mult].momentum->py = plab[1];
      ps->in[ps->mult].momentum->pz = plab[2];
      ps->in[ps->mult].momentum->E  = sqrt(plab[0]*plab[0]+plab[1]*plab[1]
                                          +plab[2]*plab[2]+amass*amass);
      ps->mult++;
   }
   pes->in[0].runNo = runNo;
   pes->in[0].eventNo = eventNo;
   return 0;
}

/*-------------------------
 * getseeds_
 *-------------------------
 */
int getseeds_(int *iseed1, int *iseed2)
{
	/* This checks to see if thisInputStream already
	   contains random number seeds for this event. If
	   it does, then those values are copied into the
	   iseed1 and iseed2 variables. If not, then the
	   contents of iseed1 and iseed2 are left unchanged.
	*/
	if(thisInputEvent == NULL)return 0;
	if(thisInputEvent->physicsEvents == NULL)return 0;
	if(thisInputEvent->physicsEvents->mult<1)return 0;
	s_PhysicsEvent_t *pe = &thisInputEvent->physicsEvents->in[0];
	if(pe->reactions == NULL)return 0;
	if(pe->reactions->mult<1)return 0;
	s_Random_t *rnd = pe->reactions->in[0].random;
	if(rnd == NULL || rnd==HDDM_NULL){
		/* No seeds stored in event. Return */
		return 0;
	}else{
		/* Seeds found in event, copy them back to caller for use */
		*iseed1 = rnd->seed1;
		*iseed2 = rnd->seed2;
		return 1;
	}
}

/*-------------------------
 * storeseeds_
 *-------------------------
 */
int storeseeds_(int *iseed1, int *iseed2)
{
	/* This copies the given seed values into
	   thisInputStream, overwriting any values that
	   already exist there.
	*/
	if(thisInputEvent == NULL)return 0;
	if(thisInputEvent->physicsEvents == NULL)return 0;
	if(thisInputEvent->physicsEvents->mult<1)return 0;
	s_PhysicsEvent_t *pe = &thisInputEvent->physicsEvents->in[0];
	if(pe->reactions == NULL)return 0;
	if(pe->reactions->mult<1)return 0;
	s_Random_t *rnd = pe->reactions->in[0].random;
	if(rnd == NULL || rnd==HDDM_NULL){
		/* No seeds stored in event. Add them */
		rnd = pe->reactions->in[0].random = make_s_Random();
		rnd->seed3 = 709975946 + pe->eventNo;
		rnd->seed4 = 912931182 + pe->eventNo;
	}

	rnd->seed1 = *iseed1;
	rnd->seed2 = *iseed2;
	
	return 0;
}

/*-------------------------
 * closeInput
 *-------------------------
 */
int closeInput ()
{
	/* Close the HDDM input file */
   if (thisInputStream)
   {
      close_s_HDDM(thisInputStream);
      thisInputStream = 0;
   }
   return 0;
}


/* entry points from Fortran */

int openinput_ (char* filename)
{
   int retcode = openInput(strtok(filename," "));
   return retcode;
}

int skipinput_ (int* count)
{
   return skipInput(*count);
}

int nextinput_ ()
{
   return nextInput();
}

int loadinput_ ()
{
   return loadInput();
}

int storeinput_ (int* runNo, int* eventNo, int* ntracks)
{
   return storeInput(*runNo,*eventNo,*ntracks);
}

int closeinput_ ()
{
   return closeInput();
}

int extractrunnumber_(int *runNo){
  return extractRunNumber(runNo);
}
