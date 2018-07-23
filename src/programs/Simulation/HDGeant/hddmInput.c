/*
 * hddmInput  - functions to handle Monte Carlo generator input to HDGeant
 *              through the standard hddm i/o mechanism.
 *
 * Interface:
 *      openInput(filename) - open input stream to file <filename>
 *      skipInput(count) - skip next <count> events on open input file
 *      nextInput() - advance to next event on open input stream
 *      loadInput() - push current input event to Geant kine structures
 *      storeInput() - pop current input event from Geant kine structures
 *      closeInput() - close currently open input stream
 *
 * Richard Jones
 * University of Connecticut
 * July 13, 2001
 *
 * Usage Notes:
 * 1) Most Monte Carlo generators do not care where the vertex is placed
 *    inside the target, and specify only the final-state particles'
 *    momenta.  In this case the vertex position has to be randomized by
 *    the simulation within the beam/target overlap volume.  If the vertex
 *    position from the generator is (0,0,0) then the simulation vertex is
 *    generated uniformly inside the cylinder specified by TARGET_LENGTH,
 *    BEAM_DIAMETER, and TARGET_CENTER defined below.
 * 2) The start time for the event in HDGeant is defined to be the
 *    instant the beam photon passes through the midplane of the target,
 *    or would have passed through the midplane if it had gotten that far.
 *
 * Revision history:
 *
 * > Aug 17, 2007 - David Lawrence
 * Fill in id, parentid, pdgtype, and mech fields of reactions objects in HDDM.
 * Mostly zeros, but it makes it clear the fields are invalid and allows
 * cleaner printing.
 *
 * > Nov 17, 2006 - Richard Jones
 *      Added code to load_event that sets the Geant tofg parameter so
 *      that the start time of the event conforms to note (2) above.
 * 
 * > Apr 10, 2006 - David Lawrence
 *      Added comments to explain a little what each of these routines is 
 *      doing.  No functional changes.
 *
 * >  Dec 15, 2004 - Richard Jones
 *      Changed former behaviour of simulation to overwrite the vertex
 *      coordinates from the input record, if the simulation decides to
 *      override the input values.  At present this happens whenever the
 *      input record specifies 0,0,0 for the vertex, but in the future it
 *      may be decided to let the simulator determine the vertex position
 *      in other cases.  Since it is not part of the simulation proper, the
 *      decision was made to store this information in the reaction tag.
 */

#define TARGET_LENGTH 29.9746
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

float beam_momentum[4];
float target_momentum[4];

float get_beam_momentum_(const int *comp) {
   return beam_momentum[*comp];
}

float get_target_momentum_(const int *comp) {
   return target_momentum[*comp];
}

int extractRunNumber(int *runNo) {
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
      return 9;     /* input stream was never opened */
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
int loadInput (int override_run_number, int myInputRunNo)
{
   /* Extracts the "thrown" particle 4-vectors and types from the
    * current HDDM buffer "thisInputEvent" and creates a vertex for
    * them (gsvert) and defines the GEANT (gskine) for tracking.
    */
   s_Reactions_t* reacts;
   int reactCount, ir;
   int runNo = (override_run_number>0)?(myInputRunNo):thisInputEvent->physicsEvents->in[0].runNo;
   int eventNo = thisInputEvent->physicsEvents->in[0].eventNo;
   seteventid_(&runNo,&eventNo);
   reacts = thisInputEvent ->physicsEvents->in[0].reactions;
   if (reacts == 0)
      return 1;
   reactCount = reacts->mult;
   for (ir = 0; ir < reactCount; ir++)
   {
      s_Reaction_t* react = &reacts->in[ir];
      s_Beam_t* beam = react->beam;
      s_Target_t* target = react->target;
      s_Vertices_t* verts = react->vertices;
      int vertCount = verts->mult;
      s_Vertex_t* vert;
      float zero = 0;
      float v0[4];
      int iv;

      if (vertCount == 0)
      {
         continue;
      }

      vert = &verts->in[0];
      v0[0] = vert->origin->vx;
      v0[1] = vert->origin->vy;
      v0[2] = vert->origin->vz;
      v0[3] = vert->origin->t;
      if ((v0[0] == 0) && (v0[1] == 0) && (v0[2] == 0))
      {
         v0[0] = 1;
         v0[1] = 1;
         v0[2] = TARGET_CENTER;
         while (v0[0]*v0[0] + v0[1]*v0[1] > 0.25)
         {
            int len = 3;
            grndm_(v0,&len);
            v0[0] -= 0.5;
            v0[1] -= 0.5;
            v0[2] -= 0.5;
         }
         v0[0] *= BEAM_DIAMETER;
         v0[1] *= BEAM_DIAMETER;
         v0[2] *= TARGET_LENGTH;
         v0[2] += TARGET_CENTER;
         v0[3] = (v0[3] == 0)? settofg_(v0,&zero) * 1e9 : 0;
      }
      else
      {
         v0[3] = (v0[3] == 0)? settofg_(v0,&zero) * 1e9 : 0;
         v0[0] = 0;
         v0[1] = 0;
         v0[2] = 0;
      }

      if (beam != NULL && beam != (s_Beam_t*)&hddm_s_nullTarget)
      {
         beam_momentum[0] = beam->momentum->E;
         beam_momentum[1] = beam->momentum->px;
         beam_momentum[2] = beam->momentum->py;
         beam_momentum[3] = beam->momentum->pz;
      }
      else
      {
         beam_momentum[0] = 0;
         beam_momentum[1] = 0;
         beam_momentum[2] = 0;
         beam_momentum[3] = 0;
      }

      if (target != NULL && target != (s_Target_t*)&hddm_s_nullTarget)
      {
         target_momentum[0] = target->momentum->E;
         target_momentum[1] = target->momentum->px;
         target_momentum[2] = target->momentum->py;
         target_momentum[3] = target->momentum->pz;
      }
      else
      {
         target_momentum[0] = 0;
         target_momentum[1] = 0;
         target_momentum[2] = 0;
         target_momentum[3] = 0;
      }

      for (iv = 0; iv < vertCount; iv++)
      {
         int ntbeam = 0;
         int nttarg = 0;
         int nubuf = 0;
         float ubuf;
         int nvtx;
         float v[4];
         s_Products_t* prods;
         int prodCount, ip;
         vert = &verts->in[iv];
         v[0] = vert->origin->vx += v0[0];
         v[1] = vert->origin->vy += v0[1];
         v[2] = vert->origin->vz += v0[2];
         v[3] = vert->origin->t += v0[3];
         settofg_(v, &v[3]);
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
      ps->in[ps->mult].pdgtype = 0; /* don't bother with the PDG type here */
      ps->in[ps->mult].id = itra;   /* unique value for this particle within the event */
      gidSet(itra, itra);           /* assume same value for geant id */
      ps->in[ps->mult].parentid = 0;/* All internally generated particles have no parent */
      ps->in[ps->mult].mech = 0;    /* maybe this should be set to something? */
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
 * storeBeam
 *-------------------------
 */
int storeBeam (float vect[7], float t0)
{
   /* This is called from gukine in the case where the user wants to
    * halt simulation and save the present (single) track in the
    * Monte Carlo reactions header, perhaps for simulation later.
    * The original vertex information is moved into the beam tag
    * and then overwritten with the state of the current track.
    * This function assumes that storeInput has already been called
    * at least once for this event.
    */

   s_PhysicsEvents_t* pes;
   s_Reactions_t* rs;
   s_Beam_t* bs;
   s_Vertices_t* vs;
   s_Origin_t* or;
   s_Products_t* ps;
   int nvtx, ntbeam, nttarg, itra, nubuf;
   float vert[3], plab[3], tofg, ubuf[10];
   Particle_t kind;
   int ilast;

   pes = thisInputEvent->physicsEvents;
   if (pes == 0 || pes == HDDM_NULL || pes->mult == 0)
      return 0;
   rs = pes->in[0].reactions;
   if (rs == 0 || rs == HDDM_NULL || rs->mult == 0)
      return 0;
   vs = rs->in[0].vertices;
   if (vs == 0 || vs == HDDM_NULL || vs->mult == 0)
      return 0;
   ps = vs->in[0].products;
   if (ps == 0 || ps == HDDM_NULL || ps->mult == 0)
      return 0;
   bs = rs->in[0].beam;
   if (bs == HDDM_NULL)
      bs = make_s_Beam();
   ilast = ps->mult - 1;
   bs->type = ps->in[ilast].type;
   if (bs->momentum != HDDM_NULL)
      FREE(bs->momentum);
   bs->momentum = ps->in[ilast].momentum;
   ps->in[ilast].momentum = HDDM_NULL;
   if (bs->polarization != HDDM_NULL)
      FREE(bs->polarization);
   bs->polarization = ps->in[ilast].polarization;
   ps->in[ilast].polarization = HDDM_NULL;
   if (bs->properties != HDDM_NULL)
      FREE(bs->properties);
   bs->properties = ps->in[ilast].properties;
   ps->in[ilast].properties = HDDM_NULL;
   rs->in[0].beam = bs;
   for (itra = 1; itra <= 1; itra++) {
      char chnpar[99];
      int itrtyp;
      float amass,charge,tlife;
      gfkine_(&itra,vert,plab,&kind,&nvtx,ubuf,&nubuf);
      gfpart_(&kind,chnpar,&itrtyp,&amass,&charge,&tlife,ubuf,&nubuf);
      gfvert_(&nvtx,vert,&ntbeam,&nttarg,&tofg,ubuf,&nubuf);
      or = vs->in[0].origin;
      if (or == HDDM_NULL) {
         vs->in[0].origin = or = make_s_Origin();
      }
      or->vx = vect[0];
      or->vy = vect[1];
      or->vz = vect[2];
      or->t = t0 * 1e9;
      ps->in[ilast].type = kind;
      ps->in[ilast].pdgtype = 22;   /* assume a beam photon */
      ps->in[ilast].id = itra;
      ps->in[ilast].parentid = 0;
      ps->in[ilast].mech = 0;
      ps->in[ilast].momentum = make_s_Momentum();
      ps->in[ilast].momentum->px = vect[6] * vect[3];
      ps->in[ilast].momentum->py = vect[6] * vect[4];
      ps->in[ilast].momentum->pz = vect[6] * vect[5];
      ps->in[ilast].momentum->E  = vect[6];
   }
   return 1;
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
   if (thisInputEvent == NULL)
      return 0;
   if (thisInputEvent->physicsEvents == NULL)
      return 0;
   if (thisInputEvent->physicsEvents->mult < 1)
      return 0;
   s_PhysicsEvent_t *pe = &thisInputEvent->physicsEvents->in[0];
   if (pe->reactions == NULL)
      return 0;
   if (pe->reactions->mult < 1)
      return 0;
   s_Random_t *rnd = pe->reactions->in[0].random;
   if (rnd == NULL || rnd == HDDM_NULL) {
      /* No seeds stored in event. Return */
      return 0;
   }
   else {
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
   if (thisInputEvent == NULL)
      return 0;
   if (thisInputEvent->physicsEvents == NULL)
      return 0;
   if (thisInputEvent->physicsEvents->mult < 1)
      return 0;
   s_PhysicsEvent_t *pe = &thisInputEvent->physicsEvents->in[0];
   if (pe->reactions == NULL)
      return 0;
   if (pe->reactions->mult < 1)
      return 0;
   s_Random_t *rnd = pe->reactions->in[0].random;
   if (rnd == NULL || rnd == HDDM_NULL) {
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

int loadinput_ (int *override_run_number,int *myInputRunNo)
{
  return loadInput(*override_run_number,*myInputRunNo);
}

int storeinput_ (int* runNo, int* eventNo, int* ntracks)
{
   return storeInput(*runNo,*eventNo,*ntracks);
}

int storebeam_ (float* vect, float* t)
{
   float t0 = *t;
   return storeBeam(vect, t0);
}

int closeinput_ ()
{
   return closeInput();
}

int extractrunnumber_(int *runNo){
  return extractRunNumber(runNo);
}
