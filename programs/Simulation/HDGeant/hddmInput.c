/*
 * hddmInput  - functions to handle Monte Carlo generator input to HDGeant
 *		through the standard hddm i/o mechanism.
 *
 * Interface:
 *	openInput(filename) - open input stream to file <filename>
 *	skipInput(count) - skip next <count> events on open input file
 *      nextInput() - advance to next event on open input stream
 *      loadInput() - push current input event to Geant kine structures
 *	getInput() - return pointer to top structure of current event
 *	closeInput() - close currently open input stream
 *
 * Richard Jones
 * University of Connecticut
 * July 13, 2001
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <hddm_s.h>

s_iostream* thisInputStream = 0;
s_HDDM_t* thisInputEvent = 0;

int openInput (char* filename)
{
   thisInputStream = open_s_HDDM(filename);
   return (thisInputStream == 0);
}

int skipInput (int count)
{
   int* buff = malloc(1000000);
   while (thisInputStream && (count > 0))
   {
      int ret = fread(buff, sizeof(int), 1, thisInputStream->fd);
      if (ret == 1)
      {
         int nw = *buff;
         ret = fread(buff, sizeof(int), nw, thisInputStream->fd);
         --count;
      }
      else
      {
         break;
      }
   }
   free(buff);
   return count;
}

int nextInput ()
{
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

int loadInput ()
{
   s_Reactions_t* reacts;
   int reactCount, ir;

   reacts = thisInputEvent ->physicsEvent ->reactions;
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
         int ntbeam = 0;
         int nttarg = 0;
         int nubuf = 0;
         float ubuf;
         int nvtx;
         s_Products_t* prods;
         int prodCount, ip;
         s_Vertex_t* vert = &verts->in[iv];

#define TARGET_LENGTH 30
#define BEAM_DIAMETER 0.5
#define TARGET_CENTER 65

#ifdef USE_MC_VERTEX
         v[0] = vert->origin->vx;
         v[1] = vert->origin->vy;
         v[2] = vert->origin->vz;
#else
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
#endif
 
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
            p[0] = prod->momentum->px;
            p[1] = prod->momentum->py;
            p[2] = prod->momentum->pz;
            if (prod->decayVertex == 0)
            {
               gskine_(p, &kind, &nvtx, &ubuf, &nubuf, &ntrk);
            }
         }
      }
   }
}

int closeInput ()
{
   if (thisInputStream)
   {
      close_s_HDDM(thisInputStream);
      thisInputStream = 0;
   }
   return 0;
}

s_HDDM_t* getInput()
{
   return thisInputEvent;
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

int closeinput_ ()
{
   return closeInput();
}
