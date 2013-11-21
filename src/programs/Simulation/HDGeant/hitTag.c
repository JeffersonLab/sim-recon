/*
 * hitTag - registers hits for the tagger focal plane counters
 *
 *	This is a part of the hits package for the
 *	HDGeant simulation program for Hall D.
 *
 *	version 1.0 	-Richard Jones November 16, 2006
 *
 * Programmer's Notes:
 * -------------------
 * 1) There is no tagger in the HDGeant simulation so no tagger hits are
 *    generated during tracking.  This hitTagger() function is called at
 *    event initialization time to register the tagged photon that is
 *    supposed to have caused the event. 
 * 2) Only microscope hits are produced in this version.
 * 3) In the simulation of physics events (external generator) with
 *    background enabled, pickTagger() produces a list of tagger hits
 *    that includes the original photon from the generator plus all
 *    of the background photons.  Note that this includes many photons
 *    that never reached the GlueX target because they were stopped
 *    at the collimator.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <HDDM/hddm_s.h>
#include <geant3.h>
#include <bintree.h>

#define TWO_HIT_RESOL   25.
#define MAX_HITS        5000
#define C_CM_PER_NS     29.9792458
#define REF_TIME_Z_CM   65.

#define TAGMS_E_RANGE	1.280
#define TAGMS_E_MIN 	8.000
#define TAGMS_CHANNELS	128
#define TAGMS_T_RANGE   500.
#define TAGMS_T_MIN     -200.

binTree_t* taggerTree = 0;
static int channelCount = 0;

/* register hits during event initialization (from gukine) */

void hitTagger (float xin[4], float xout[4],
                float pin[5], float pout[5], float dEsum,
                int track, int stack, int history)
{
   double chan = TAGMS_CHANNELS*(pin[3]-TAGMS_E_MIN)/TAGMS_E_RANGE+1;
   double E = TAGMS_E_MIN+(chan+0.5)*TAGMS_E_RANGE/TAGMS_CHANNELS;
   double t = xin[3]*1e9-(xin[2]-REF_TIME_Z_CM)/C_CM_PER_NS;
   t = (fabs(t) < 1e-3)? 0 : t;

   if (chan < 1 || chan > TAGMS_CHANNELS) return;   /* tagger energy limits */

   /* post the hit to the hits tree, mark slab as hit */

   {
      int nhit;
      s_TaggerHits_t* hits;
      int mark = chan;
      void** twig = getTwig(&taggerTree, mark);
      if (*twig == 0)
      {
         s_Tagger_t* tag = *twig = make_s_Tagger();
         s_MicroChannels_t* channels = make_s_MicroChannels(1);
         hits = make_s_TaggerHits(MAX_HITS);
         hits->mult = 0;
         channels->in[0].taggerHits = hits;
         channels->in[0].column = chan;
         channels->in[0].row = 0;
         channels->in[0].E = E;
         channels->mult = 1;
         tag->microChannels = channels;
         channelCount++;
      }
      else
      {
         s_Tagger_t* tag = *twig;
         hits = tag->microChannels->in[0].taggerHits;
      }
   
      if (hits != HDDM_NULL)
      {
         for (nhit = 0; nhit < hits->mult; nhit++)
         {
            if (fabs(hits->in[nhit].t - t) < TWO_HIT_RESOL)
            {
               break;
            }
         }
         if (nhit < hits->mult)         /* ignore second hit */
         {
         }
         else if (nhit < MAX_HITS)         /* create new hit */
         {
            hits->in[nhit].t = t;
            hits->mult++;
         }
         else
         {
            fprintf(stderr,"HDGeant error in hitTagger: ");
            fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
         }
      }
   }
}

/* entry point from fortran */

void hittagger_ (float* xin, float* xout,
                 float* pin, float* pout, float* dEsum,
                 int* track, int* stack, int* history)
{
   hitTagger(xin,xout,pin,pout,*dEsum,*track,*stack,*history);
}


/* pick and package the hits for shipping */

s_Tagger_t* pickTagger ()
{
   s_Tagger_t* box;
   s_Tagger_t* item;

   if (channelCount == 0)
   {
      return HDDM_NULL;
   }

   box = make_s_Tagger();
   box->microChannels = make_s_MicroChannels(channelCount);
   while (item = (s_Tagger_t*) pickTwig(&taggerTree))
   {
      s_MicroChannels_t* channels = item->microChannels;
      int channel;
      for (channel=0; channel < channels->mult; ++channel)
      {
         s_TaggerHits_t* hits = channels->in[channel].taggerHits;

         /* constraint t values to lie within time range */
         int i;
         int iok=0;
         for (iok=i=0; i < hits->mult; i++)
         {
            if ((hits->in[i].t >= TAGMS_T_MIN) &&
                (hits->in[i].t <= (TAGMS_T_MIN+TAGMS_T_RANGE)))
            {
               if (iok < i)
               {
                  hits->in[iok] = hits->in[i];
               }
               ++iok;
            }
         }
         if (iok)
         {
            hits->mult = iok;
            int m = box->microChannels->mult++;
            box->microChannels->in[m] = channels->in[0];
         }
         else if (hits != HDDM_NULL)
         {
            FREE(hits);
         }
      }
      if (channels != HDDM_NULL)
      {
         FREE(channels);
      }
      FREE(item);
   }

   channelCount = 0;

   if ((box->microChannels != HDDM_NULL) &&
       (box->microChannels->mult == 0))
   {
      FREE(box->microChannels);
      box->microChannels = HDDM_NULL;
   }
   if (box->microChannels->mult == 0)
   {
      FREE(box);
      box = HDDM_NULL;
   }
   return box;
}
