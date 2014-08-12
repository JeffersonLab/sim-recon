/*
 * hitTag - registers hits for the tagger focal plane counters
 *
 *	This is a part of the hits package for the
 *	HDGeant simulation program for Hall D.
 *
 * version 1.0 	-Richard Jones November 16, 2006
 * version 2.0 	-Richard Jones July 1, 2014
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
 *
 * update July 1, 2014 (version 2.0)
 * ---------------------------------
 * 1) Read the tagger channel energy bounds from the ccdb instead of
 *    hard-wiring them here.
 * 2) Add hits in both the fixed_array and microscope detectors.
 * 3) Fix the bug that forced the E value written into the hits structure
 *    to always contain the exact simulated beam photon energy, instead
 *    of the mean value for the hit tagger channel. Now only the mean
 *    photon energy for the hit channel is recorded.
 * 4) The recorded photon energy from the tagger is computed from the
 *    endpoint energy in the ccdb multiplied by the scaled_energy_range
 *    array values.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <HDDM/hddm_s.h>
#include <geant3.h>
#include <bintree.h>
#include <calibDB.h>

#define BEAM_BUCKET_SPACING_NS  2.
#define MICRO_TWO_HIT_RESOL     25.
#define MICRO_MAX_HITS          5000
#define FIXED_TWO_HIT_RESOL     25.
#define FIXED_MAX_HITS          5000
#define C_CM_PER_NS             29.9792458
#define REF_TIME_Z_CM           65.
#define TAG_T_MIN_NS            -20
#define TAG_T_MAX_NS            +20

float endpoint_energy_GeV = 0;
static int micro_nchannels = 102;
float* micro_channel_Elow = 0;
float* micro_channel_Ehigh = 0;
static int fixed_nchannels = 274;
float* fixed_channel_Elow = 0;
float* fixed_channel_Ehigh = 0;

binTree_t* microTree = 0;
binTree_t* fixedTree = 0;
static int microCount = 0;
static int fixedCount = 0;
static int loadDone = 0;

/* register hits during event initialization (from gukine) */

void hitTagger (float xin[4], float xout[4],
                float pin[5], float pout[5], float dEsum,
                int track, int stack, int history)
{
   int micro_chan;
   int fixed_chan;
   double E = pin[3];
   double t = xin[3]*1e9-(xin[2]-REF_TIME_Z_CM)/C_CM_PER_NS;
   t = floor(t/BEAM_BUCKET_SPACING_NS+0.5)*BEAM_BUCKET_SPACING_NS;

   if (loadDone==0){
     /* read tagger set endpoint energy from calibdb */
     {
       char dbname[] = "/PHOTON_BEAM/endpoint_energy::mc";
       unsigned int ndata = 1;
       if (GetCalib(dbname, &ndata, &endpoint_energy_GeV)) {
	 fprintf(stderr,"HDGeant error in hitTagger: %s %s\n",
		 "failed to read photon beam endpoint energy",
		 "from calibdb, cannot continue.");
	 exit(-2);
       }
     }

     /* read microscope channel energy bounds from calibdb */
     {
       char dbname[] = "/PHOTON_BEAM/microscope/scaled_energy_range::mc";
       int ndata = micro_nchannels;
       micro_channel_Elow = malloc(ndata*sizeof(float));
       micro_channel_Ehigh = malloc(ndata*sizeof(float));
 
       GetColumn(dbname, &ndata, micro_channel_Elow,"xlow");
       GetColumn(dbname, &ndata, micro_channel_Ehigh,"xhigh");
       
       int i;
       for (i=0; i < ndata; ++i) {
	 micro_channel_Elow[i] *= endpoint_energy_GeV;
	 micro_channel_Ehigh[i] *= endpoint_energy_GeV;
       }
     }
     
      /* read fixed array channel energy bounds from calibdb */
     {
       char dbname[] = "/PHOTON_BEAM/hodoscope/scaled_energy_range::mc";
       int ndata = fixed_nchannels;
       fixed_channel_Elow = malloc(ndata*sizeof(float));
       fixed_channel_Ehigh = malloc(ndata*sizeof(float));

       GetColumn(dbname, &ndata, fixed_channel_Elow,"xlow");
       GetColumn(dbname, &ndata, fixed_channel_Ehigh,"xhigh");

       int i;
       for (i=0; i < ndata; ++i) {
	 fixed_channel_Elow[i] *= endpoint_energy_GeV;
	 fixed_channel_Ehigh[i] *= endpoint_energy_GeV;
	 
       }
     }

     fprintf(stderr,"TAGGER: ALL parameters loaded from Data Base\n");
     loadDone = 1;
   }

   /* look up hit tagger channel, if any */
   micro_chan = -1;
   if (E < endpoint_energy_GeV) {
      int i;
      for (i=0; i < micro_nchannels; ++i) {
         if ( E < micro_channel_Ehigh[i] &&
              E > micro_channel_Elow[i] )
         {
            E = (micro_channel_Ehigh[i] + micro_channel_Elow[i])/2;
            micro_chan = i;
            break;
         }
      }
   }
   fixed_chan = -1;
   if (micro_chan == -1) {
      int i;
      for (i=0; i < fixed_nchannels; ++i) {
         if ( E < fixed_channel_Ehigh[i] &&
              E > fixed_channel_Elow[i] )
         {
            E = (fixed_channel_Elow[i] + fixed_channel_Ehigh[i])/2;
            fixed_chan = i;
            break;
         }
      }
   }

   /* post the hit to the microscope hits tree, mark channel as hit */

   if (micro_chan > -1) {
      int nhit;
      s_TaggerHits_t* hits;
      int mark = micro_chan + 1000;
      void** twig = getTwig(&microTree, mark);
      if (*twig == 0)
      {
         s_Tagger_t* tag = *twig = make_s_Tagger();
         s_MicroChannels_t* channels = make_s_MicroChannels(1);
         hits = make_s_TaggerHits(MICRO_MAX_HITS);
         hits->mult = 0;
         channels->in[0].taggerHits = hits;
         channels->in[0].column = micro_chan;
         channels->in[0].row = 0;
         channels->in[0].E = E;
         channels->mult = 1;
         tag->microChannels = channels;
         microCount++;
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
            if (fabs(hits->in[nhit].t - t) < MICRO_TWO_HIT_RESOL)
            {
               break;
            }
         }
         if (nhit < hits->mult)         /* ignore second hit */
         {
         }
         else if (nhit < MICRO_MAX_HITS)   /* create new hit */
         {
            hits->in[nhit].t = t;
            hits->mult++;
         }
         else
         {
            fprintf(stderr,"HDGeant error in hitTagger: ");
            fprintf(stderr,"max hit count %d exceeded, truncating!\n",
                    MICRO_MAX_HITS);
         }
      }
   }

   /* post the hit to the fixed array hits tree, mark channel as hit */

   if (fixed_chan > -1) {
      int nhit;
      s_TaggerHits_t* hits;
      int mark = fixed_chan + 1000;
      void** twig = getTwig(&fixedTree, mark);
      if (*twig == 0)
      {
         s_Tagger_t* tag = *twig = make_s_Tagger();
         s_FixedChannels_t* channels = make_s_FixedChannels(1);
         hits = make_s_TaggerHits(FIXED_MAX_HITS);
         hits->mult = 0;
         channels->in[0].taggerHits = hits;
         channels->in[0].channel = fixed_chan;
         channels->in[0].E = E;
         channels->mult = 1;
         tag->fixedChannels = channels;
         fixedCount++;
      }
      else
      {
         s_Tagger_t* tag = *twig;
         hits = tag->fixedChannels->in[0].taggerHits;
      }
   
      if (hits != HDDM_NULL)
      {
         for (nhit = 0; nhit < hits->mult; nhit++)
         {
            if (fabs(hits->in[nhit].t - t) < FIXED_TWO_HIT_RESOL)
            {
               break;
            }
         }
         if (nhit < hits->mult)         /* ignore second hit */
         {
         }
         else if (nhit < FIXED_MAX_HITS)   /* create new hit */
         {
            hits->in[nhit].t = t;
            hits->mult++;
         }
         else
         {
            fprintf(stderr,"HDGeant error in hitTagger: ");
            fprintf(stderr,"max hit count %d exceeded, truncating!\n",
                    FIXED_MAX_HITS);
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

   if (microCount == 0 && fixedCount == 0)
   {
      return HDDM_NULL;
   }

   box = make_s_Tagger();

   box->microChannels = make_s_MicroChannels(microCount);
   while ((item = (s_Tagger_t*) pickTwig(&microTree)))
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
            if ((hits->in[i].t >= TAG_T_MIN_NS) &&
                (hits->in[i].t <= TAG_T_MAX_NS))
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

   box->fixedChannels = make_s_FixedChannels(fixedCount);
   while ((item = (s_Tagger_t*) pickTwig(&fixedTree)))
   {
      s_FixedChannels_t* channels = item->fixedChannels;
      int channel;
      for (channel=0; channel < channels->mult; ++channel)
      {
         s_TaggerHits_t* hits = channels->in[channel].taggerHits;

         /* constraint t values to lie within time range */
         int i;
         int iok=0;
         for (iok=i=0; i < hits->mult; i++)
         {
            if ((hits->in[i].t >= TAG_T_MIN_NS) &&
                (hits->in[i].t <= TAG_T_MAX_NS))
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
            int m = box->fixedChannels->mult++;
            box->fixedChannels->in[m] = channels->in[0];
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

   microCount = 0;
   fixedCount = 0;

   if ((box->microChannels != HDDM_NULL) &&
       (box->microChannels->mult == 0))
   {
      FREE(box->microChannels);
      box->microChannels = HDDM_NULL;
   }
   if ((box->fixedChannels != HDDM_NULL) &&
       (box->fixedChannels->mult == 0))
   {
      FREE(box->fixedChannels);
      box->fixedChannels = HDDM_NULL;
   }
   if (box->microChannels->mult == 0 &&
       box->fixedChannels->mult == 0)
   {
      FREE(box);
      box = HDDM_NULL;
   }
   return box;
}
