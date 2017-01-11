/*
 * hitPSC - registers hits for Pair Spectrometer Coarse paddles
 *
 *        This is a part of the hits package for the
 *        HDGeant simulation program for Hall D.
 *
 *        version 1.0         -Simon Taylor, Oct 16, 2014
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <HDDM/hddm_s.h>
#include <geant3.h>
#include <bintree.h>
#include <gid_map.h>
#include "calibDB.h"
extern s_HDDM_t* thisInputEvent;

//static float ATTEN_LENGTH    = 150.;
//static float C_EFFECTIVE     = 15.;
static float TWO_HIT_RESOL   = 25.;
static float THRESH_MEV      = 0.010;

// the coarse PS has two arms (north/south) of 8 modules each
#define NUM_MODULES_PER_ARM 8 

// Comment by RTJ:
// When I introduced the convenience constant MAX_HITS,
// I never intended it to be a tunable simulation parameter.
// Do not use it as such.  Do NOT MODIFY it, or the way
// it functions in the algorithm.  If you want to truncate
// the hit list, do it in mcsmear.
#define MAX_HITS 100
 
binTree_t* pscTree = 0;
static int paddleCount = 0;
static int pointCount = 0;
static int initialized = 0;


/* register hits during tracking (from gustep) */

void hitPSC(float xin[4], float xout[4],float pin[5], float pout[5], float dEsum,
           int track, int stack, int history, int ipart)
{
   float x[3], t;
   float dx[3], dr;
   float dEdx;
   float xlocal[3];

   if (!initialized) {
     //  Get calibration constants ...  
     
     initialized = 1;
   }

   x[0] = (xin[0] + xout[0])/2;
   x[1] = (xin[1] + xout[1])/2;
   x[2] = (xin[2] + xout[2])/2;
   t    = (xin[3] + xout[3])/2 * 1e9;
   transformCoord(x,"global",xlocal,"local");
   dx[0] = xin[0] - xout[0];
   dx[1] = xin[1] - xout[1];
   dx[2] = xin[2] - xout[2];
   dr = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
   if (dr > 1e-3)
   {
      dEdx = dEsum/dr;
   }
   else
   {
      dEdx = 0;
   }

   if (history == 0)
   {
      int mark = (1<<30) + pointCount;
      void** twig = getTwig(&pscTree, mark);
      if (*twig == 0)
      {
         s_PairSpectrometerCoarse_t* psc = *twig = make_s_PairSpectrometerCoarse();
         s_PscTruthPoints_t* points = make_s_PscTruthPoints(1);
         psc->pscTruthPoints = points;
         int a = thisInputEvent->physicsEvents->in[0].reactions->in[0].vertices->in[0].products->mult;
         points->in[0].primary = (track <= a && stack == 0);
         points->in[0].track = track;
         points->in[0].t = t;
         points->in[0].z = x[2];
         points->in[0].x = x[0];
         points->in[0].y = x[1];
         points->in[0].px = pin[0]*pin[4];
         points->in[0].py = pin[1]*pin[4];
         points->in[0].pz = pin[2]*pin[4];
         points->in[0].E = pin[3];
         points->in[0].dEdx = dEdx;
         points->in[0].ptype = ipart;
         points->in[0].arm = getmodule_wrapper_() / NUM_MODULES_PER_ARM;
         points->in[0].module = getmodule_wrapper_() % NUM_MODULES_PER_ARM;
         points->in[0].trackID = make_s_TrackID();
         points->in[0].trackID->itrack = gidGetId(track);
         points->mult = 1;
         pointCount++;
      }
   }

   /* post the hit to the hits tree, mark module as hit */
   if (dEsum > 0)
   {
      int nhit;
      s_PscTruthHits_t* hits;
      int module = getmodule_wrapper_();
      int mark = module;
      void** twig = getTwig(&pscTree, mark);
      if (*twig == 0)
      {
         s_PairSpectrometerCoarse_t* psc = *twig = make_s_PairSpectrometerCoarse();
         s_PscPaddles_t* paddles = make_s_PscPaddles(1);
         paddles->mult = 1;
         paddles->in[0].arm = module / NUM_MODULES_PER_ARM;
         paddles->in[0].module = module % NUM_MODULES_PER_ARM;
         paddles->in[0].pscTruthHits = hits = make_s_PscTruthHits(MAX_HITS);
         psc->pscPaddles = paddles;
         paddleCount++;
      }
      else
      {
         s_PairSpectrometerCoarse_t* psc = *twig;
         hits = psc->pscPaddles->in[0].pscTruthHits;
      }

      for (nhit = 0; nhit < hits->mult; nhit++)
      {
         if (fabs(hits->in[nhit].t - t) < TWO_HIT_RESOL)
         {
            break;
         }
      }
      if (nhit < hits->mult)                /* merge with former hit */
      {
         if (t < hits->in[nhit].t)
         {
            hits->in[nhit].ptype = ipart;
            hits->in[nhit].itrack = gidGetId(track);
         }
         hits->in[nhit].t = 
                 (hits->in[nhit].t * hits->in[nhit].dE + t * dEsum) /
                 (hits->in[nhit].dE + dEsum);
                        hits->in[nhit].dE += dEsum;
      }
      else if (nhit < MAX_HITS)                /* create new hit */
      {
         hits->in[nhit].t = t;
         hits->in[nhit].dE = dEsum;
         hits->in[nhit].ptype = ipart;
         hits->in[nhit].itrack = gidGetId(track);
         hits->mult++;
      }
      else
      {
         fprintf(stderr,"HDGeant error in hitPSC: ");
         fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
         exit(2);
      }
   }
}

/* entry point from fortran */

void hitpsc_(float* xin, float* xout,
                   float* pin, float* pout, float* dEsum,
                   int* track, int* stack, int* history, int* ipart)
{
   hitPSC(xin,xout,pin,pout,*dEsum,*track,*stack,*history,*ipart);
}


/* pick and package the hits for shipping */

s_PairSpectrometerCoarse_t* pickPsc ()
{
   s_PairSpectrometerCoarse_t* box;
   s_PairSpectrometerCoarse_t* item;

   if ((paddleCount == 0) && (pointCount == 0))
   {
      return HDDM_NULL;
   }

   box = make_s_PairSpectrometerCoarse();
   box->pscPaddles = make_s_PscPaddles(paddleCount);
   box->pscTruthPoints = make_s_PscTruthPoints(pointCount);
   while ((item = (s_PairSpectrometerCoarse_t*) pickTwig(&pscTree)))
   {
      s_PscPaddles_t* paddles = item->pscPaddles;
      int paddle;
      s_PscTruthPoints_t* points = item->pscTruthPoints;
      int point;

      for (paddle=0; paddle < paddles->mult; ++paddle)
      {
         int m = box->pscPaddles->mult;

         s_PscTruthHits_t* hits = paddles->in[paddle].pscTruthHits;

         /* compress out the hits below threshold */
         int i,iok;
         for (iok=i=0; i < hits->mult; i++)
         {
            if (hits->in[i].dE >= THRESH_MEV/1e3)
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
            box->pscPaddles->in[m] = paddles->in[paddle];
            box->pscPaddles->mult++;
         }
         else if (hits != HDDM_NULL)
         {
            FREE(hits);
         }
      }
      if (paddles != HDDM_NULL)
      {
         FREE(paddles);
      }

      int last_track = -1;
      double last_t = 1e9;
      for (point=0; point < points->mult; ++point)
      {
         if (points->in[point].trackID->itrack > 0 &&
            (points->in[point].track != last_track ||
             fabs(points->in[point].t - last_t) > 0.1))
         {
            int m = box->pscTruthPoints->mult++;
            box->pscTruthPoints->in[m] = item->pscTruthPoints->in[point];
            last_track = points->in[point].track;
            last_t = points->in[point].t;
         }
      }
      if (points != HDDM_NULL)
      {
         FREE(points);
      }
      FREE(item);
   }

   paddleCount = pointCount = 0;

   if ((box->pscPaddles != HDDM_NULL) &&
       (box->pscPaddles->mult == 0))
   {
      FREE(box->pscPaddles);
      box->pscPaddles = HDDM_NULL;
   }
   if ((box->pscTruthPoints != HDDM_NULL) &&
       (box->pscTruthPoints->mult == 0))
   {
      FREE(box->pscTruthPoints);
      box->pscTruthPoints = HDDM_NULL;
   }
   if ((box->pscPaddles->mult == 0) &&
       (box->pscTruthPoints->mult == 0))
   {
      FREE(box);
      box = HDDM_NULL;
   }
   return box;
}
