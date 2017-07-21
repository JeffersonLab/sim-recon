/*
 * hitTPOL - registers hits for triplet polarimeter silicon detector
 *
 *        This is a part of the hits package for the
 *        HDGeant simulation program for Hall D.
 *
 *        version 1.0         -Richard Jones, Jan 14, 2017
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

static float TWO_HIT_RESOL   = 1000.;
static float THRESH_MEV      = 0.010;

// Comment by RTJ:
// When I introduced the convenience constant MAX_HITS,
// I never intended it to be a tunable simulation parameter.
// Do not use it as such.  Do NOT MODIFY it, or the way
// it functions in the algorithm.  If you want to truncate
// the hit list, do it in mcsmear.
#define MAX_HITS 100
 
binTree_t* tpolTree = 0;
static int sectorCount = 0;
static int pointCount = 0;
static int initialized = 0;


/* register hits during tracking (from gustep) */

void hitTPOL(float xin[4], float xout[4],float pin[5], float pout[5], float dEsum,
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

   int itrack = (stack == 0)? gidGetId(track) : -1;

   if (history == 0)
   {
      int mark = (1<<30) + pointCount;
      void** twig = getTwig(&tpolTree, mark);
      if (*twig == 0)
      {
         s_TripletPolarimeter_t* tpol = *twig = make_s_TripletPolarimeter();
         s_TpolTruthPoints_t* points = make_s_TpolTruthPoints(1);
         tpol->tpolTruthPoints = points;
         int a = thisInputEvent->physicsEvents->in[0].reactions->in[0].vertices->in[0].products->mult;
         points->in[0].primary = (track <= a && stack == 0);
         points->in[0].track = track;
         points->in[0].t = t;
         points->in[0].r = sqrt(x[0]*x[0] + x[1]*x[1]);
         points->in[0].phi = atan2(x[1], x[0]);
         points->in[0].px = pin[0]*pin[4];
         points->in[0].py = pin[1]*pin[4];
         points->in[0].pz = pin[2]*pin[4];
         points->in[0].E = pin[3];
         points->in[0].dEdx = dEdx;
         points->in[0].ptype = ipart;
         points->in[0].trackID = make_s_TrackID();
         points->in[0].trackID->itrack = itrack;
         points->mult = 1;
         pointCount++;
      }
   }

   /* post the hit to the hits tree, mark module as hit */
   if (dEsum > 0)
   {
      int nhit;
      s_TpolTruthHits_t* hits;
      int ringno = 0; //getring_wrapper_();
      int sectno = getsector_wrapper_();
      int mark = sectno;
      void** twig = getTwig(&tpolTree, mark);
      if (*twig == 0)
      {
         s_TripletPolarimeter_t* tpol = *twig = make_s_TripletPolarimeter();
         s_TpolSectors_t* sectors = make_s_TpolSectors(1);
         sectors->mult = 1;
         sectors->in[0].ring = ringno;
         sectors->in[0].sector = sectno;
         sectors->in[0].tpolTruthHits = hits = make_s_TpolTruthHits(MAX_HITS);
         tpol->tpolSectors = sectors;
         sectorCount++;
      }
      else
      {
         s_TripletPolarimeter_t* tpol = *twig;
         hits = tpol->tpolSectors->in[0].tpolTruthHits;
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
            hits->in[nhit].itrack = itrack;
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
         hits->in[nhit].itrack = itrack;
         hits->mult++;
      }
      else
      {
         fprintf(stderr,"HDGeant error in hitTPOL: ");
         fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
         exit(2);
      }
   }
}

/* entry point from fortran */

void hittpol_(float* xin, float* xout,
                   float* pin, float* pout, float* dEsum,
                   int* track, int* stack, int* history, int* ipart)
{
   hitTPOL(xin,xout,pin,pout,*dEsum,*track,*stack,*history,*ipart);
}


/* pick and package the hits for shipping */

s_TripletPolarimeter_t* pickTpol ()
{
   s_TripletPolarimeter_t* box;
   s_TripletPolarimeter_t* item;

   if ((sectorCount == 0) && (pointCount == 0))
   {
      return HDDM_NULL;
   }

   box = make_s_TripletPolarimeter();
   box->tpolSectors = make_s_TpolSectors(sectorCount);
   box->tpolTruthPoints = make_s_TpolTruthPoints(pointCount);
   while ((item = (s_TripletPolarimeter_t*) pickTwig(&tpolTree)))
   {
      s_TpolSectors_t* sectors = item->tpolSectors;
      int sector;
      s_TpolTruthPoints_t* points = item->tpolTruthPoints;
      int point;

      for (sector=0; sector < sectors->mult; ++sector)
      {
         int m = box->tpolSectors->mult;

         s_TpolTruthHits_t* hits = sectors->in[sector].tpolTruthHits;

         /* compress out the hits below threshold */
         int i,iok;
         for (iok=i=0; i < hits->mult; i++)
         {
            if (hits->in[i].dE > THRESH_MEV/1e3)
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
            box->tpolSectors->in[m] = sectors->in[sector];
            box->tpolSectors->mult++;
         }
         else if (hits != HDDM_NULL)
         {
            FREE(hits);
         }
      }
      if (sectors != HDDM_NULL)
      {
         FREE(sectors);
      }

      for (point=0; point < points->mult; ++point)
      {
         int track = points->in[point].track;
         double t = points->in[point].t;
         int m = box->tpolTruthPoints->mult;
         if (points->in[point].trackID->itrack < 0 ||
            (m > 0 &&  box->tpolTruthPoints->in[m-1].track == track &&
             fabs(box->tpolTruthPoints->in[m-1].t - t) < 0.5))
         {
            FREE(points->in[point].trackID);
            continue;
         }
         box->tpolTruthPoints->in[m] = item->tpolTruthPoints->in[point];
         box->tpolTruthPoints->mult++;
      }
      if (points != HDDM_NULL)
      {
         FREE(points);
      }
      FREE(item);
   }

   sectorCount = pointCount = 0;

   if ((box->tpolSectors != HDDM_NULL) &&
       (box->tpolSectors->mult == 0))
   {
      FREE(box->tpolSectors);
      box->tpolSectors = HDDM_NULL;
   }
   if ((box->tpolTruthPoints != HDDM_NULL) &&
       (box->tpolTruthPoints->mult == 0))
   {
      FREE(box->tpolTruthPoints);
      box->tpolTruthPoints = HDDM_NULL;
   }
   if ((box->tpolSectors->mult == 0) &&
       (box->tpolTruthPoints->mult == 0))
   {
      FREE(box);
      box = HDDM_NULL;
   }
   return box;
}
