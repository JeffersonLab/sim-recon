/*
 * hitPS - registers hits for Pair Spectrometer
 *
 *	This is a part of the hits package for the
 *	HDGeant simulation program for Hall D.
 *
 *	version 1.0 	-Simon Taylor, Oct 16, 2014
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

// the fine PS has two arms (north/south) of 145 columns each
#define NUM_COLUMN_PER_ARM 145 

// Comment by RTJ:
// When I introduced the convenience constant MAX_HITS,
// I never intended it to be a tunable simulation parameter.
// Do not use it as such.  Do NOT MODIFY it, or the way
// it functions in the algorithm.  If you want to truncate
// the hit list, do it in mcsmear.
#define MAX_HITS 100
 
binTree_t* psTree = 0;
static int tileCount = 0;
static int pointCount = 0;
static int initialized = 0;


/* register hits during tracking (from gustep) */

void hitPS(float xin[4], float xout[4],float pin[5], float pout[5], float dEsum,
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
      void** twig = getTwig(&psTree, mark);
      if (*twig == 0)
      {
         s_PairSpectrometerFine_t* ps = *twig = make_s_PairSpectrometerFine();
         s_PsTruthPoints_t* points = make_s_PsTruthPoints(1);
         ps->psTruthPoints = points;
         int a = thisInputEvent->physicsEvents->in[0].reactions->in[0].vertices->in[0].products->mult;
         points->in[0].primary = (stack <= a);
         points->in[0].track = track;
         points->in[0].t = t;
         points->in[0].z = x[2];
         points->in[0].x = x[0];
	 points->in[0].y = x[1];
         points->in[0].phi = atan2(x[1],x[0]);
         points->in[0].px = pin[0]*pin[4];
         points->in[0].py = pin[1]*pin[4];
         points->in[0].pz = pin[2]*pin[4];
         points->in[0].E = pin[3];
         points->in[0].dEdx = dEdx;
         points->in[0].ptype = ipart;
	 // the fine PS has two arms: North/South (0/1)
         points->in[0].arm = getcolumn_wrapper_() / NUM_COLUMN_PER_ARM;
         points->in[0].column = getcolumn_wrapper_() % NUM_COLUMN_PER_ARM;
	 points->in[0].trackID = make_s_TrackID();
	 points->in[0].trackID->itrack = gidGetId(track);
	 points->mult = 1;
         pointCount++;


      }
   }

   /* post the hit to the hits tree, mark column as hit */
   if (dEsum > 0)
   {
      int nhit;
      s_PsTruthHits_t* hits;
      int column = getcolumn_wrapper_();
      int mark = column;
      void** twig = getTwig(&psTree, mark);
      if (*twig == 0)
      {
         s_PairSpectrometerFine_t* ps = *twig = make_s_PairSpectrometerFine();
         s_PsTiles_t* tiles = make_s_PsTiles(1);
         tiles->mult = 1;
	 // the fine PS has two arms: North/South (0/1)
         tiles->in[0].arm = column / NUM_COLUMN_PER_ARM;
         tiles->in[0].column = column % NUM_COLUMN_PER_ARM;
         tiles->in[0].psTruthHits = hits = make_s_PsTruthHits(MAX_HITS);
         ps->psTiles = tiles;
         tileCount++;
      }
      else
      {
         s_PairSpectrometerFine_t* ps = *twig;
         hits = ps->psTiles->in[0].psTruthHits;
      }

      for (nhit = 0; nhit < hits->mult; nhit++)
      {
         if (fabs(hits->in[nhit].t - t) < TWO_HIT_RESOL)
         {
            break;
         }
      }
      if (nhit < hits->mult)		/* merge with former hit */
      {
         if (t < hits->in[nhit].t)
         {
            hits->in[nhit].ptype = ipart;
            hits->in[nhit].itrack = gidGetId(track);
         }
         hits->in[nhit].t = 
                 (hits->in[nhit].t * hits->in[nhit].dE + t * dEsum) /
                 (hits->in[nhit].dE += dEsum);
      }
      else if (nhit < MAX_HITS)		/* create new hit */
      {
         hits->in[nhit].t = t;
         hits->in[nhit].dE = dEsum;
         hits->in[nhit].ptype = ipart;
         hits->in[nhit].itrack = gidGetId(track);
         hits->mult++;
      }
      else
      {
         fprintf(stderr,"HDGeant error in hitPS: ");
         fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
         exit(2);
      }
   }
}

/* entry point from fortran */

void hitps_(float* xin, float* xout,
                   float* pin, float* pout, float* dEsum,
                   int* track, int* stack, int* history, int* ipart)
{
   hitPS(xin,xout,pin,pout,*dEsum,*track,*stack,*history,*ipart);
}


/* pick and package the hits for shipping */

s_PairSpectrometerFine_t* pickPs ()
{
   s_PairSpectrometerFine_t* box;
   s_PairSpectrometerFine_t* item;

   if ((tileCount == 0) && (pointCount == 0))
   {
      return HDDM_NULL;
   }

   box = make_s_PairSpectrometerFine();
   box->psTiles = make_s_PsTiles(tileCount);
   box->psTruthPoints = make_s_PsTruthPoints(pointCount);
   while ((item = (s_PairSpectrometerFine_t*) pickTwig(&psTree)))
   {
      s_PsTiles_t* tiles = item->psTiles;
      int tile;
      s_PsTruthPoints_t* points = item->psTruthPoints;
      int point;

      for (tile=0; tile < tiles->mult; ++tile)
      {
         int m = box->psTiles->mult;

         s_PsTruthHits_t* hits = tiles->in[tile].psTruthHits;

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
            box->psTiles->in[m] = tiles->in[tile];
            box->psTiles->mult++;
         }
         else if (hits != HDDM_NULL)
         {
            FREE(hits);
         }
      }
      if (tiles != HDDM_NULL)
      {
         FREE(tiles);
      }

      for (point=0; point < points->mult; ++point)
      {
         int m = box->psTruthPoints->mult++;
         box->psTruthPoints->in[m] = item->psTruthPoints->in[point];
      }
      if (points != HDDM_NULL)
      {
         FREE(points);
      }
      FREE(item);
   }

   tileCount = pointCount = 0;

   if ((box->psTiles != HDDM_NULL) &&
       (box->psTiles->mult == 0))
   {
      FREE(box->psTiles);
      box->psTiles = HDDM_NULL;
   }
   if ((box->psTruthPoints != HDDM_NULL) &&
       (box->psTruthPoints->mult == 0))
   {
      FREE(box->psTruthPoints);
      box->psTruthPoints = HDDM_NULL;
   }
   if ((box->psTiles->mult == 0) &&
       (box->psTruthPoints->mult == 0))
   {
      FREE(box);
      box = HDDM_NULL;
   }
   return box;
}
