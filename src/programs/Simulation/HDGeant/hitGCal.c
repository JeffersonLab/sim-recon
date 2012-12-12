/*
 * hitGCal - registers hits for gap calorimeter
 *
 *	This is a part of the hits package for the
 *	HDGeant simulation program for Hall D.
 *
 *	version 1.0 	-Richard Jones July 16, 2001
 *
 * changes: Wed Jun 20 13:19:56 EDT 2007 B. Zihlmann 
 *          add ipart to the function hitGapEMcal
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <hddm_s.h>
#include <geant3.h>
#include <bintree.h>
extern s_HDDM_t* thisInputEvent;

//#define ATTEN_LENGTH	100.
#define ATTEN_LENGTH	1e6
#define C_EFFECTIVE	15.
#define WIDTH_OF_BLOCK  4.
#define LENGTH_OF_BLOCK 45.
#define TWO_HIT_RESOL   75.
#define MAX_HITS        100
#define THRESH_MEV      30.
#define ACTIVE_RADIUS   120.e6
#define CENTRAL_ROW     29
#define CENTRAL_COLUMN  29


binTree_t* gapEMcalTree = 0;
static int cellCount = 0;
static int showerCount = 0;


/* register hits during tracking (from gustep) */

void hitGapEMcal (float xin[4], float xout[4],
                      float pin[5], float pout[5], float dEsum,
                      int track, int stack, int history, int ipart)
{
   float x[3], t;
   float xgcal[3];
   float zeroHat[] = {0,0,0};

   x[0] = (xin[0] + xout[0])/2;
   x[1] = (xin[1] + xout[1])/2;
   x[2] = (xin[2] + xout[2])/2;
   t    = (xin[3] + xout[3])/2 * 1e9;
   transformCoord(zeroHat,"local",xgcal,"gCAL");

   /* post the hit to the truth tree */

   if ((history == 0) && (pin[3] > THRESH_MEV/1e3))
   {
      s_GcalTruthShowers_t* showers;
      float r = sqrt(xin[0]*xin[0]+xin[1]*xin[1]);
      float phi = atan2(xin[1],xin[0]);
      int mark = (1<<30) + showerCount;
      void** twig = getTwig(&gapEMcalTree, mark);
      if (*twig == 0)
      {
         s_GapEMcal_t* cal = *twig = make_s_GapEMcal();
         cal->gcalTruthShowers = showers = make_s_GcalTruthShowers(1);
        int a = thisInputEvent->physicsEvents->in[0].reactions->in[0].vertices->in[0].products->mult;
         showers->in[0].primary = (stack <= a);
         showers->in[0].track = track;
         showers->in[0].z = xin[2];
         showers->in[0].r = r;
         showers->in[0].phi = phi;
         showers->in[0].t = xin[3]*1e9;
         showers->in[0].px = pin[0]*pin[4];
         showers->in[0].py = pin[1]*pin[4];
         showers->in[0].pz = pin[2]*pin[4];
         showers->in[0].E = pin[3];
         showers->in[0].ptype = ipart;
         showers->mult = 1;
         showerCount++;
      }
   }

   /* post the hit to the hits tree, mark block as hit */

   if (dEsum > 0)
   {
      int nhit;
      s_GcalHits_t* hits;
      int module = getmodule_wrapper_();
      float dist = LENGTH_OF_BLOCK-xgcal[2];
      float dEcorr = dEsum * exp(-dist/ATTEN_LENGTH);
      float tcorr = t + dist/C_EFFECTIVE;
      int mark = ((module+1)<<16);
      void** twig = getTwig(&gapEMcalTree, mark);
      if (*twig == 0)
      {
         s_GapEMcal_t* cal = *twig = make_s_GapEMcal();
         s_GcalCells_t* cells = make_s_GcalCells(1);
         cells->mult = 1;
         cells->in[0].module = module;
         cells->in[0].gcalHits = hits = make_s_GcalHits(MAX_HITS);
         cal->gcalCells = cells;
         cellCount++;
      }
      else
      {
         s_GapEMcal_t* cal = *twig;
         hits = cal->gcalCells->in[0].gcalHits;
      }

      for (nhit = 0; nhit < hits->mult; nhit++)
      {
         if (fabs(hits->in[nhit].t - tcorr) < TWO_HIT_RESOL)
         {
            break;
         }
      }
      if (nhit < hits->mult)		/* merge with former hit */
      {
         hits->in[nhit].t =
                       (hits->in[nhit].t * hits->in[nhit].E + tcorr*dEcorr)
                     / (hits->in[nhit].E += dEcorr);
      }
      else if (nhit < MAX_HITS)         /* create new hit */
      {
         hits->in[nhit].t = tcorr;
         hits->in[nhit].E = dEcorr;
         hits->mult++;
      }
      else
      {
         fprintf(stderr,"HDGeant error in hitgapEMcal: ");
         fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
         exit(2);
      }
   }
}

/* entry point from fortran */

void hitgapemcal_(float* xin, float* xout,
                      float* pin, float* pout, float* dEsum,
                      int* track, int* stack, int* history, int* ipart)
{
   hitGapEMcal(xin,xout,pin,pout,*dEsum,*track,*stack,*history, *ipart);
}


/* pick and package the hits for shipping */

s_GapEMcal_t* pickGapEMcal ()
{
   s_GapEMcal_t* box;
   s_GapEMcal_t* item;

#if TESTING_CAL_CONTAINMENT
   double Etotal = 0;
#endif
   if ((cellCount == 0) && (showerCount == 0))
   {
      return HDDM_NULL;
   }

   box = make_s_GapEMcal();
   box->gcalCells = make_s_GcalCells(cellCount);
   box->gcalTruthShowers = make_s_GcalTruthShowers(showerCount);
   while (item = (s_GapEMcal_t*) pickTwig(&gapEMcalTree))
   {
      s_GcalCells_t* cells = item->gcalCells;
      int cell;
      s_GcalTruthShowers_t* showers = item->gcalTruthShowers;
      int shower;
      for (cell=0; cell < cells->mult; ++cell)
      {
	 int m = box->gcalCells->mult;
         int mok = 0;

         s_GcalHits_t* hits = cells->in[cell].gcalHits;

         /* compress out the hits below threshold */
         int i,iok;
         for (iok=i=0; i < hits->mult; i++)
         {
            if (hits->in[i].E >= THRESH_MEV/1e3)
            {
#if TESTING_CAL_CONTAINMENT
  Etotal += hits->in[i].E;
#endif
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
            box->gcalCells->in[m] = cells->in[cell];
            box->gcalCells->mult++;
         }
         else if (hits != HDDM_NULL)
         {
            FREE(hits);
         }
         if (hits != HDDM_NULL)
         {
            FREE(hits);
         }
      }

      for (shower=0; shower < showers->mult; ++shower)
      {
         int m = box->gcalTruthShowers->mult++;
         box->gcalTruthShowers->in[m] = showers->in[shower];
      }
      if (cells != HDDM_NULL)
      {
         FREE(cells);
      }
      if (showers != HDDM_NULL)
      {
         FREE(showers);
      }
      FREE(item);
   }

   cellCount = showerCount = 0;

   if ((box->gcalCells != HDDM_NULL) &&
       (box->gcalCells->mult == 0))
   {
      FREE(box->gcalCells);
      box->gcalCells = HDDM_NULL;
   }
   if ((box->gcalTruthShowers != HDDM_NULL) &&
       (box->gcalTruthShowers->mult == 0))
   {
      FREE(box->gcalTruthShowers);
      box->gcalTruthShowers = HDDM_NULL;
   }
   if ((box->gcalCells->mult == 0) &&
       (box->gcalTruthShowers->mult == 0))
   {
      FREE(box);
      box = HDDM_NULL;
   }
#if TESTING_CAL_CONTAINMENT
  printf("GCal energy sum: %f\n",Etotal);
#endif
   return box;
}
