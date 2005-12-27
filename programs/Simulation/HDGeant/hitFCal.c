/*
 * hitFCal - registers hits for forward calorimeter
 *
 *	This is a part of the hits package for the
 *	HDGeant simulation program for Hall D.
 *
 *	version 1.0 	-Richard Jones July 16, 2001
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <hddm_s.h>
#include <geant3.h>
#include <bintree.h>

#define TWO_HIT_RESOL   75.
#define MAX_HITS        100
#define THRESH_MEV      30.


binTree_t* forwardEMcalTree = 0;
static int blockCount = 0;
static int showerCount = 0;


/* register hits during tracking (from gustep) */

void hitForwardEMcal (float xin[4], float xout[4],
                      float pin[5], float pout[5], float dEsum,
                      int track, int stack)
{
   float x[3], t;
   float xfcal[3];
   float zeroHat[] = {0,0,0};

   if (dEsum == 0) return;              /* only seen if it deposits energy */

   x[0] = (xin[0] + xout[0])/2;
   x[1] = (xin[1] + xout[1])/2;
   x[2] = (xin[2] + xout[2])/2;
   t    = (xin[3] + xout[3])/2 * 1e9;
   transformCoord(zeroHat,"local",xfcal,"FCAL");

   /* post the hit to the hits tree, mark block as hit */
   {
      int nhit;
      s_FcalHits_t* hits;
      int row = getrow_();
      int column = getcolumn_();
      int mark = (row<<16) + column;
      void** twig = getTwig(&forwardEMcalTree, mark);
      if (*twig == 0)
      {
         s_ForwardEMcal_t* cal = *twig = make_s_ForwardEMcal();
         s_FcalBlocks_t* blocks = make_s_FcalBlocks(1);
         blocks->mult = 1;
         blocks->in[0].row = row;
         blocks->in[0].column = column;
         blocks->in[0].fcalHits = hits = make_s_FcalHits(MAX_HITS);
         cal->fcalBlocks = blocks;
         blockCount++;
      }
      else
      {
         s_ForwardEMcal_t* cal = *twig;
         hits = cal->fcalBlocks->in[0].fcalHits;
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
         hits->in[nhit].t =
                       (hits->in[nhit].t * hits->in[nhit].E + t*dEsum)
                     / (hits->in[nhit].E += dEsum);
      }
      else if (nhit < MAX_HITS)         /* create new hit */
      {
         hits->in[nhit].t = t;
         hits->in[nhit].E = dEsum;
         hits->mult++;
      }
      else
      {
         fprintf(stderr,"HDGeant error in hitforwardEMcal: ");
         fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
         exit(2);
      }
   }

   /* post the hit to the truth tree, once per primary track */
   {
      s_FcalTruthShowers_t* showers;
      int mark = (1<<30) + track;
      void** twig = getTwig(&forwardEMcalTree, mark);
      if (*twig == 0)
      {
         s_ForwardEMcal_t* cal = *twig = make_s_ForwardEMcal();
         cal->fcalTruthShowers = showers = make_s_FcalTruthShowers(1);
         showers->in[0].primary = (stack == 0);
         showers->in[0].track = track;
         showers->in[0].t = t;
         showers->in[0].x = x[0];
         showers->in[0].y = x[1];
         showers->in[0].z = x[2];
         showers->in[0].E = dEsum;
         showers->mult = 1;
         showerCount++;
      }
      else
      {
         showers = ((s_ForwardEMcal_t*) *twig)->fcalTruthShowers;
         showers->in[0].x = (showers->in[0].x * showers->in[0].E + x[0]*dEsum)
                          / (showers->in[0].E + dEsum);
         showers->in[0].y = (showers->in[0].y * showers->in[0].E + x[1]*dEsum)
                          / (showers->in[0].E + dEsum);
         showers->in[0].z = (showers->in[0].z * showers->in[0].E + x[2]*dEsum)
                          / (showers->in[0].E + dEsum);
         showers->in[0].t = (showers->in[0].t * showers->in[0].E + t*dEsum)
                          / (showers->in[0].E += dEsum);
      }
   }
}

/* entry point from fortran */

void hitforwardemcal_(float* xin, float* xout,
                      float* pin, float* pout, float* dEsum,
                      int* track, int* stack)
{
   hitForwardEMcal(xin,xout,pin,pout,*dEsum,*track,*stack);
}


/* pick and package the hits for shipping */

s_ForwardEMcal_t* pickForwardEMcal ()
{
   s_ForwardEMcal_t* box;
   s_ForwardEMcal_t* item;

   if ((blockCount == 0) && (showerCount == 0))
   {
      return HDDM_NULL;
   }

   box = make_s_ForwardEMcal();
   box->fcalBlocks = make_s_FcalBlocks(blockCount);
   box->fcalTruthShowers = make_s_FcalTruthShowers(showerCount);
   while (item = (s_ForwardEMcal_t*) pickTwig(&forwardEMcalTree))
   {
      s_FcalBlocks_t* blocks = item->fcalBlocks;
      s_FcalTruthShowers_t* showers = item->fcalTruthShowers;
      if (blocks != HDDM_NULL)
      {
         s_FcalHits_t* hits = blocks->in[0].fcalHits;
	 int m = box->fcalBlocks->mult;
         int mok = 0;

         /* compress out the hits below threshold */
         if (hits != HDDM_NULL)
         {
            int i;
            for (i=0; i < hits->mult; i++)
            {
               if (hits->in[i].E >= THRESH_MEV/1e3)
               {
                  if (mok < i)
                  {
                     hits->in[mok] = hits->in[i];
                  }
                  ++mok;
               }
            }
            hits->mult = mok;
            if (mok == 0)
            {
               blocks->in[0].fcalHits = HDDM_NULL;
               FREE(hits);
            }
         }

         if (mok)
         {
            box->fcalBlocks->in[m] = blocks->in[0];
            box->fcalBlocks->mult++;
         }
         FREE(blocks);
      }
      else if (showers != HDDM_NULL)
      {
         int m = box->fcalTruthShowers->mult++;
         box->fcalTruthShowers->in[m] = showers->in[0];
         FREE(showers);
      }
      FREE(item);
   }
   blockCount = showerCount = 0;

   return box;
}
