/*
 * hitStart - registers hits for Start counter
 *
 *	This is a part of the hits package for the
 *	HDGeant simulation program for Hall D.
 *
 *	version 1.0 	-Richard Jones July 16, 2001
 */

#include <stdlib.h>
#include <stdio.h>

#include <hddm_s.h>
#include <geant3.h>
#include <bintree.h>

#define Z0		-45.
#define ATTEN_LENGTH	150.
#define C_EFFECTIVE	15.
#define TWO_HIT_RESOL	25.
#define MAX_HITS 	10

binTree_t* startCntrTree = 0;
static int paddleCount = 0;
static int pointCount = 0;


/* register hits during tracking (from gustep) */

void hitStartCntr (float xin[4], float xout[4],
                   float pin[5], float pout[5], float dEsum, int track)
{
   float x[3], t;
   float dx[3], dr;
   float dEdx;
   float xlocal[3];

   if (dEsum == 0) return;		/* only seen if it deposits energy */

   x[0] = (xin[0] + xout[0])/2;
   x[1] = (xin[1] + xout[1])/2;
   x[2] = (xin[2] + xout[2])/2;
   t    = (xin[3] + xout[3])/2 * 1e9;
   getlocalcoord_(x,xlocal,"VRTX",4);
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

   /* post the hit to the hits tree, mark sector as hit */
   {
      int nhit;
      s_Hits_t* hits;
      int count = getcount_();
      int sector = getsector_();
      float phim = (sector - 0.5) * 2*M_PI/count;
      float dzup = xlocal[2] - Z0;
      float tcorr = t + dzup/C_EFFECTIVE;
      float dEcorr = dEsum * exp(-dzup/ATTEN_LENGTH);
      void** twig = getTwig(&startCntrTree, sector);
      if (*twig == 0)
      {
         s_StartCntr_t* vtx = *twig = make_s_StartCntr();
         vtx->paddles = make_s_Paddles(1);
         vtx->paddles->mult = 1;
         vtx->paddles->in[0].phim = phim;
         vtx->paddles->in[0].hits = hits = make_s_Hits(MAX_HITS);
         paddleCount++;
      }
      else
      {
         s_StartCntr_t* vtx = *twig;
         hits = vtx->paddles->in[0].hits;
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
                      (hits->in[nhit].t * hits->in[nhit].dE + tcorr * dEcorr)
                    / (hits->in[nhit].dE += dEcorr);
      }
      else if (nhit < MAX_HITS)		/* create new hit */
      {
         hits->in[nhit].t = tcorr ;
         hits->in[nhit].dE = dEcorr;
         hits->mult++;
      }
      else
      {
         fprintf(stderr,"HDGeant error in hitStart: ");
         fprintf(stderr,"maximum hit count %d exceeded, quitting!\n",MAX_HITS);
         exit(2);
      }
   }

   /* post the hit to the truth tree, once per primary track */
   {
      int mark = (track << 16);
      void** twig = getTwig(&startCntrTree, mark);
      if (*twig == 0)
      {
         s_StartCntr_t* vtx = *twig = make_s_StartCntr();
         s_StartPoints_t* points = make_s_StartPoints(1);
         vtx->startPoints = points;
         points->in[0].t = t;
         points->in[0].z = x[2];
         points->in[0].r = sqrt(x[0]*x[0] + x[1]*x[1]);
         points->in[0].phi = atan2(x[1],x[0]);
         points->in[0].dEdx = dEdx;
         points->mult = 1;
         pointCount++;
      }
   }
}

/* entry point from fortran */

void hitstartcntr_(float* xin, float* xout,
                   float* pin, float* pout, float* dEsum, int* track)
{
   hitStartCntr(xin,xout,pin,pout,*dEsum,*track);
}


/* pick and package the hits for shipping */

s_StartCntr_t* pickStartCntr ()
{
   s_StartCntr_t* box;
   s_StartCntr_t* item;

   if ((paddleCount == 0) && (pointCount == 0))
   {
      return 0;
   }

   box = make_s_StartCntr();
   box->paddles = make_s_Paddles(paddleCount);
   box->startPoints = make_s_StartPoints(pointCount);
   while (item = (s_StartCntr_t*) pickTwig(&startCntrTree))
   {
      if (item->paddles)
      {
         int m = box->paddles->mult++;
         box->paddles->in[m] = item->paddles->in[0];
         FREE(item->paddles);
      }
      else if (item->startPoints)
      {
         int m = box->startPoints->mult++;
         box->startPoints->in[m] = item->startPoints->in[0];
         FREE(item->startPoints);
      }
      FREE(item);
   }
   paddleCount = pointCount = 0;
   return box;
}
