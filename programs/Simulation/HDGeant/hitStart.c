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
#include <hittree.h>

#define MAX_HITS 	10
#define TWO_HIT_RESOL	25.

hitTree_t* startCntrTree = 0;

void hitStartCntr (float xin[4], float xout[4],
                   float pin[5], float pout[5], float dEsum, int track)
{
   float x[3], t;
   float dx[3], dr;
   float dEdx;
   float xlocal[3];
   const int one = 1;

   if (dEsum == 0) return;		/* only seen if it deposits energy */

   gmtod_(x,xlocal,&one);
   x[0] = (xin[0] + xout[0])/2;
   x[1] = (xin[1] + xout[1])/2;
   x[2] = (xin[2] + xout[2])/2;
   t    = (xin[3] + xout[3])/2 * 1e9;
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
      int sector = getsector_();
      void** twig = getTwig(&startCntrTree, sector);
      if (twig == 0)
      {
         s_StartCntr_t* vtx = *twig = make_s_StartCntr();
         vtx->paddles = make_s_Paddles(1);
         vtx->paddles->mult = 1;
         vtx->paddles->in[0].hits = hits = make_s_Hits(MAX_HITS);
      }
      else
      {
         s_StartCntr_t* vtx = *twig;
         hits = vtx->paddles->in[0].hits;
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
         hits->in[nhit].t = (hits->in[nhit].t * hits->in[nhit].dE + t*dEsum)
                          / (hits->in[nhit].dE += dEsum);
      }
      else if (nhit < MAX_HITS)		/* create new hit */
      {
         hits->in[nhit].t = t;
         hits->in[nhit].dE = dEsum;
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
      if (twig == 0)
      {
         s_StartCntr_t* vtx = *twig = make_s_StartCntr();
         s_StartPoints_t* points = make_s_StartPoints(1);
         vtx->startPoints = points;
         points->in[0].z = x[2];
         points->in[0].r = sqrt(x[0]*x[0] + x[1]*x[1]);
         points->in[0].phi = atan2(x[1],x[0]);
         points->in[0].dEdx = dEdx;
         points->mult = 1;
      }
   }
}

/* entry point from fortran */

void hitstartcntr_(float* xin, float* xout,
                   float* pin, float* pout, float* dEsum, int* track)
{
   hitStartCntr(xin,xout,pin,pout,*dEsum,*track);
}
