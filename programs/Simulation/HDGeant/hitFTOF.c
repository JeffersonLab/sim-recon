/*
 * hitFTOF - registers hits for forward Time-Of-Flight
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

#define MAX_HITS        10
#define Y0		-125.
#define TWO_HIT_RESOL   25.
#define SLAB_WIDTH	5.0
#define ATTEN_LENGTH	150
#define C_EFFECTIVE	15

hitTree_t* forwardTOFTree = 0;

void hitForwardTOF (float xin[4], float xout[4],
                    float pin[5], float pout[5], float dEsum, int track)
{
   float x[3], t;
   float dx[3], dr;
   float dEdx;
   float xlocal[3];
   const int one = 1;

   if (dEsum == 0) return;              /* only seen if it deposits energy */

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

   /* post the hit to the hits tree, mark slab as hit */
   {
      int nhit;
      s_Hits_t* leftHits;
      s_Hits_t* rightHits;
      float tleft  = t - (xlocal[0]/C_EFFECTIVE);
      float tright = t + (xlocal[0]/C_EFFECTIVE);
      float dEleft  = dEsum * exp(+xlocal[0]/ATTEN_LENGTH);
      float dEright = dEsum * exp(-xlocal[0]/ATTEN_LENGTH);
      int slab = (xlocal[1] - Y0)/SLAB_WIDTH +1;
      void** twig = getTwig(&forwardTOFTree, slab);
      if (twig == 0)
      {
         s_ForwardTOF_t* tof = *twig = make_s_ForwardTOF();
         tof->slabs = make_s_Slabs(1);
         tof->slabs->mult = 1;
         tof->slabs->in[0].y = (slab -1)*SLAB_WIDTH + Y0;
         tof->slabs->in[0].left->hits = leftHits = make_s_Hits(MAX_HITS);
         tof->slabs->in[0].right->hits = rightHits = make_s_Hits(MAX_HITS);
      }
      else
      {
         s_ForwardTOF_t* tof = *twig;
         leftHits = tof->slabs->in[0].left->hits;
         rightHits = tof->slabs->in[0].right->hits;
      }

      for (nhit = 0; nhit < leftHits->mult; nhit++)
      {
         if (fabs(leftHits->in[nhit].t - t) < TWO_HIT_RESOL)
         {
            break;
         }
      }
      if (nhit < leftHits->mult)         /* merge with former hit */
      {
         leftHits->in[nhit].t = 
               (leftHits->in[nhit].t * leftHits->in[nhit].dE + tleft * dEleft)
             / (leftHits->in[nhit].dE += dEleft);
         rightHits->in[nhit].t = 
               (rightHits->in[nhit].t * rightHits->in[nhit].dE + tright * dEsum)
             / (rightHits->in[nhit].dE += dEright);
      }
      else if (nhit < MAX_HITS)         /* create new hit */
      {
         leftHits->in[nhit].t = tleft;
         leftHits->in[nhit].dE = dEleft;
         leftHits->mult++;
         rightHits->in[nhit].t = tright;
         rightHits->in[nhit].dE = dEright;
         rightHits->mult++;
      }
      else
      {
         fprintf(stderr,"HDGeant error in hitForwardTOF: ");
         fprintf(stderr,"maximum hit count %d exceeded, quitting!\n",MAX_HITS);
         exit(2);
      }
   }

   /* post the hit to the truth tree, once per primary track */
   {
      int mark = (track << 16);
      void** twig = getTwig(&forwardTOFTree, mark);
      if (twig == 0)
      {
         s_ForwardTOF_t* tof = *twig = make_s_ForwardTOF();
         s_TofPoints_t* points = make_s_TofPoints(1);
         points->in[0].x = x[0];
         points->in[0].y = x[1];
         points->mult = 1;
      }
   }
}

/* entry point from fortran */

void hitforwardtof_ (float* xin, float* xout,
                     float* pin, float* pout, float* dEsum, int* track)
{
   hitForwardTOF(xin,xout,pin,pout,*dEsum,*track);
}
