/*
 * hitFTOF - registers hits for forward Time-Of-Flight
 *
 *	This is a part of the hits package for the
 *	HDGeant simulation program for Hall D.
 *
 *	version 1.0 	-Richard Jones July 16, 2001
 *
 * Programmer's Notes:
 * -------------------
 * 1) In applying the attenuation to light propagating down to both ends
 *    of the counters, there has to be some point where the attenuation
 *    factor is 1.  I chose it to be the midplane, so that in the middle
 *    of the counter both ends see the unattenuated dE values.  Closer to
 *    either end, that end has a larger dE value and the opposite end a
 *    lower dE value than the actual deposition.
 * 2) In applying the propagation delay to light propagating down to the
 *    ends of the counters, there has to be some point where the timing
 *    offset is 0.  I chose it to be the midplane, so that for hits in
 *    the middle of the counter the t values measure time-of-flight from
 *    the t=0 of the event.  For hits closer to one end, that end sees
 *    a t value smaller than its true time-of-flight, and the other end
 *    sees a value correspondingly larger.  The average is the true tof.
 */

#include <stdlib.h>
#include <stdio.h>

#include <hddm_s.h>
#include <geant3.h>
#include <bintree.h>

#define ATTEN_LENGTH	150
#define C_EFFECTIVE	15
#define TWO_HIT_RESOL   25.
#define MAX_HITS        100

binTree_t* forwardTOFTree = 0;
static int slabCount = 0;
static int pointCount = 0;


/* register hits during tracking (from gustep) */

void hitForwardTOF (float xin[4], float xout[4],
                    float pin[5], float pout[5], float dEsum, int track)
{
   float x[3], t;
   float dx[3], dr;
   float dEdx;
   float xlocal[3];
   float xftof[3];
   float zeroHat[] = {0,0,0};

   if (dEsum == 0) return;              /* only seen if it deposits energy */

   x[0] = (xin[0] + xout[0])/2;
   x[1] = (xin[1] + xout[1])/2;
   x[2] = (xin[2] + xout[2])/2;
   t    = (xin[3] + xout[3])/2 * 1e9;
   transformCoord(x,"global",xlocal,"FTOF");
   transformCoord(zeroHat,"local",xftof,"FTOF");
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
      int row = getrow_();
      int column = getcolumn_();
      float dxleft = xlocal[0];
      float dxright = -xlocal[0];
      float tleft  = (column == 2) ? 0 : t + dxleft/C_EFFECTIVE;
      float tright = (column == 1) ? 0 : t + dxright/C_EFFECTIVE;
      float dEleft  = (column == 2) ? 0 : dEsum * exp(-dxleft/ATTEN_LENGTH);
      float dEright = (column == 1) ? 0 : dEsum * exp(-dxright/ATTEN_LENGTH);
      int mark = (row << 8) + column;
      void** twig = getTwig(&forwardTOFTree, mark);
      if (*twig == 0)
      {
         s_ForwardTOF_t* tof = *twig = make_s_ForwardTOF();
         tof->slabs = make_s_Slabs(1);
         tof->slabs->mult = 1;
         tof->slabs->in[0].y = xftof[1];
         tof->slabs->in[0].left = make_s_Left();
         tof->slabs->in[0].left->hits = leftHits = make_s_Hits(MAX_HITS);
         tof->slabs->in[0].right = make_s_Right();
         tof->slabs->in[0].right->hits = rightHits = make_s_Hits(MAX_HITS);
         slabCount++;
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
      if (*twig == 0)
      {
         s_ForwardTOF_t* tof = *twig = make_s_ForwardTOF();
         s_TofPoints_t* points = make_s_TofPoints(1);
         tof->tofPoints = points;
         points->in[0].x = x[0];
         points->in[0].y = x[1];
         points->in[0].t = t;
         points->mult = 1;
         pointCount++;
      }
   }
}

/* entry point from fortran */

void hitforwardtof_ (float* xin, float* xout,
                     float* pin, float* pout, float* dEsum, int* track)
{
   hitForwardTOF(xin,xout,pin,pout,*dEsum,*track);
}


/* pick and package the hits for shipping */

s_ForwardTOF_t* pickForwardTOF ()
{
   s_ForwardTOF_t* box;
   s_ForwardTOF_t* item;

   if ((slabCount == 0) && (pointCount == 0))
   {
      return 0;
   }

   box = make_s_ForwardTOF();
   box->slabs = make_s_Slabs(slabCount);
   box->tofPoints = make_s_TofPoints(pointCount);
   while (item = (s_ForwardTOF_t*) pickTwig(&forwardTOFTree))
   {
      if (item->slabs)
      {
         int m = box->slabs->mult++;
         box->slabs->in[m] = item->slabs->in[0];
         FREE(item->slabs);
      }
      else if (item->tofPoints)
      {
         int m = box->tofPoints->mult++;
         box->tofPoints->in[m] = item->tofPoints->in[0];
         FREE(item->tofPoints);
      }
      FREE(item);
   }
   slabCount = pointCount = 0;
   return box;
}
