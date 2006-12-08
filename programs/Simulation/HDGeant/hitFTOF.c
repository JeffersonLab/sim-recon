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
#include <math.h>

#include <hddm_s.h>
#include <geant3.h>
#include <bintree.h>

#define ATTEN_LENGTH	150
#define C_EFFECTIVE	15
#define TWO_HIT_RESOL   25.
#define MAX_HITS        100
#define THRESH_MEV      0.8

binTree_t* forwardTOFTree = 0;
static int counterCount = 0;
static int pointCount = 0;


/* register hits during tracking (from gustep) */

void hitForwardTOF (float xin[4], float xout[4],
                    float pin[5], float pout[5], float dEsum,
                    int track, int stack, int history)
{
   float x[3], t;
   float dx[3], dr;
   float dEdx;
   float xlocal[3];
   float xftof[3];
   float zeroHat[] = {0,0,0};
   int plane = getplane_();

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

   /* post the hit to the truth tree */

   if ((history == 0) && (plane == 0))
   {
      int mark = (1<<30) + pointCount;
      void** twig = getTwig(&forwardTOFTree, mark);
      if (*twig == 0)
      {
         s_ForwardTOF_t* tof = *twig = make_s_ForwardTOF();
         s_FtofTruthPoints_t* points = make_s_FtofTruthPoints(1);
         tof->ftofTruthPoints = points;
         points->in[0].primary = (stack == 0);
         points->in[0].track = track;
         points->in[0].x = x[0];
         points->in[0].y = x[1];
         points->in[0].z = x[2];
         points->in[0].t = t;
         points->mult = 1;
         pointCount++;
      }
   }

   /* post the hit to the hits tree, mark slab as hit */

   if (dEsum > 0)
   {
      int nhit;
      s_FtofNorthHits_t* northHits;
      s_FtofSouthHits_t* southHits;
      int row = getrow_();
      int column = getcolumn_();
      float dxnorth = xlocal[0];
      float dxsouth = -xlocal[0];
      float tnorth  = (column == 2) ? 0 : t + dxnorth/C_EFFECTIVE;
      float tsouth = (column == 1) ? 0 : t + dxsouth/C_EFFECTIVE;
      float dEnorth  = (column == 2) ? 0 : dEsum * exp(-dxnorth/ATTEN_LENGTH);
      float dEsouth = (column == 1) ? 0 : dEsum * exp(-dxsouth/ATTEN_LENGTH);
      float ycenter = (fabs(xftof[1]) < 1e-4) ? 0 : xftof[1];
      int mark = (plane<<20) + (row<<10) + column;
      void** twig = getTwig(&forwardTOFTree, mark);
      if (*twig == 0)
      {
         s_ForwardTOF_t* tof = *twig = make_s_ForwardTOF();
         s_FtofCounters_t* counters = make_s_FtofCounters(1);
         counters->mult = 1;
         counters->in[0].plane = plane;
         counters->in[0].bar = row;
         northHits = HDDM_NULL;
         southHits = HDDM_NULL;
         if (column == 0 || column == 1)
         {
           counters->in[0].ftofNorthHits = northHits
                                        = make_s_FtofNorthHits(MAX_HITS);
         }
         if (column == 0 || column == 2)
         {
           counters->in[0].ftofSouthHits = southHits
                                         = make_s_FtofSouthHits(MAX_HITS);
         }
         tof->ftofCounters = counters;
         counterCount++;
      }
      else
      {
         s_ForwardTOF_t* tof = *twig;
         northHits = tof->ftofCounters->in[0].ftofNorthHits;
         southHits = tof->ftofCounters->in[0].ftofSouthHits;
      }

      if (northHits != HDDM_NULL)
      {
         for (nhit = 0; nhit < northHits->mult; nhit++)
         {
            if (fabs(northHits->in[nhit].t - t) < TWO_HIT_RESOL)
            {
               break;
            }
         }
         if (nhit < northHits->mult)         /* merge with former hit */
         {
            northHits->in[nhit].t = 
               (northHits->in[nhit].t * northHits->in[nhit].dE 
                                      + tnorth * dEnorth)
                / (northHits->in[nhit].dE += dEnorth);
         }
         else if (nhit < MAX_HITS)         /* create new hit */
         {
            northHits->in[nhit].t = tnorth;
            northHits->in[nhit].dE = dEnorth;
            northHits->mult++;
         }
         else
         {
            fprintf(stderr,"HDGeant error in hitForwardTOF: ");
            fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
         }
      }

      if (southHits != HDDM_NULL)
      {
         for (nhit = 0; nhit < southHits->mult; nhit++)
         {
            if (fabs(southHits->in[nhit].t - t) < TWO_HIT_RESOL)
            {
               break;
            }
         }
         if (nhit < southHits->mult)         /* merge with former hit */
         {
            southHits->in[nhit].t = 
             (southHits->in[nhit].t * southHits->in[nhit].dE
             + tsouth * dEsouth) / (southHits->in[nhit].dE += dEsouth);
         }
         else if (nhit < MAX_HITS)         /* create new hit */
         {
            southHits->in[nhit].t = tsouth;
            southHits->in[nhit].dE = dEsouth;
            southHits->mult++;
         }
         else
         {
            fprintf(stderr,"HDGeant error in hitForwardTOF: ");
            fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
         }
      }
   }
}

/* entry point from fortran */

void hitforwardtof_ (float* xin, float* xout,
                     float* pin, float* pout, float* dEsum,
                     int* track, int* stack, int* history)
{
   hitForwardTOF(xin,xout,pin,pout,*dEsum,*track,*stack,*history);
}


/* pick and package the hits for shipping */

s_ForwardTOF_t* pickForwardTOF ()
{
   s_ForwardTOF_t* box;
   s_ForwardTOF_t* item;

   if ((counterCount == 0) && (pointCount == 0))
   {
      return HDDM_NULL;
   }

   box = make_s_ForwardTOF();
   box->ftofCounters = make_s_FtofCounters(counterCount);
   box->ftofTruthPoints = make_s_FtofTruthPoints(pointCount);
   while (item = (s_ForwardTOF_t*) pickTwig(&forwardTOFTree))
   {
      s_FtofCounters_t* counters = item->ftofCounters;
      int counter;
      s_FtofTruthPoints_t* points = item->ftofTruthPoints;
      int point;

      for (counter=0; counter < counters->mult; ++counter)
      {
         s_FtofNorthHits_t* northHits = counters->in[counter].ftofNorthHits;
         s_FtofSouthHits_t* southHits = counters->in[counter].ftofSouthHits;

      /* compress out the hits below threshold */
         int iok,i;
         int mok=0;
         for (iok=i=0; i < northHits->mult; i++)
         {
            if (northHits->in[i].dE >= THRESH_MEV/1e3)
            {
               if (iok < i)
               {
                  northHits->in[iok] = northHits->in[i];
               }
               ++mok;
               ++iok;
            }
         }
         if (iok)
         {
            northHits->mult = iok;
         }
         else if (northHits != HDDM_NULL)
         {
            counters->in[counter].ftofNorthHits = HDDM_NULL;
            FREE(northHits);
         }

         for (iok=i=0; i < southHits->mult; i++)
         {
            if (southHits->in[i].dE >= THRESH_MEV/1e3)
            {
               if (iok < i)
               {
                  southHits->in[iok] = southHits->in[i];
               }
               ++mok;
               ++iok;
            }
         }
         if (iok)
         {
            southHits->mult = iok;
         }
         else if (southHits != HDDM_NULL) 
         {
            counters->in[counter].ftofSouthHits = HDDM_NULL;
            FREE(southHits);
         }
         if (mok)
         {
            int m = box->ftofCounters->mult++;
            box->ftofCounters->in[m] = counters->in[counter];
         }
      }
      if (counters != HDDM_NULL)
      {
         FREE(counters);
      }

      for (point=0; point < points->mult; ++point)
      {
         int m = box->ftofTruthPoints->mult++;
         box->ftofTruthPoints->in[m] = points->in[point];
      }
      if (points != HDDM_NULL)
      {
         FREE(points);
      }
      FREE(item);
   }

   counterCount = pointCount = 0;

   if ((box->ftofCounters != HDDM_NULL) &&
       (box->ftofCounters->mult == 0))
   {
      FREE(box->ftofCounters);
      box->ftofCounters = HDDM_NULL;
   }
   if ((box->ftofTruthPoints != HDDM_NULL) &&
       (box->ftofTruthPoints->mult == 0))
   {
      FREE(box->ftofTruthPoints);
      box->ftofTruthPoints = HDDM_NULL;
   }
   if ((box->ftofCounters->mult == 0) &&
       (box->ftofTruthPoints->mult == 0))
   {
      FREE(box);
      box = HDDM_NULL;
   }
   return box;
}
