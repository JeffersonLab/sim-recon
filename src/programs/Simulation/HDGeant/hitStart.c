/*
 * hitStart - registers hits for Start counter
 *
 *	This is a part of the hits package for the
 *	HDGeant simulation program for Hall D.
 *
 *	version 1.0 	-Richard Jones July 16, 2001
 *
 * changes: Wed Jun 20 13:19:56 EDT 2007 B. Zihlmann 
 *          add ipart to the function hitStartCntr
 *
 * Programmer's Notes:
 * -------------------
 * 1) In applying the attenuation to light propagating down to the end
 *    of the counters, there has to be some point where the attenuation
 *    factor is 1.  I chose it to be the midplane, so that in the middle
 *    of the counters the attenuation factor is 1.
 * 2) In applying the propagation delay to light propagating down to the
 *    end of the counters, there has to be some point where the timing
 *    offset is 0.  I chose it to be the midplane, so that for hits in
 *    the middle of the counter the t values measure time-of-flight from
 *    the t=0 of the event.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <hddm_s.h>
#include <geant3.h>
#include <bintree.h>
#include "calibDB.h"

static float ATTEN_LENGTH    = 150.;
static float C_EFFECTIVE     = 15.;
static float TWO_HIT_RESOL   = 25.;
static int   MAX_HITS 	     = 100;
static float THRESH_MEV      = 0.150;
static float LIGHT_GUIDE     = 140.;
static float ANGLE_COR       = 1.17;
static float BENT_REGION     = 50.;

binTree_t* startCntrTree = 0;
static int paddleCount = 0;
static int pointCount = 0;
static int initialized = 0;


/* register hits during tracking (from gustep) */

void hitStartCntr (float xin[4], float xout[4],
                   float pin[5], float pout[5], float dEsum,
                   int track, int stack, int history, int ipart)
{
   float x[3], t;
   float dx[3], dr;
   float dEdx;
   float xlocal[3];
   float xvrtx[3];

   if (!initialized) {

    mystr_t strings[50];
    float values[50];
    int nvalues = 50;
    int status = GetConstants("START_COUNTER/start_parms", &nvalues, values, strings);

    if (!status) {
      int ncounter = 0;
      int i;
      for ( i=0;i<(int)nvalues;i++){
        //printf("%d %s \n",i,strings[i].str);
        if (!strcmp(strings[i].str,"START_ATTEN_LENGTH")) {
          ATTEN_LENGTH  = values[i];
          ncounter++;
        }
        if (!strcmp(strings[i].str,"START_C_EFFECTIVE")) {
          C_EFFECTIVE  = values[i];
          ncounter++;
        }
        if (!strcmp(strings[i].str,"START_TWO_HIT_RESOL")) {
          TWO_HIT_RESOL  = values[i];
          ncounter++;
        }
        if (!strcmp(strings[i].str,"START_MAX_HITS")) {
          MAX_HITS  = (int)values[i];
          ncounter++;
        }
        if (!strcmp(strings[i].str,"START_THRESH_MEV")) {
          TWO_HIT_RESOL  = values[i];
          ncounter++;
        }
        if (!strcmp(strings[i].str,"START_LIGHT_GUIDE")) {
          LIGHT_GUIDE  = values[i];
          ncounter++;
        }
        if (!strcmp(strings[i].str,"START_ANGLE_COR")) {
          ANGLE_COR  = values[i];
          ncounter++;
        }
        if (!strcmp(strings[i].str,"START_BENT_REGION")) {
          BENT_REGION  = values[i];
          ncounter++;
        }	
      }
      if (ncounter==8){
        printf("START: ALL parameters loaded from Data Base\n");
      } else if (ncounter<8){
        printf("START: NOT ALL necessary parameters found in Data Base %d out of 8\n",ncounter);
      } else {
        printf("START: SOME parameters found more than once in Data Base\n");
      }
    }
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

   float dbent  = 0.0;
   float dpath  = 0.0;
   if(xlocal[2] >= BENT_REGION){
     dbent = ( xlocal[2] - BENT_REGION )*ANGLE_COR;
     dpath = LIGHT_GUIDE + BENT_REGION + dbent;
   } else {
     dpath  = LIGHT_GUIDE + xlocal[2];
   }

   float dEcorr = dEsum * exp(-dpath/ATTEN_LENGTH);
   float tcorr  = t + dpath/C_EFFECTIVE;


   //   printf("x_gl, z_gl, x_l, z_l %f %f %f\n",
   //  	  xin[0],xin[1],xin[2]);

   //   printf("x_gl, z_gl, x_l, z_l %f %f %f %f %f %f %f\n",
   //  	  x[0],x[1],x[2], xlocal[0],xlocal[1],xlocal[2],dpath);


   /* post the hit to the truth tree */

   if (history == 0)
   {
      int mark = (1<<30) + pointCount;
      void** twig = getTwig(&startCntrTree, mark);
      if (*twig == 0)
      {
         s_StartCntr_t* stc = *twig = make_s_StartCntr();
         s_StcTruthPoints_t* points = make_s_StcTruthPoints(1);
         stc->stcTruthPoints = points;
         points->in[0].primary = (stack == 0);
         points->in[0].track = track;
         points->in[0].t = t;
         points->in[0].z = x[2];
         points->in[0].r = sqrt(x[0]*x[0]+x[1]*x[1]);
         points->in[0].phi = atan2(x[1],x[0]);
         points->in[0].px = pin[0]*pin[4];
         points->in[0].py = pin[1]*pin[4];
         points->in[0].pz = pin[2]*pin[4];
         points->in[0].E = pin[3];
         points->in[0].dEdx = dEcorr;
         points->in[0].ptype = ipart;
         points->in[0].sector = getsector_();
         points->mult = 1;
         pointCount++;


      }
   }

   /* post the hit to the hits tree, mark sector as hit */

   //   if( (ipart==8) && (x[2]<90.)){
   //     printf("x_gl, z_gl, x_l, z_l %f %f %f %f %f %f\n",
   //	    x[0],x[1],x[2], xlocal[0],xlocal[1],xlocal[2]);
   //   }


   if (dEsum > 0)
   {
      int nhit;
      s_StcTruthHits_t* hits;
      int sector = getsector_();
      float phim = atan2(xvrtx[1],xvrtx[0]);


      
      //      printf("x_gl, z_gl, x_l, z_l %f %f %f %f %f %f\n",
      //  x[0],x[1],x[2], xlocal[0],xlocal[1],xlocal[2]);



      //      float dpath = xlocal[2]+(10.2-xlocal[0])*0.4;
      //      float tcorr = t + dpath/C_EFFECTIVE;
      //      float dEcorr = dEsum * exp(-dpath/ATTEN_LENGTH);
      int mark = sector;
      void** twig = getTwig(&startCntrTree, mark);
      if (*twig == 0)
      {
         s_StartCntr_t* stc = *twig = make_s_StartCntr();
         s_StcPaddles_t* paddles = make_s_StcPaddles(1);
         paddles->mult = 1;
         paddles->in[0].sector = sector;
         paddles->in[0].stcTruthHits = hits = make_s_StcTruthHits(MAX_HITS);
         stc->stcPaddles = paddles;
         paddleCount++;
      }
      else
      {
         s_StartCntr_t* stc = *twig;
         hits = stc->stcPaddles->in[0].stcTruthHits;
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
         fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
         exit(2);
      }
   }
}

/* entry point from fortran */

void hitstartcntr_(float* xin, float* xout,
                   float* pin, float* pout, float* dEsum,
                   int* track, int* stack, int* history, int* ipart)
{
   hitStartCntr(xin,xout,pin,pout,*dEsum,*track,*stack,*history,*ipart);
}


/* pick and package the hits for shipping */

s_StartCntr_t* pickStartCntr ()
{
   s_StartCntr_t* box;
   s_StartCntr_t* item;

   if ((paddleCount == 0) && (pointCount == 0))
   {
      return HDDM_NULL;
   }

   box = make_s_StartCntr();
   box->stcPaddles = make_s_StcPaddles(paddleCount);
   box->stcTruthPoints = make_s_StcTruthPoints(pointCount);
   while (item = (s_StartCntr_t*) pickTwig(&startCntrTree))
   {
      s_StcPaddles_t* paddles = item->stcPaddles;
      int paddle;
      s_StcTruthPoints_t* points = item->stcTruthPoints;
      int point;

      for (paddle=0; paddle < paddles->mult; ++paddle)
      {
         int m = box->stcPaddles->mult;

         s_StcTruthHits_t* hits = paddles->in[paddle].stcTruthHits;

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
            box->stcPaddles->in[m] = paddles->in[paddle];
            box->stcPaddles->mult++;
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

      for (point=0; point < points->mult; ++point)
      {
         int m = box->stcTruthPoints->mult++;
         box->stcTruthPoints->in[m] = item->stcTruthPoints->in[point];
      }
      if (points != HDDM_NULL)
      {
         FREE(points);
      }
      FREE(item);
   }

   paddleCount = pointCount = 0;

   if ((box->stcPaddles != HDDM_NULL) &&
       (box->stcPaddles->mult == 0))
   {
      FREE(box->stcPaddles);
      box->stcPaddles = HDDM_NULL;
   }
   if ((box->stcTruthPoints != HDDM_NULL) &&
       (box->stcTruthPoints->mult == 0))
   {
      FREE(box->stcTruthPoints);
      box->stcTruthPoints = HDDM_NULL;
   }
   if ((box->stcPaddles->mult == 0) &&
       (box->stcTruthPoints->mult == 0))
   {
      FREE(box);
      box = HDDM_NULL;
   }
   return box;
}
