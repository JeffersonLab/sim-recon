/*
 * hitStart - registers hits for Start counter
 *
 *    This is a part of the hits package for the
 *    HDGeant simulation program for Hall D.
 *
 *    version 1.0     -Richard Jones July 16, 2001
 *
 * changes: Wed Jun 20 13:19:56 EDT 2007 B. Zihlmann 
 *          add ipart to the function hitStartCntr
 * changes: Tue Aug 25 17:49:21 EDT 2015 E. Pooser
 *          1) Change ANGLE_COR from 1.038 to 1.054 (this corresponds to the
 *          correct 18.5 deg bend towards the beam line in the nose region)
 *          2) Add channel by channel corrections for the propagation time and
 *          attenuation in which constants were determined from beam data and 
 *          bench data (taken at FIU) respectively
 *          
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

#include <HDDM/hddm_s.h>
#include <geant3.h>
#include <bintree.h>
#include <gid_map.h>
#include "calibDB.h"
extern s_HDDM_t* thisInputEvent;

static float ATTEN_LENGTH    = 150.;
static float C_EFFECTIVE     = 15.;
static float TWO_HIT_RESOL   = 25.;
static int   START_MAX_HITS  = 100;
static float THRESH_MEV      = 0.150;
static float LIGHT_GUIDE     = 0.;
//static float ANGLE_COR       = 1.038;
static float ANGLE_COR       = 1.054;
static float BENT_REGION     = 39.465;
static float STRAIGHT_LENGTH = 39.465;
static float BEND_LENGTH     = 3.592375;
static float NOSE_LENGTH     = 15.536625;

static float SC_STRAIGHT_ATTENUATION_A[30], SC_STRAIGHT_ATTENUATION_B[30], SC_STRAIGHT_ATTENUATION_C[30];
static float SC_BENDNOSE_ATTENUATION_A[30], SC_BENDNOSE_ATTENUATION_B[30], SC_BENDNOSE_ATTENUATION_C[30];
static float SC_STRAIGHT_PROPAGATION_A[30], SC_STRAIGHT_PROPAGATION_B[30];
static float SC_BEND_PROPAGATION_A[30],     SC_BEND_PROPAGATION_B[30];
static float SC_NOSE_PROPAGATION_A[30],     SC_NOSE_PROPAGATION_B[30];

static int   NCHANNELS = 30;

// Comment by RTJ:
// When I introduced the convenience constant MAX_HITS,
// I never intended it to be a tunable simulation parameter.
// Do not use it as such.  Do NOT MODIFY it, or the way
// it functions in the algorithm.  If you want to truncate
// the hit list, do it in mcsmear.
#define MAX_HITS 100
 
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

   if (!initialized) {

    mystr_t strings[50];
    float values[50];
    int nvalues = 50;
    int status = GetConstants("START_COUNTER/start_parms", &nvalues, values, strings);

    if (!status) {
      int ncounter = 0;
      int i;
      for ( i=0;i<(int)nvalues;i++){
        //printf("%d %s \n", i, strings[i].str);
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
          START_MAX_HITS  = (int)values[i];
          ncounter++;
        }
        if (!strcmp(strings[i].str,"START_THRESH_MEV")) {
          THRESH_MEV  = values[i];
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

    // Attenuations correction constants for straight section
    int sc_straight_attenuation_a = GetColumn("START_COUNTER/attenuation_factor", &NCHANNELS, SC_STRAIGHT_ATTENUATION_A, "SC_STRAIGHT_ATTENUATION_A");
    if (sc_straight_attenuation_a) 
      printf("ERROR LOADING SC_STRAIGHT_ATTENUATION_A from START_COUNTER/attenuation_factor");
    int sc_straight_attenuation_b = GetColumn("START_COUNTER/attenuation_factor", &NCHANNELS, SC_STRAIGHT_ATTENUATION_B, "SC_STRAIGHT_ATTENUATION_B");
    if (sc_straight_attenuation_b) 
      printf("ERROR LOADING SC_STRAIGHT_ATTENUATION_B from START_COUNTER/attenuation_factor");
    int sc_straight_attenuation_c = GetColumn("START_COUNTER/attenuation_factor", &NCHANNELS, SC_STRAIGHT_ATTENUATION_C, "SC_STRAIGHT_ATTENUATION_C");
    if (sc_straight_attenuation_c) 
      printf("ERROR LOADING SC_STRAIGHT_ATTENUATION_C from START_COUNTER/attenuation_factor");

    // Attenuation correction constants for bend/nose section
    int sc_bendnose_attenuation_a = GetColumn("START_COUNTER/attenuation_factor", &NCHANNELS, SC_BENDNOSE_ATTENUATION_A, "SC_BENDNOSE_ATTENUATION_A");
    if (sc_bendnose_attenuation_a) 
      printf("ERROR LOADING SC_BENDNOSE_ATTENUATION_A from START_COUNTER/attenuation_factor");
    int sc_bendnose_attenuation_b = GetColumn("START_COUNTER/attenuation_factor", &NCHANNELS, SC_BENDNOSE_ATTENUATION_B, "SC_BENDNOSE_ATTENUATION_B");
    if (sc_bendnose_attenuation_b) 
      printf("ERROR LOADING SC_BENDNOSE_ATTENUATION_B from START_COUNTER/attenuation_factor");
    int sc_bendnose_attenuation_c = GetColumn("START_COUNTER/attenuation_factor", &NCHANNELS, SC_BENDNOSE_ATTENUATION_C, "SC_BENDNOSE_ATTENUATION_C");
    if (sc_bendnose_attenuation_c) 
      printf("ERROR LOADING SC_BENDNOSE_ATTENUATION_C from START_COUNTER/attenuation_factor");

    // Propagation time correction constants for straight section
    int sc_straight_propagation_a = GetColumn("START_COUNTER/propagation_speed", &NCHANNELS, SC_STRAIGHT_PROPAGATION_A, "SC_STRAIGHT_PROPAGATION_A");
    if (sc_straight_propagation_a) 
      printf("ERROR LOADING SC_STRAIGHT_PROPAGATION_A from START_COUNTER/propagation_speed");
    int sc_straight_propagation_b = GetColumn("START_COUNTER/propagation_speed", &NCHANNELS, SC_STRAIGHT_PROPAGATION_B, "SC_STRAIGHT_PROPAGATION_B");
    if (sc_straight_propagation_b) 
      printf("ERROR LOADING SC_STRAIGHT_PROPAGATION_B from START_COUNTER/propagation_speed");

    // Propagation time correction constants for bend section
    int sc_bend_propagation_a = GetColumn("START_COUNTER/propagation_speed", &NCHANNELS, SC_BEND_PROPAGATION_A, "SC_BEND_PROPAGATION_A");
    if (sc_bend_propagation_a) 
      printf("ERROR LOADING SC_BEND_PROPAGATION_A from START_COUNTER/propagation_speed");
    int sc_bend_propagation_b = GetColumn("START_COUNTER/propagation_speed", &NCHANNELS, SC_BEND_PROPAGATION_B, "SC_BEND_PROPAGATION_B");
    if (sc_bend_propagation_b) 
      printf("ERROR LOADING SC_BEND_PROPAGATION_B from START_COUNTER/propagation_speed");

    // Propagation time correction constants for nose section
    int sc_nose_propagation_a = GetColumn("START_COUNTER/propagation_speed", &NCHANNELS, SC_NOSE_PROPAGATION_A, "SC_NOSE_PROPAGATION_A");
    if (sc_nose_propagation_a) 
      printf("ERROR LOADING SC_NOSE_PROPAGATION_A from START_COUNTER/propagation_speed");
    int sc_nose_propagation_b = GetColumn("START_COUNTER/propagation_speed", &NCHANNELS, SC_NOSE_PROPAGATION_B, "SC_NOSE_PROPAGATION_B");
    if (sc_nose_propagation_b) 
      printf("ERROR LOADING SC_NOSE_PROPAGATION_B from START_COUNTER/propagation_speed");

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

   /* float dbent  = 0.0; */
   /* float dpath  = 0.0; */
   /* if(xlocal[2] >= BENT_REGION){ */
   /*   dbent = ( xlocal[2] - BENT_REGION )*ANGLE_COR; */
   /*   dpath = BENT_REGION + dbent; */
   /* } else { */
   /*   dpath  = xlocal[2]; */
   /* } */

   /* float dEcorr = dEsum * exp(-dpath/ATTEN_LENGTH); */
   /* float tcorr  = t + dpath/C_EFFECTIVE; */


   //   printf("x_gl, z_gl, x_l, z_l %f %f %f\n",
   //        xin[0],xin[1],xin[2]);

   //   printf("x_gl, z_gl, x_l, z_l %f %f %f %f %f %f %f\n",
   //        x[0],x[1],x[2], xlocal[0],xlocal[1],xlocal[2],dpath);


   /* post the hit to the truth tree */

   int itrack = (stack == 0)? gidGetId(track) : -1;

   if (history == 0)
   {
      int mark = (1<<30) + pointCount;
      void** twig = getTwig(&startCntrTree, mark);
      if (*twig == 0)
      {
         s_StartCntr_t* stc = *twig = make_s_StartCntr();
         s_StcTruthPoints_t* points = make_s_StcTruthPoints(1);
         stc->stcTruthPoints = points;
         int a = thisInputEvent->physicsEvents->in[0].reactions->in[0].vertices->in[0].products->mult;
         points->in[0].primary = (track <= a && stack == 0);
         points->in[0].track = track;
         points->in[0].t = t;
         points->in[0].z = x[2];
         points->in[0].r = sqrt(x[0]*x[0]+x[1]*x[1]);
         points->in[0].phi = atan2(x[1],x[0]);
         points->in[0].px = pin[0]*pin[4];
         points->in[0].py = pin[1]*pin[4];
         points->in[0].pz = pin[2]*pin[4];
         points->in[0].E = pin[3];
         points->in[0].dEdx = dEdx;
         points->in[0].ptype = ipart;
         points->in[0].sector = getsector_wrapper_();
         points->in[0].trackID = make_s_TrackID();
         points->in[0].trackID->itrack = itrack;
         points->mult = 1;
         pointCount++;
      }
   }

   /* post the hit to the hits tree, mark sector as hit */

   //   if( (ipart==8) && (x[2]<90.)){
   //     printf("x_gl, z_gl, x_l, z_l %f %f %f %f %f %f\n",
   //        x[0],x[1],x[2], xlocal[0],xlocal[1],xlocal[2]);
   //   }


   if (dEsum > 0)
   {
      int nhit;
      s_StcTruthHits_t* hits;
      int sector = getsector_wrapper_();
       
      //      printf("x_gl, z_gl, x_l, z_l %f %f %f %f %f %f\n",
      //  x[0],x[1],x[2], xlocal[0],xlocal[1],xlocal[2]);
      
      float dbent  = 0.0;
      float dpath  = 0.0;
      if(xlocal[2] >= BENT_REGION){
          dbent = ( xlocal[2] - BENT_REGION )*ANGLE_COR;
          dpath = BENT_REGION + dbent;
      } else {
          dpath  = xlocal[2];
      }

      /* float dEcorr = dEsum * exp(-dpath/ATTEN_LENGTH); */
      /* float tcorr  = t + dpath/C_EFFECTIVE; */

      /* printf("\n Sector %d Fired \n t = %.5f \n dEsum = %.5f \n dpath = %.5f \n",  */
      /*          sector, t, dEsum, dpath); */

      int sector_index = sector - 1;
      float dEcorr = 9.9E+9;
      float tcorr  = 9.9E+9;
      
      if (xlocal[2] <= STRAIGHT_LENGTH)
    {
      dEcorr = dEsum * exp(dpath*SC_STRAIGHT_ATTENUATION_B[sector_index]);
      tcorr  = t + dpath * SC_STRAIGHT_PROPAGATION_B[sector_index] + SC_STRAIGHT_PROPAGATION_A[sector_index];

      /* printf("HIT OCCURED IN STRAIGHT SECTION \n"); */
      /* printf("Attenuation Corrections: A = %.5f, B = %.5f, C = %.5f \n", SC_STRAIGHT_ATTENUATION_A[sector_index], SC_STRAIGHT_ATTENUATION_B[sector_index], SC_STRAIGHT_ATTENUATION_C[sector_index]); */
      /* printf("Time Corrections: B = %.5f, A = %.5f \n", SC_STRAIGHT_PROPAGATION_B[sector_index], SC_STRAIGHT_PROPAGATION_A[sector_index]);  */
    }
      else if (xlocal[2] > STRAIGHT_LENGTH && xlocal[2] <= (STRAIGHT_LENGTH + BEND_LENGTH))
    {
      dEcorr = dEsum * ((SC_BENDNOSE_ATTENUATION_A[sector_index] * exp(dpath*SC_BENDNOSE_ATTENUATION_B[sector_index]) + SC_BENDNOSE_ATTENUATION_C[sector_index]) / 
                  SC_STRAIGHT_ATTENUATION_A[sector_index]);
      tcorr  = t + dpath * SC_BEND_PROPAGATION_B[sector_index] + SC_BEND_PROPAGATION_A[sector_index];

      /* printf("HIT OCCURED IN BEND SECTION \n"); */
      /* printf("Attenuation Corrections: A = %.5f, B = %.5f, C = %.5f \n", SC_BENDNOSE_ATTENUATION_A[sector_index], SC_BENDNOSE_ATTENUATION_B[sector_index], SC_BENDNOSE_ATTENUATION_C[sector_index]); */
      /* printf("Time Corrections: B = %.5f,  A = %.5f \n", SC_BEND_PROPAGATION_B[sector_index], SC_BEND_PROPAGATION_A[sector_index]);  */
    }
      else if (xlocal[2] > (STRAIGHT_LENGTH + BEND_LENGTH) && xlocal[2] <= (STRAIGHT_LENGTH + BEND_LENGTH + NOSE_LENGTH))
    {
      dEcorr = dEsum * ((SC_BENDNOSE_ATTENUATION_A[sector_index] * exp(dpath*SC_BENDNOSE_ATTENUATION_B[sector_index]) + SC_BENDNOSE_ATTENUATION_C[sector_index]) / 
                  SC_STRAIGHT_ATTENUATION_A[sector_index]);
      
      tcorr  = t + dpath * SC_NOSE_PROPAGATION_B[sector_index] + SC_NOSE_PROPAGATION_A[sector_index];

      /* printf("HIT OCCURED IN NOSE SECTION \n"); */
      /* printf("Attenuation Corrections: A = %.5f, B = %.5f, C = %.5f \n", SC_BENDNOSE_ATTENUATION_A[sector_index], SC_BENDNOSE_ATTENUATION_B[sector_index], SC_BENDNOSE_ATTENUATION_C[sector_index]); */
      /* printf("Time Corrections: B = %.5f,  A = %.5f \n", SC_NOSE_PROPAGATION_B[sector_index], SC_NOSE_PROPAGATION_A[sector_index]); */ 
    }
      else return;

      /* printf("tcorr = %.5f \n dEcorr = %.5f \n", tcorr, dEcorr); */
    
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
      if (nhit < hits->mult)        /* merge with former hit */
      {
         if (tcorr < hits->in[nhit].t)
         {
            hits->in[nhit].ptype = ipart;
            hits->in[nhit].itrack = itrack;
         }
         hits->in[nhit].t = 
                 (hits->in[nhit].t * hits->in[nhit].dE + tcorr * dEcorr) /
                 (hits->in[nhit].dE + dEcorr);
            hits->in[nhit].dE += dEcorr;
      }
      else if (nhit < MAX_HITS)        /* create new hit */
      {
         hits->in[nhit].t = tcorr ;
         hits->in[nhit].dE = dEcorr;
         hits->in[nhit].ptype = ipart;
         hits->in[nhit].itrack = itrack;
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
   while ((item = (s_StartCntr_t*) pickTwig(&startCntrTree)))
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

      int last_track = -1;
      double last_t = 1e9;
      for (point=0; point < points->mult; ++point)
      {
         if (points->in[point].trackID->itrack > 0 &&
            (points->in[point].track != last_track ||
             fabs(points->in[point].t - last_t) > 0.1))
         {
            int m = box->stcTruthPoints->mult++;
            box->stcTruthPoints->in[m] = item->stcTruthPoints->in[point];
            last_track = points->in[point].track;
            last_t = points->in[point].t;
         }
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
