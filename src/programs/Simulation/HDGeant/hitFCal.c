/*
 * hitFCal - registers hits for forward calorimeter
 *
 *	This is a part of the hits package for the
 *	HDGeant simulation program for Hall D.
 *
 *	version 1.0 	-Richard Jones July 16, 2001
 *
 * changes: Wed Jun 20 13:19:56 EDT 2007 B. Zihlmann 
 *          add ipart to the function hitForwardEMcal
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#include <hddm_s.h>
#include <geant3.h>
#include <bintree.h>

#include "calibDB.h"


static float ATTEN_LENGTH1    =	100.;   //Outer Glass
static float  ATTEN_LENGTH2   =  73.1;  //Inner Glass
static float  RHG_RADIUS      =  30.;   //Radius of inner glass
static float  C_EFFECTIVE     =  15.;   // cm/ns
static float  WIDTH_OF_BLOCK  =   4.;   //cm
static float  LENGTH_OF_BLOCK =  45.;   //cm
static float  TWO_HIT_RESOL   =  75.;   // ns
static int    MAX_HITS        =  100;   // maximum hits per block
static float  THRESH_MEV      =   5.;
static float  ACTIVE_RADIUS   = 120.;
static int    CENTRAL_ROW     =  29;
static int    CENTRAL_COLUMN  =  29;


binTree_t* forwardEMcalTree = 0;
static int blockCount = 0;
static int showerCount = 0;
static int initialized = 0;


/* register hits during tracking (from gustep) */

void hitForwardEMcal (float xin[4], float xout[4],
                      float pin[5], float pout[5], float dEsum,
                      int track, int stack, int history, int ipart)
{
   float x[3], t;
   float xfcal[3];
   float zeroHat[] = {0,0,0};

   if (!initialized){
     
     mystr_t strings[50];
     float values[50];
     int nvalues = 50;
     int status = GetConstants("FCAL/fcal_parms", &nvalues, values, strings);
     
     if (!status) {
       
       int ncounter = 0;
       int i;
       for ( i=0;i<(int)nvalues;i++){
	 //printf("%d %s \n",i,strings[i].str);
	 if (!strcmp(strings[i].str,"FCAL_ATTEN_LENGTH1")) {
	   ATTEN_LENGTH1  = values[i];
	   ncounter++;
	 }
         if (!strcmp(strings[i].str,"FCAL_ATTEN_LENGTH2")) {
	   ATTEN_LENGTH2  = values[i];
	   ncounter++;
	 }
	 if (!strcmp(strings[i].str,"FCAL_RHG_RADIUS")) {
	   RHG_RADIUS  = values[i];
	   ncounter++;
	 }
	 if (!strcmp(strings[i].str,"FCAL_C_EFFECTIVE")) {
	   C_EFFECTIVE  = values[i];
	   ncounter++;
	 }
	 if (!strcmp(strings[i].str,"FCAL_WIDTH_OF_BLOCK")) {
	   WIDTH_OF_BLOCK  = values[i];
	   ncounter++;
	 }
	 if (!strcmp(strings[i].str,"FCAL_LENGTH_OF_BLOCK")) {
	   LENGTH_OF_BLOCK  = values[i];
	   ncounter++;
	 }
	 if (!strcmp(strings[i].str,"FCAL_TWO_HIT_RESOL")) {
	   TWO_HIT_RESOL  = values[i];
	   ncounter++;
	 }
	 if (!strcmp(strings[i].str,"FCAL_MAX_HITS")) {
	   MAX_HITS  = (int)values[i];
	   ncounter++;
	 }
	 if (!strcmp(strings[i].str,"FCAL_THRESH_MEV")) {
	   THRESH_MEV  = values[i];
	   ncounter++;
	 }
	 if (!strcmp(strings[i].str,"FCAL_ACTIVE_RADIUS")) {
	   ACTIVE_RADIUS  = values[i];
	   ncounter++;
	 }
	 if (!strcmp(strings[i].str,"FCAL_CENTRAL_ROW")) {
	   CENTRAL_ROW  = (int)values[i];
	   ncounter++;
	 }
	 if (!strcmp(strings[i].str,"FCAL_CENTRAL_COLUMN")) {
	   CENTRAL_COLUMN  = (int)values[i];
	   ncounter++;
	 }
       }
       if (ncounter==12){
	 printf("FCAL: ALL parameters loaded from Data Base\n");
       } else if (ncounter<12){
	 printf("FCAL: NOT ALL necessary parameters found in Data Base %d out of 12\n",ncounter);
       } else {
	 printf("FCAL: SOME parameters found more than once in Data Base\n");
       }
     }
     initialized = 1;
   }
   
   x[0] = (xin[0] + xout[0])/2;
   x[1] = (xin[1] + xout[1])/2;
   x[2] = (xin[2] + xout[2])/2;
   t    = (xin[3] + xout[3])/2 * 1e9;
   transformCoord(x,"global",xfcal,"FCAL");

   /* post the hit to the truth tree */

   if ((history == 0) && (pin[3] > THRESH_MEV/1e3))
   {
      s_FcalTruthShowers_t* showers;
      int mark = (1<<30) + showerCount;
      void** twig = getTwig(&forwardEMcalTree, mark);
      if (*twig == 0)
      {
         s_ForwardEMcal_t* cal = *twig = make_s_ForwardEMcal();
         cal->fcalTruthShowers = showers = make_s_FcalTruthShowers(1);
         showers->in[0].primary = (stack == 0);
         showers->in[0].track = track;
         showers->in[0].t = xin[3]*1e9;
         showers->in[0].x = xin[0];
         showers->in[0].y = xin[1];
         showers->in[0].z = xin[2];
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
      s_FcalTruthHits_t* hits;
      int row = getrow_();
      int column = getcolumn_();
      
      float dist = 0.5*LENGTH_OF_BLOCK-xfcal[2];
      float dEcorr;
      float y0 = (row - CENTRAL_ROW)*WIDTH_OF_BLOCK;
      float x0 = (column - CENTRAL_COLUMN)*WIDTH_OF_BLOCK;
      float rad = sqrt(x0*x0+y0*y0);
    

      
      if(rad>RHG_RADIUS){

	dEcorr = dEsum * exp(-dist/ATTEN_LENGTH1);

      }
      
      else{
	
	dEcorr = dEsum * exp(-dist/ATTEN_LENGTH2);
      
      }

      float tcorr = t + dist/C_EFFECTIVE;
      int mark = ((row+1)<<16) + (column+1);
      void** twig = getTwig(&forwardEMcalTree, mark);
      if (*twig == 0)
      {
         s_ForwardEMcal_t* cal = *twig = make_s_ForwardEMcal();
         s_FcalBlocks_t* blocks = make_s_FcalBlocks(1);
         blocks->mult = 1;
         blocks->in[0].row = row;
         blocks->in[0].column = column;
         blocks->in[0].fcalTruthHits = hits = make_s_FcalTruthHits(MAX_HITS);
         cal->fcalBlocks = blocks;
         blockCount++;
      }
      else
      {
         s_ForwardEMcal_t* cal = *twig;
         hits = cal->fcalBlocks->in[0].fcalTruthHits;
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
         fprintf(stderr,"HDGeant error in hitforwardEMcal: ");
         fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
         exit(2);
      }
   }
}

/* entry point from fortran */

void hitforwardemcal_(float* xin, float* xout,
                      float* pin, float* pout, float* dEsum,
                      int* track, int* stack, int* history, int* ipart)
{
   hitForwardEMcal(xin,xout,pin,pout,*dEsum,*track,*stack,*history, *ipart);
}


/* pick and package the hits for shipping */

s_ForwardEMcal_t* pickForwardEMcal ()
{
   s_ForwardEMcal_t* box;
   s_ForwardEMcal_t* item;

#if TESTING_CAL_CONTAINMENT
   double Etotal = 0;
#endif
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
      int block;
      s_FcalTruthShowers_t* showers = item->fcalTruthShowers;
      int shower;
      for (block=0; block < blocks->mult; ++block)
      {
         int row = blocks->in[block].row;
         int column = blocks->in[block].column;
         float y0 = (row - CENTRAL_ROW)*WIDTH_OF_BLOCK;
         float x0 = (column - CENTRAL_COLUMN)*WIDTH_OF_BLOCK;
         float dist = sqrt(x0*x0+y0*y0);

         s_FcalTruthHits_t* hits = blocks->in[block].fcalTruthHits;

         /* compress out the hits outside the active region */
         if (dist < ACTIVE_RADIUS)
         {
	    int m = box->fcalBlocks->mult;

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
               box->fcalBlocks->in[m] = blocks->in[block];
               box->fcalBlocks->mult++;
            }
            else if (hits != HDDM_NULL)
            {
               FREE(hits);
            }
         }
         else if (hits != HDDM_NULL)
         {
            FREE(hits);
         }
      }

      for (shower=0; shower < showers->mult; ++shower)
      {
         int m = box->fcalTruthShowers->mult++;
         box->fcalTruthShowers->in[m] = showers->in[shower];
      }
      if (blocks != HDDM_NULL)
      {
         FREE(blocks);
      }
      if (showers != HDDM_NULL)
      {
         FREE(showers);
      }
      FREE(item);
   }

   blockCount = showerCount = 0;

   if ((box->fcalBlocks != HDDM_NULL) &&
       (box->fcalBlocks->mult == 0))
   {
      FREE(box->fcalBlocks);
      box->fcalBlocks = HDDM_NULL;
   }
   if ((box->fcalTruthShowers != HDDM_NULL) &&
       (box->fcalTruthShowers->mult == 0))
   {
      FREE(box->fcalTruthShowers);
      box->fcalTruthShowers = HDDM_NULL;
   }
   if ((box->fcalBlocks->mult == 0) &&
       (box->fcalTruthShowers->mult == 0))
   {
      FREE(box);
      box = HDDM_NULL;
   }
#if TESTING_CAL_CONTAINMENT
  printf("FCal energy sum: %f\n",Etotal/0.614);
#endif
   return box;
}
