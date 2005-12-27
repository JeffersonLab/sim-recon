/*
 * hitBCal - registers hits for barrel calorimeter
 *
 *	This is a part of the hits package for the
 *	HDGeant simulation program for Hall D.
 *
 *	version 1.0 	-Richard Jones July 16, 2001
 *
 * Programmer's Notes:
 * -------------------
 * 1) In applying the attenuation to light propagating down to both ends
 *    of the modules, there has to be some point where the attenuation
 *    factor is 1.  I chose it to be the midplane, so that in the middle
 *    of the module both ends see the unattenuated E values.  Closer to
 *    either end, that end has a larger E value and the opposite end a
 *    lower E value than the actual deposition.
 * 2) In applying the propagation delay to light propagating down to the
 *    ends of the modules, there has to be some point where the timing
 *    offset is 0.  I chose it to be the midplane, so that for hits in
 *    the middle of the module the t values measure time-of-flight from
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

#define ATTEN_LENGTH	150.
#define C_EFFECTIVE	15.
#define THRESH_MEV	10.
#define TWO_HIT_RESOL	50.
#define MAX_HITS 	100

binTree_t* barrelEMcalTree = 0;
static int cellCount = 0;
static int showerCount = 0;


/* register hits during tracking (from gustep) */

void hitBarrelEMcal (float xin[4], float xout[4],
                     float pin[5], float pout[5], float dEsum,
                     int track, int stack)
{
   float x[3], t;
   float dx[3], dr;
   float xlocal[3];
   float xbcal[3];
   float xHat[] = {1,0,0};

   if (dEsum == 0) return;              /* only seen if it deposits energy */

   x[0] = (xin[0] + xout[0])/2;
   x[1] = (xin[1] + xout[1])/2;
   x[2] = (xin[2] + xout[2])/2;
   t    = (xin[3] + xout[3])/2 * 1e9;
   transformCoord(x,"global",xlocal,"BCAL");
   transformCoord(xHat,"local",xbcal,"BCAL");

   /* post the hit to the hits tree, mark sector as hit */
   {
      int nshot;
      s_BcalUpstreamHits_t* upshots;
      s_BcalDownstreamHits_t* downshots;
      int sector = getsector_();
      int layer  = getlayer_();
      int module = getmodule_();
      float phim = atan2(xbcal[1],xbcal[0]);
      float radius=sqrt(xlocal[0]*xlocal[0] + xlocal[1]*xlocal[1]);
      float dzup = xlocal[2];
      float dzdown = -xlocal[2];
      float tup = t + dzup/C_EFFECTIVE;
      float tdown = t + dzdown/C_EFFECTIVE;
      float Eup = dEsum * exp(-dzup/ATTEN_LENGTH);
      float Edown = dEsum * exp(-dzdown/ATTEN_LENGTH);
      int mark = (module<<16)+ (layer<<9) + sector;
      
      void** twig = getTwig(&barrelEMcalTree, mark);
      if (*twig == 0)
      {
         s_BarrelEMcal_t* bcal = *twig = make_s_BarrelEMcal();
         s_BcalCells_t* cells = make_s_BcalCells(1);
         cells->mult = 1;
         cells->in[0].module = module;
         cells->in[0].layer = layer;
         cells->in[0].sector = sector;
         cells->in[0].bcalUpstreamHits = upshots
                                       = make_s_BcalUpstreamHits(MAX_HITS);
         cells->in[0].bcalDownstreamHits = downshots
                                       = make_s_BcalDownstreamHits(MAX_HITS);
         bcal->bcalCells = cells;
         cellCount++;
      }
      else
      {
         s_BarrelEMcal_t* bcal = *twig;
         upshots = bcal->bcalCells->in[0].bcalUpstreamHits;
         downshots = bcal->bcalCells->in[0].bcalDownstreamHits;
      }

      for (nshot = 0; nshot < upshots->mult; nshot++)
      {
         if (fabs(upshots->in[nshot].t - tup) < TWO_HIT_RESOL)
         {
            break;
         }
      }
      if (nshot < upshots->mult)		/* merge with former hit */
      {
         upshots->in[nshot].t =
                  (upshots->in[nshot].t * upshots->in[nshot].E + tup * Eup)
                / (upshots->in[nshot].E += Eup);
         downshots->in[nshot].t =
             (downshots->in[nshot].t * downshots->in[nshot].E + tdown * Edown)
           / (downshots->in[nshot].E += Edown);
      }
      else if (nshot < MAX_HITS)		/* create new hit */
      {
         upshots->in[nshot].t = tup;
         upshots->in[nshot].E = Eup;
         upshots->mult++;
         downshots->in[nshot].t = tdown;
         downshots->in[nshot].E = Edown;
         downshots->mult++;
      }
      else
      {
         fprintf(stderr,"HDGeant error in hitBarrelEMcal: ");
         fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
      }
   }


   /* post the hit to the truth tree, once per primary track */
   {
      s_BcalTruthShowers_t* showers;
      float r = sqrt(x[0]*x[0]+x[1]*x[1]);
      float phi = atan2(x[1],x[0]);
      int mark = (1<<30) + track;
      void** twig = getTwig(&barrelEMcalTree, mark);
      if (*twig == 0) 
      {
         s_BarrelEMcal_t* bcal = *twig = make_s_BarrelEMcal();
         bcal->bcalTruthShowers = showers = make_s_BcalTruthShowers(1);
         showers->in[0].primary = (stack == 0);
         showers->in[0].track = track;
         showers->in[0].z = x[2];
         showers->in[0].r = r;
         showers->in[0].phi = phi;
         showers->in[0].t = t;
         showers->in[0].E = dEsum;
         showers->mult = 1;
         showerCount++;
      }
      else
      {
         showers = ((s_BarrelEMcal_t*) *twig)->bcalTruthShowers;
         showers->in[0].z =
                        (showers->in[0].z * showers->in[0].E + x[2]*dEsum)
                      / (showers->in[0].E + dEsum);
         showers->in[0].r =
                        (showers->in[0].r * showers->in[0].E + r*dEsum)
                      / (showers->in[0].E + dEsum);
         showers->in[0].phi =
                            (showers->in[0].phi * showers->in[0].E + phi*dEsum)
                          / (showers->in[0].E + dEsum);
         showers->in[0].t = (showers->in[0].t * showers->in[0].E + t*dEsum)
                          / (showers->in[0].E += dEsum);
      }
   }
}

/* entry point from fortran */

void hitbarrelemcal_(float* xin, float* xout,
                     float* pin, float* pout, float* dEsum,
                     int* track, int* stack)
{
   hitBarrelEMcal(xin,xout,pin,pout,*dEsum,*track,*stack);
}




/* pick and package the hits for shipping */

s_BarrelEMcal_t* pickBarrelEMcal ()
{
   s_BarrelEMcal_t* box;
   s_BarrelEMcal_t* item;

   if ((cellCount == 0) && (showerCount == 0))
   {
      return HDDM_NULL;
   }

   box = make_s_BarrelEMcal();
   box->bcalCells = make_s_BcalCells(cellCount);
   box->bcalTruthShowers = make_s_BcalTruthShowers(showerCount);
   while (item = (s_BarrelEMcal_t*) pickTwig(&barrelEMcalTree))
   {
      s_BcalCells_t* cells = item->bcalCells;
      s_BcalTruthShowers_t* showers = item->bcalTruthShowers;

      if (cells != HDDM_NULL)
      {
	 int m = box->bcalCells->mult;
         int mok = 0;

         /* compress out the hits below threshold */
         s_BcalUpstreamHits_t* upshots = cells->in[0].bcalUpstreamHits;
         s_BcalDownstreamHits_t* downshots = cells->in[0].bcalDownstreamHits;
          
         if (upshots != HDDM_NULL)
         {         
            int i,iok;
            for (iok=i=0; i < upshots->mult; i++)
            {
               if (upshots->in[i].E >= THRESH_MEV/1e3)
               {
                  if (iok < i)
                  {
                     upshots->in[iok] = upshots->in[i];
                  }
                  ++iok;
                  ++mok;
               }
            }
            upshots->mult = iok;
            if (iok == 0)
            {
               cells->in[0].bcalUpstreamHits = HDDM_NULL;
               FREE(upshots);
            }
         }
         if (downshots != HDDM_NULL)
         {
            int i,iok;
            for (iok=i=0; i < downshots->mult; i++)
            {
               if (downshots->in[i].E >= THRESH_MEV/1e3)
               {
                  if (iok < i)
                  {
                     downshots->in[iok] = downshots->in[i];
                  }
                  ++iok;
                  ++mok;
               }
            }
            downshots->mult = iok;
            if (iok == 0)
            {
               cells->in[0].bcalDownstreamHits = HDDM_NULL;
               FREE(downshots);
            }
         }

         if (mok)
         {
            box->bcalCells->in[m] = cells->in[0];
            box->bcalCells->mult++;
         }
         FREE(cells);
      }

      else if (showers != HDDM_NULL)
      {
         int m = box->bcalTruthShowers->mult++;
         box->bcalTruthShowers->in[m] = showers->in[0];
         FREE(showers);
      }

      FREE(item);
   }

   cellCount = showerCount = 0;
   return box;
}
