/*
 * hitBCal - registers hits for barrel calorimeter
 *
 *	This is a part of the hits package for the
 *	HDGeant simulation program for Hall D.
 *
 *	version 1.0 	-Richard Jones July 16, 2001
 *
 * changes: Wed Jun 20 13:19:56 EDT 2007 B. Zihlmann 
 *          add ipart to the function hitBarrelEMcal
 * changes: Fri Jul 13 08:54:48 EDT 2007 M. Shepherd
 *          remove attenuation, condense up and dowstream hits
 *          pass up true z position and hit time
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <hddm_s.h>
#include <geant3.h>
#include <bintree.h>

#define THRESH_MEV	     1.
#define TWO_HIT_RESOL	50.
#define MAX_HITS 	   100

binTree_t* barrelEMcalTree = 0;
static int cellCount = 0;
static int showerCount = 0;

/* register hits during tracking (from gustep) */

void hitBarrelEMcal (float xin[4], float xout[4],
                     float pin[5], float pout[5], float dEsum,
                     int track, int stack, int history, int ipart)
{
   float x[3], t;
   float dx[3], dr;
   float xlocal[3];
   float xbcal[3];
   float xHat[] = {1,0,0};

   x[0] = (xin[0] + xout[0])/2;
   x[1] = (xin[1] + xout[1])/2;
   x[2] = (xin[2] + xout[2])/2;
   t    = (xin[3] + xout[3])/2 * 1e9;
   transformCoord(x,"global",xlocal,"BCAL");
   transformCoord(xHat,"local",xbcal,"BCAL");

	/* Under certain conditions the time in xout[3] will
	   be invalid (unusually large). Check for this and
		use only the in time in these cases
	*/
	if(xout[3] > 1.0){
		 t = xin[3] * 1e9;
		printf("xout[3]=%g\n", xout[3]);
	}

   /* post the hit to the truth tree */

   if ((history == 0) && (pin[3] > THRESH_MEV/1e3)) {
      s_BcalTruthShowers_t* showers;
      float r = sqrt(xin[0]*xin[0]+xin[1]*xin[1]);
      float phi = atan2(xin[1],xin[0]);
      int mark = (1<<30) + showerCount;
      void** twig = getTwig(&barrelEMcalTree, mark);
      if (*twig == 0) 
      {
         s_BarrelEMcal_t* bcal = *twig = make_s_BarrelEMcal();
         bcal->bcalTruthShowers = showers = make_s_BcalTruthShowers(1);
         showers->in[0].primary = (stack == 0);
         showers->in[0].track = track;
         showers->in[0].z = xin[2];
         showers->in[0].r = r;
         showers->in[0].phi = phi;
         showers->in[0].t = xin[3]*1e9;
         showers->in[0].px = pin[0]*pin[4];
         showers->in[0].py = pin[1]*pin[4];
         showers->in[0].pz = pin[2]*pin[4];
         showers->in[0].E = pin[3];
         showers->in[0].ptype = ipart;
         showers->mult = 1;
         showerCount++;
      }
   }

   /* post the hit to the hits tree, mark sector as hit */

   if (dEsum > 0)
   {
      int nshot;
      s_BcalHits_t* hits;
      int sector = getsector_();
      int layer  = getlayer_();
      int module = getmodule_();
      float phim = atan2(xbcal[1],xbcal[0]);
      float radius=sqrt(xlocal[0]*xlocal[0] + xlocal[1]*xlocal[1]);
      float zLocal = xlocal[2];
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
         cells->in[0].bcalHits = hits
                               = make_s_BcalHits(MAX_HITS);
         bcal->bcalCells = cells;
         cellCount++;
      }
      else
      {
         s_BarrelEMcal_t* bcal = *twig;
         hits = bcal->bcalCells->in[0].bcalHits;
      }

      for (nshot = 0; nshot < hits->mult; nshot++)
      {
         if (fabs(hits->in[nshot].t - t) < TWO_HIT_RESOL)
         {
            break;
         }
      }
      if (nshot < hits->mult)		/* merge with former hit */
      {
         hits->in[nshot].t =
                  (hits->in[nshot].t * hits->in[nshot].E + t * dEsum)
                / (hits->in[nshot].E + dEsum);
         hits->in[nshot].zLocal =
                  (hits->in[nshot].zLocal * hits->in[nshot].E + zLocal * dEsum)
                / (hits->in[nshot].E + dEsum);
         hits->in[nshot].E += dEsum;
      }
      else if (nshot < MAX_HITS)		/* create new hit */
      {
         hits->in[nshot].t = t;
         hits->in[nshot].E = dEsum;
         hits->in[nshot].zLocal = zLocal;
         hits->mult++;
      }
      else
      {
         fprintf(stderr,"HDGeant error in hitBarrelEMcal: ");
         fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
      }
   }
}

/* entry point from fortran */

void hitbarrelemcal_(float* xin, float* xout,
                     float* pin, float* pout, float* dEsum,
                     int* track, int* stack, int* history, int* ipart)
{
   hitBarrelEMcal(xin,xout,pin,pout,*dEsum,*track,*stack,*history, *ipart);
}




/* pick and package the hits for shipping */

s_BarrelEMcal_t* pickBarrelEMcal ()
{
   s_BarrelEMcal_t* box;
   s_BarrelEMcal_t* item;
#if TESTING_CAL_CONTAINMENT
  double Etotal = 0;
#endif

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
      int cell;
      s_BcalTruthShowers_t* showers = item->bcalTruthShowers;
      int shower;
      for (cell=0; cell < cells->mult; ++cell)
      {
	 int m = box->bcalCells->mult;
         int mok = 0;

         s_BcalHits_t* hits = cells->in[cell].bcalHits;
          
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
               ++mok;
            }
         }
         if (hits != HDDM_NULL)
         {
            hits->mult = iok;
            if (iok == 0)
            {
               cells->in[cell].bcalHits = HDDM_NULL;
               FREE(hits);
            }
         }

         if (mok)
         {
            box->bcalCells->in[m] = cells->in[cell];
            box->bcalCells->mult++;
         }
      }
      if (cells != HDDM_NULL)
      {
         FREE(cells);
      }

      for (shower=0; shower < showers->mult; ++shower)
      {
         int m = box->bcalTruthShowers->mult++;
         box->bcalTruthShowers->in[m] = showers->in[shower];
      }
      if (showers != HDDM_NULL)
      {
         FREE(showers);
      }

      FREE(item);
   }

   cellCount = showerCount = 0;

   if ((box->bcalCells != HDDM_NULL) &&
       (box->bcalCells->mult == 0))
   {
      FREE(box->bcalCells);
      box->bcalCells = HDDM_NULL;
   }
   if ((box->bcalTruthShowers != HDDM_NULL) &&
       (box->bcalTruthShowers->mult == 0))
   {
      FREE(box->bcalTruthShowers);
      box->bcalTruthShowers = HDDM_NULL;
   }
   if ((box->bcalCells->mult == 0) &&
       (box->bcalTruthShowers->mult == 0))
   {
      FREE(box);
      box = HDDM_NULL;
   }
#if TESTING_CAL_CONTAINMENT
  printf("BCal energy sum: %f\n",Etotal);
#endif
   return box;
}
