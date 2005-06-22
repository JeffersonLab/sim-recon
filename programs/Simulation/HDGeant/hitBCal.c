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
static int moduleCount = 0;
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
      s_Showers_t* upshots;
      s_Showers_t* downshots;
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
      
      /* checks */
      /*  void** twig = getTwig(&barrelEMcalTree, sector);*/
      /*   printf("\nInside hitbcal, layer = %d  ",layer); */
      /* printf("\nInside hitbcal, sector = %d  ",sector); */
      /* printf("\nInside hitbcal, mark1 = %d  ",mark); */

      void** twig = getTwig(&barrelEMcalTree, mark);

      /*   printf("\nInside hitbcal, layer = %d  ",layer); */

      if (*twig == 0)
      {
         s_BarrelEMcal_t* cal = *twig = make_s_BarrelEMcal();
         cal->mods = make_s_Mods(1);
         cal->mods->mult = 1;
         cal->mods->in[0].module = module;
         cal->mods->in[0].shells  = make_s_Shells(1);
         cal->mods->in[0].shells->mult = 1;
         cal->mods->in[0].shells->in[0].layer = layer;
         cal->mods->in[0].shells->in[0].cones     = make_s_Cones(1); 
         cal->mods->in[0].shells->in[0].cones->mult = 1;
         cal->mods->in[0].shells->in[0].cones->in[0].sector = sector;
         cal->mods->in[0].shells->in[0].cones->in[0].upstream = make_s_Upstream();
         cal->mods->in[0].shells->in[0].cones->in[0].downstream = make_s_Downstream();
         cal->mods->in[0].shells->in[0].cones->in[0].upstream->showers =
         upshots = make_s_Showers(MAX_HITS);
         cal->mods->in[0].shells->in[0].cones->in[0].downstream->showers =
         downshots = make_s_Showers(MAX_HITS);
         moduleCount++;
      }
      else
      {
         s_BarrelEMcal_t* cal = *twig;
         upshots = cal->mods->in[0].shells->in[0].cones->in[0].upstream->showers;
         downshots = cal->mods->in[0].shells->in[0].cones->in[0].downstream->showers;
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
      s_BarrelShowers_t* showers;
      float r = sqrt(x[0]*x[0]+x[1]*x[1]);
      float phi = atan2(x[1],x[0]);
      /*    int mark = layer<<10 + sector+track;*/
            int mark = (track << 25);
	    /*  printf("\nInside hitbcal, mark2 = %d  ",mark);*/
      void** twig = getTwig(&barrelEMcalTree, mark);
      if (*twig == 0) 
      {
         s_BarrelEMcal_t* cal = *twig = make_s_BarrelEMcal();
         cal->barrelShowers = showers = make_s_BarrelShowers(1);
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
         showers = ((s_BarrelEMcal_t*) *twig)->barrelShowers;
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

   if ((moduleCount == 0) && (showerCount == 0))
   {
      return 0;
   }

   /*    printf("moduleCount ,showerCount = %4d%4d !\n",moduleCount,showerCount);*/

   box = make_s_BarrelEMcal();
   box->mods = make_s_Mods(48);
   box->barrelShowers = make_s_BarrelShowers(showerCount);
   while (item = (s_BarrelEMcal_t*) pickTwig(&barrelEMcalTree))
   {
      if (item->mods)
      {
	int n = box->mods->mult;
        int module = item->mods->in[0].module;

         if (n==0)
         {
            box->mods->in[0] = item->mods->in[0];
            box->mods->in[0].shells = make_s_Shells(9);
            box->mods->mult++;
         }

         else if (module>box->mods->in[n-1].module)
	   { 
           if (box->mods->in[n-1].shells->mult == 0)
            {
               FREE(box->mods->in[--n].shells);
               box->mods->mult--;
            }
            box->mods->in[n] = item->mods->in[0];
            box->mods->in[n].shells = make_s_Shells(9);
            box->mods->mult++;
           }
         else
	 {
	   n--;
	 }   


        float E1 = item->mods->in[0].shells->in[0].cones->in[0].
        upstream->showers->in[0].E;
        float E2 = item->mods->in[0].shells->in[0].cones->in[0].
        downstream->showers->in[0].E;
         float E = E1+E2;

         int layer =  item->mods->in[0].shells->in[0].layer;
         int m = box->mods->in[n].shells->mult;
         if (m == 0)
         {
            box->mods->in[n].shells->in[0] = item->mods->in[0].shells->in[0];
            box->mods->in[n].shells->in[0].cones = make_s_Cones(4);
            box->mods->in[n].shells->mult++;
         }

         else if (layer > box->mods->in[n].shells->in[m-1].layer)

         {
	 

            if (box->mods->in[n].shells->in[m-1].cones->mult == 0)
            {
               FREE(box->mods->in[n].shells->in[--m].cones);
               box->mods->in[n].shells->mult--;
            }
            box->mods->in[n].shells->in[m] = item->mods->in[0].shells->in[0];
            box->mods->in[n].shells->in[m].cones = make_s_Cones(4);
            box->mods->in[n].shells->mult++;
         }
         else
         {
            m--;
         }
         if (E > THRESH_MEV/1000)
         {
            int mm = box->mods->in[n].shells->in[m].cones->mult++;
            box->mods->in[n].shells->in[m].cones->in[mm] = 
                                      item->mods->in[0].shells->in[0].
            cones->in[0];
         }
         else
         {
            FREE(item->mods->in[0].shells->in[0].cones->in[0].upstream->showers);
            FREE(item->mods->in[0].shells->in[0].cones->in[0].downstream->
            showers);
            FREE(item->mods->in[0].shells->in[0].cones->in[0].upstream);
            FREE(item->mods->in[0].shells->in[0].cones->in[0].downstream);
         }

         FREE(item->mods->in[0].shells->in[0].cones);
         FREE(item->mods->in[0].shells);
         FREE(item->mods); 
      }
      else if (item->barrelShowers)
      {
         int m = box->barrelShowers->mult++;
         box->barrelShowers->in[m] = item->barrelShowers->in[0];
         FREE(item->barrelShowers);
      }
      FREE(item);
   }

   /* Reduce the last <shell>...</shell> if empty */
   /*      {
      int m = box->shells->mult;
      while ((m > 0) && box->shells->in[m-1].cones->mult == 0)
      {
         FREE(box->shells->in[--m].cones);
         box->shells->mult--;
      }
   }
   */

   /*   
   {
    int n = box->mods->mult;
    
    while (n>0)
    {
    int m = box->mods->in[n-1].shells->mult;
    while (n>0 && m>0 && box->mods->in[n-1].shells->in[m-1].cones->mult == 0)   
       {
             FREE(box->mods->in[n-1].shells->in[--m].cones);
             box->mods->in[n-1].shells->mult--;	 
       }
    --n;
    } 
     
   }
   */  


   
   {
    int n = box->mods->mult;
    
    while (n>0)
    {
    int m = box->mods->in[n-1].shells->mult;
    while (n>0 && m>0 && box->mods->in[n-1].shells->in[m-1].cones->mult == 0)   
       {
             FREE(box->mods->in[n-1].shells->in[m-1].cones);
             box->mods->in[n-1].shells->mult--;	 
	     --m;
       }
    --n;
    } 
     
   }
   



   /*   
     {
          int n = box->mods->mult;
          while (n>0 && box->mods->in[n-1].shells->mult==0)
	    {
	      FREE(box->mods->in[--n].shells);
	      box->mods->mult--;
            }    
     }
   
   */


   
     {
          int n = box->mods->mult;
          while (n>0 && box->mods->in[n-1].shells->mult==0)
	    {
	      FREE(box->mods->in[n-1].shells);
	      box->mods->mult--;
	      --n;
            }    

     }
   
   
   

   moduleCount = showerCount = 0;
   return box;
}

