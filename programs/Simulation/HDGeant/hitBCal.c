/*
 * hitBCal - registers hits for barrel calorimeter
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

#define Z0		15.
#define TWO_HIT_RESOL	25.
#define ATTEN_LENGTH	150.
#define C_EFFECTIVE	15.
#define MAX_HITS 	99

hitTree_t* barrelEMcalTree = 0;

void hitBarrelEMcal (float xin[4], float xout[4],
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

   /* post the hit to the hits tree, mark sector as hit */
   {
      int nshot;
      s_Showers_t* upshots;
      s_Showers_t* downshots;
      int sector = getsector_();
      float tup = t + xlocal[2] * C_EFFECTIVE;
      float tdown = t - xlocal[2] * C_EFFECTIVE;
      float Eup = dEsum * exp(-xlocal[2] / ATTEN_LENGTH);
      float Edown = dEsum * exp(+xlocal[2] / ATTEN_LENGTH);
      void** twig = getTwig(&barrelEMcalTree, sector);
      if (twig == 0)
      {
         s_BarrelEMcal_t* cal = *twig = make_s_BarrelEMcal();
         cal->modules = make_s_Modules(1);
         cal->modules->mult = 1;
         cal->modules->in[0].upstream = make_s_Upstream();
         cal->modules->in[0].downstream = make_s_Downstream();
         cal->modules->in[0].upstream->showers =
         upshots = make_s_Showers(MAX_HITS);
         cal->modules->in[0].downstream->showers =
         downshots = make_s_Showers(MAX_HITS);
      }
      else
      {
         s_BarrelEMcal_t* cal = *twig;
         upshots = cal->modules->in[0].upstream->showers;
         downshots = cal->modules->in[0].downstream->showers;
      }

      for (nshot = 0; nshot < upshots->mult; nshot++)
      {
         if (fabs(upshots->in[nshot].t - t) < TWO_HIT_RESOL)
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
         upshots->in[nshot].t =
                    (upshots->in[nshot].t * upshots->in[nshot].E + tup * Eup)
                  / (upshots->in[nshot].E += Eup);
         upshots->mult++;
         downshots->in[nshot].t = 
            (downshots->in[nshot].t * downshots->in[nshot].E +  tdown * Edown)
          / (downshots->in[nshot].E += Edown);
         downshots->mult++;
      }
      else
      {
         fprintf(stderr,"HDGeant error in hitBarrelEMcal: ");
         fprintf(stderr,"maximum hit count %d exceeded, quitting!\n",MAX_HITS);
         exit(2);
      }
   }

   /* post the hit to the truth tree, once per primary track */
   {
      s_BarrelShowers_t* showers;
      float phi = atan2(x[1],x[0]);
      int mark = (track << 16);
      void** twig = getTwig(&barrelEMcalTree, mark);
      if (twig == 0)
      {
         s_BarrelEMcal_t* cal = *twig = make_s_BarrelEMcal();
         cal->barrelShowers = showers = make_s_BarrelShowers(1);
         showers->in[0].z = x[2];
         showers->in[0].phi = phi;
         showers->in[0].t = t;
         showers->in[0].E = dEsum;
         showers->mult = 1;
      }
      else
      {
         showers->in[0].z = (showers->in[0].z * showers->in[0].E + x[2]*dEsum)
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
                     float* pin, float* pout, float* dEsum, int* track)
{
   hitBarrelEMcal(xin,xout,pin,pout,*dEsum,*track);
}
