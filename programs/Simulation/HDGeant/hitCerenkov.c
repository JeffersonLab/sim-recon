/*
 * hitCerenkov - registers hits for Cerenkov counter
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
#define TWO_HIT_RESOL   25.
#define PHI_SECTORS	24
#define N0_FIGURE	60
#define REFR_INDEX	1.0017

hitTree_t* cerenkovTree = 0;

void hitCerenkov (float xin[4], float xout[4],
                  float pin[5], float pout[5], float dEsum, int track)
{
   float x[3], t;
   float dx[3], dr;
   float dEdx;
   float xlocal[3];
   const int one = 1;
   double Eave = (pin[3] + pout[3])/2;
   double pave = (pin[4] + pout[4])/2;
   double beta = pave/Eave;
   double costheta = 1/(beta*REFR_INDEX);

   if (costheta > 1) return;		/* no light below threshold */

   gmtod_(x,xlocal,&one);
   x[0] = (xin[0] + xout[0])/2;
   x[1] = (xin[1] + xout[1])/2;
   x[2] = (xin[2] + xout[2])/2;
   t    = (xin[3] + xout[3])/2 * 1e9;
   dx[0] = xin[0] - xout[0];
   dx[1] = xin[1] - xout[1];
   dx[2] = xin[2] - xout[2];
   dr = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);


   /* post the hit to the hits tree, mark sector as hit */
   {
      int nshot;
      s_Flashes_t* shots;
      int sector = getsector_();
      double sintheta2 = 1 - costheta*costheta;
      float pe = N0_FIGURE * sintheta2 * dr;
      void** twig = getTwig(&cerenkovTree, sector);
      if (twig == 0)
      {
         s_Cerenkov_t* ckov = *twig = make_s_Cerenkov();
         ckov->sections = make_s_Sections(1);
         ckov->sections->mult = 1;
         ckov->sections->in[0].flashes = shots = make_s_Flashes(MAX_HITS);
      }
      else
      {
         s_Cerenkov_t* ckov = *twig;
         shots = ckov->sections->in[0].flashes;
      }

      for (nshot = 0; nshot < shots->mult; nshot++)
      {
         if (fabs(shots->in[nshot].t - t) < TWO_HIT_RESOL)
         {
            break;
         }
      }
      if (nshot < shots->mult)            /* merge with former hit */
      {
         shots->in[nshot].t = (shots->in[nshot].t * shots->in[nshot].pe + t*pe)
                            / (shots->in[nshot].pe += pe);
      }
      else if (nshot < MAX_HITS)         /* create new shot */
      {
         shots->in[nshot].t = t;
         shots->in[nshot].pe = pe;
         shots->mult++;
      }
      else
      {
         fprintf(stderr,"HDGeant error in hitCerenkov: ");
         fprintf(stderr,"maximum shot count %d exceeded, quitting!\n",MAX_HITS);
         exit(2);
      }
   }

   /* no truth sections here, for the moment */

}

/* entry point from fortran */

void hitcerenkov_(float* xin, float* xout,
                  float* pin, float* pout, float* dEsum, int* track)
{
   hitCerenkov(xin,xout,pin,pout,*dEsum,*track);
}
