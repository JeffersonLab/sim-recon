/*
 * hitLGD - registers hits for forward calorimeter
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

#define TWO_HIT_RESOL   75.
#define MAX_HITS        10

hitTree_t* forwardEMcalTree = 0;

void hitForwardEMcal (float xin[4], float xout[4],
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

   /* post the hit to the hits tree, mark block as hit */
   {
      int nshot;
      s_Showers_t* shots;
      int row = getrow_();
      int column = getcolumn_();
      int mark = (row << 16) + column;
      void** twig = getTwig(&forwardEMcalTree, mark);
      if (twig == 0)
      {
         s_ForwardEMcal_t* cal = *twig = make_s_ForwardEMcal();
         cal->rows = make_s_Rows(1);
         cal->rows->mult = 1;
         cal->rows->in[0].row = row;
         cal->rows->in[0].columns = make_s_Columns(1);
         cal->rows->in[0].columns->mult = 1;
         cal->rows->in[0].columns->in[0].col = column;
         cal->rows->in[0].columns->in[0].showers =
         shots = make_s_Showers(MAX_HITS);
      }
      else
      {
         s_ForwardEMcal_t* cal = *twig;
         shots = cal->rows->in[0].columns->in[0].showers;
      }

      for (nshot = 0; nshot < shots->mult; nshot++)
      {
         if (fabs(shots->in[nshot].t - t) < TWO_HIT_RESOL)
         {
            break;
         }
      }
      if (nshot < shots->mult)		/* merge with former hit */
      {
         shots->in[nshot].t =
                       (shots->in[nshot].t * shots->in[nshot].E + t*dEsum)
                     / (shots->in[nshot].E += dEsum);
      }
      else if (nshot < MAX_HITS)         /* create new hit */
      {
         shots->in[nshot].t =
                       (shots->in[nshot].t * shots->in[nshot].E + t*dEsum)
                     / (shots->in[nshot].E += dEsum);
         shots->mult++;
      }
      else
      {
         fprintf(stderr,"HDGeant error in hitforwardEMcal: ");
         fprintf(stderr,"maximum hit count %d exceeded, quitting!\n",MAX_HITS);
         exit(2);
      }
   }

   /* post the hit to the truth tree, once per primary track */
   {
      s_ForwardShowers_t* showers;
      int mark = (track << 24);
      void** twig = getTwig(&forwardEMcalTree, mark);
      if (twig == 0)
      {
         s_ForwardEMcal_t* cal = *twig = make_s_ForwardEMcal();
         cal->forwardShowers = showers = make_s_ForwardShowers(1);
         showers->in[0].t = t;
         showers->in[0].x = x[0];
         showers->in[0].y = x[1];
         showers->in[0].E = dEsum;
         showers->mult = 1;
      }
      else
      {
         showers->in[0].x = (showers->in[0].x * showers->in[0].E + x[0]*dEsum)
                          / (showers->in[0].E + dEsum);
         showers->in[0].y = (showers->in[0].y * showers->in[0].E + x[1]*dEsum)
                          / (showers->in[0].E + dEsum);
         showers->in[0].t = (showers->in[0].t * showers->in[0].E + t*dEsum)
                          / (showers->in[0].E += dEsum);
      }
   }
}

/* entry point from fortran */

void hitforwardemcal_(float* xin, float* xout,
                      float* pin, float* pout, float* dEsum, int* track)
{
   hitForwardEMcal(xin,xout,pin,pout,*dEsum,*track);
}
