/*
 * hitCDC - registers hits for Central Drift Chamber
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

#define Z0		-100.
#define TWO_HIT_RESOL	25.
#define MAX_HITS 	10

hitTree_t* centralDCTree = 0;

void hitCentralDC (float xin[4], float xout[4],
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
      int nhit;
      s_Hits_t* hits;
      int layer = getlayer_();
      int cell = getcell_();
      int ring = getring_();
      int count = getcount_();
      int sector = getsector_();
      if (layer == 0)		/* in a straw */
      {
         int mark = (ring << 20) + (sector << 8);
         void** twig = getTwig(&centralDCTree, mark);
         if (twig == 0)
         {
            float phim = (sector - 0.5) * (2*M_PI/count);
            s_CentralDC_t* cdc = *twig = make_s_CentralDC();
            cdc->rings = make_s_Rings(1);
            cdc->rings->mult = 1;
            cdc->rings->in[0].radius =
                        sqrt(xlocal[0]*xlocal[0] + xlocal[1]*xlocal[1]);
            cdc->rings->in[0].straws = make_s_Straws(1);
            cdc->rings->in[0].straws->mult = 1;
            cdc->rings->in[0].straws->in[0].phim = phim;
            cdc->rings->in[0].straws->in[0].hits =
            hits = make_s_Hits(MAX_HITS);
         }
         else
         {
            s_CentralDC_t* cdc = (s_CentralDC_t*) *twig;
            hits = cdc->rings->in[0].straws->in[0].hits;
         }
      }
      else			/* in one of the z-readout (band) layers */
      {
         float phi = atan2(xlocal[1],xlocal[0]);
         int mark = (layer << 28) + (sector << 20) + (cell << 8);
         void** twig = getTwig(&centralDCTree, mark);
         if (twig == 0)
         {
            s_CentralDC_t* cdc = *twig = make_s_CentralDC();
            cdc->cathodeCyls = make_s_CathodeCyls(1);
            cdc->cathodeCyls->mult = 1;
            cdc->cathodeCyls->in[0].radius = 
                              sqrt(xlocal[0]*xlocal[0] + xlocal[1]*xlocal[1]);
            cdc->cathodeCyls->in[0].bands = make_s_Bands(1);
            cdc->cathodeCyls->in[0].bands->mult = 1;
            cdc->cathodeCyls->in[0].bands->in[0].sector = sector;
            cdc->cathodeCyls->in[0].bands->in[0].z = xlocal[2];
            cdc->cathodeCyls->in[0].bands->in[0].hits =
            hits = make_s_Hits(MAX_HITS);
         }
         else
         {
            s_CentralDC_t* cdc = (s_CentralDC_t*) *twig;
            hits = cdc->cathodeCyls->in[0].bands->in[0].hits;
         }
      }

      for (nhit = 0; nhit < hits->mult; nhit++)
      {
         if (fabs(hits->in[nhit].t - t) < TWO_HIT_RESOL)
         {
            break;
         }
      }
      if (nhit < hits->mult)		/* merge with former hit */
      {
         hits->in[nhit].t = (hits->in[nhit].t * hits->in[nhit].dE + t*dEsum)
                          / (hits->in[nhit].dE += dEsum);
      }
      else if (nhit < MAX_HITS)		/* create new hit */
      {
         hits->in[nhit].t = t;
         hits->in[nhit].dE = dEsum;
         hits->mult++;
      }
      else
      {
         fprintf(stderr,"HDGeant error in hitCentralDC: ");
         fprintf(stderr,"max hit count %d exceeded, quitting!\n",MAX_HITS);
         exit(2);
      }

      /* post the hit to the truth tree, once per primary track per straw */

      if (layer == 0)
      {
         float phim = (sector - 0.5) * (2*M_PI/count);
         int mark = (ring << 20) + (sector << 8) + track;
         void** twig = getTwig(&centralDCTree, mark);
         if (twig == 0)
         {
            s_CentralDC_t* cdc = *twig = make_s_CentralDC();
            s_CdcPoints_t* points = make_s_CdcPoints(1);
            cdc->rings = make_s_Rings(1);
            cdc->rings->mult = 1;
            cdc->rings->in[0].radius = sqrt(x[0]*x[0] + x[1]*x[1]);
            cdc->rings->in[0].straws = make_s_Straws(1);
            cdc->rings->in[0].straws->mult = 1;
            cdc->rings->in[0].straws->in[0].phim = phim;
            cdc->rings->in[0].straws->in[0].cdcPoints = points;
            points->in[0].z = x[2];
            points->in[0].phi = atan2(x[1],x[0]);
            points->in[0].dradius = sqrt(xlocal[0]*xlocal[0] +
                                         xlocal[1]*xlocal[1]);
            points->in[0].dEdx = dEdx;
            points->mult = 1;
         }
      }
   }
}

/* entry points from fortran */

void hitcentraldc_(float* xin, float* xout,
                   float* pin, float* pout, float* dEsum, int* track)
{
   hitCentralDC(xin,xout,pin,pout,*dEsum,*track);
}
