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
#include <math.h>

#include <hddm_s.h>
#include <geant3.h>
#include <bintree.h>

#define DRIFT_SPEED	.0022
#define TWO_HIT_RESOL	250.
#define MAX_HITS 	100

binTree_t* centralDCTree = 0;
static int strawCount = 0;
static int bandCount = 0;
static int pointCount = 0;


/* register hits during tracking (from gustep) */

void hitCentralDC (float xin[4], float xout[4],
                   float pin[5], float pout[5], float dEsum, int track)
{
   float x[3], t;
   float dx[3], dr;
   float dEdx;
   float xlocal[3];
   float xcdc[3];
   float xHat[] = {1,0,0};

   if (dEsum == 0) return;              /* only seen if it deposits energy */

   x[0] = (xin[0] + xout[0])/2;
   x[1] = (xin[1] + xout[1])/2;
   x[2] = (xin[2] + xout[2])/2;
   t    = (xin[3] + xout[3])/2 * 1e9;
   transformCoord(x,"global",xlocal,"local");
   transformCoord(xHat,"local",xcdc,"CDC ");

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
      int sector = getsector_();
      float phim = atan2(xcdc[1],xcdc[0]);
      float dradius = sqrt(xlocal[0]*xlocal[0] + xlocal[1]*xlocal[1]);
      float tdrift = t + dradius/DRIFT_SPEED;
      if (layer == 0)		/* in a straw */
      {
         int mark = (ring << 20) + (sector << 8);
         void** twig = getTwig(&centralDCTree, mark);
         if (*twig == 0)
         {
            s_CentralDC_t* cdc = *twig = make_s_CentralDC();
            cdc->rings = make_s_Rings(1);
            cdc->rings->mult = 1;
            cdc->rings->in[0].radius =
                              sqrt(xcdc[0]*xcdc[0] + xcdc[1]*xcdc[1]) - 1;
            cdc->rings->in[0].straws = make_s_Straws(1);
            cdc->rings->in[0].straws->mult = 1;
            cdc->rings->in[0].straws->in[0].phim = phim;
            cdc->rings->in[0].straws->in[0].hits =
            hits = make_s_Hits(MAX_HITS);
            strawCount++;
         }
         else
         {
            s_CentralDC_t* cdc = (s_CentralDC_t*) *twig;
            hits = cdc->rings->in[0].straws->in[0].hits;
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
            hits->in[nhit].t =
                    (hits->in[nhit].t * hits->in[nhit].dE + tdrift * dEsum)
                  / (hits->in[nhit].dE += dEsum);
         }
         else if (nhit < MAX_HITS)		/* create new hit */
         {
            hits->in[nhit].t = tdrift;
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

         mark = (ring << 20) + (sector << 8) + track;
         twig = getTwig(&centralDCTree, mark);
         if (*twig == 0)
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
            points->in[0].track = track;
            points->in[0].z = x[2];
            points->in[0].phi = atan2(x[1],x[0]);
            points->in[0].dradius = dradius;
            points->in[0].dEdx = dEdx;
            points->mult = 1;
            pointCount++;
         }
      }

      else			/* in one of the z-readout (band) layers */
      {
         int mark = (layer << 28) + (sector << 20) + (cell << 8);
         void** twig = getTwig(&centralDCTree, mark);
         if (*twig == 0)
         {
            s_CentralDC_t* cdc = *twig = make_s_CentralDC();
            cdc->cathodeCyls = make_s_CathodeCyls(1);
            cdc->cathodeCyls->mult = 1;
            cdc->cathodeCyls->in[0].radius = 
                              sqrt(xlocal[0]*xlocal[0] + xlocal[1]*xlocal[1]);
            cdc->cathodeCyls->in[0].bands = make_s_Bands(1);
            cdc->cathodeCyls->in[0].bands->mult = 1;
            cdc->cathodeCyls->in[0].bands->in[0].phim = phim;
            cdc->cathodeCyls->in[0].bands->in[0].z = x[2];
            bandCount++;
         }
         else
         {
            s_CentralDC_t* cdc = (s_CentralDC_t*) *twig;
            cdc->cathodeCyls->in[0].bands->in[0].z += x[2];
            cdc->cathodeCyls->in[0].bands->in[0].z /= 2;
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


/* pick and package the hits for shipping */

s_CentralDC_t* pickCentralDC ()
{
   s_CentralDC_t* box;
   s_CentralDC_t* item;

   if ((strawCount == 0) && (bandCount == 0) && (pointCount == 0))
   {
      return 0;
   }

   box = make_s_CentralDC();
   box->cathodeCyls = make_s_CathodeCyls(5);
   box->rings = make_s_Rings(32);
   while (item = (s_CentralDC_t*) pickTwig(&centralDCTree))
   {
      if (item->cathodeCyls)
      {
         float r = item->cathodeCyls->in[0].radius;
         int m = box->cathodeCyls->mult;
         if ((m == 0) || (r > box->cathodeCyls->in[m-1].radius + 0.5))
         {
            box->cathodeCyls->in[m] = item->cathodeCyls->in[0];
            box->cathodeCyls->in[m].bands = make_s_Bands(bandCount);
            box->cathodeCyls->mult++;
         }
         else
         {
            m--;
         }
         {
            int mm = box->cathodeCyls->in[m].bands->mult++;
            box->cathodeCyls->in[m].bands->in[mm] =
                         item->cathodeCyls->in[0].bands->in[0];
         }
         FREE(item->cathodeCyls->in[0].bands);
         FREE(item->cathodeCyls);
      }
      else if (item->rings)
      {
         float r = item->rings->in[0].radius;
         int m = box->rings->mult;
         if ((m == 0) || (r > box->rings->in[m-1].radius + 0.5))
         {
            box->rings->in[m] = item->rings->in[0];
            box->rings->in[m].straws = make_s_Straws(strawCount);
            box->rings->mult++;
         }
         else
         {
            m--;
         }
         {
            float phim = item->rings->in[0].straws->in[0].phim;
            int mm = box->rings->in[m].straws->mult;
            if (item->rings->in[0].straws->in[0].hits)
            {
               box->rings->in[m].straws->in[mm] =
                                          item->rings->in[0].straws->in[0];
               box->rings->in[m].straws->mult++;
            }
            else if ((mm == 0) ||
                     (phim > box->rings->in[m].straws->in[mm-1].phim))
            {
               box->rings->in[m].straws->in[mm] =
                                          item->rings->in[0].straws->in[0];
               box->rings->in[m].straws->in[mm].cdcPoints =
                                          make_s_CdcPoints(pointCount);
               box->rings->in[m].straws->mult++;
            }
            else
            {
               mm--;
            }
            if (item->rings->in[0].straws->in[0].cdcPoints)
            {
               int mmm;
               if (box->rings->in[m].straws->in[mm].cdcPoints == 0)
               {
                  box->rings->in[m].straws->in[mm].cdcPoints =
                                            make_s_CdcPoints(pointCount);
               }
               mmm = box->rings->in[m].straws->in[mm].cdcPoints->mult++;
               box->rings->in[m].straws->in[mm].cdcPoints->in[mmm] =
                        item->rings->in[0].straws->in[0].cdcPoints->in[0];
               FREE(item->rings->in[0].straws->in[0].cdcPoints);
            }
         }
         FREE(item->rings->in[0].straws);
         FREE(item->rings);
      }
      FREE(item);
   }
   strawCount = bandCount = pointCount = 0;
   return box;
}
