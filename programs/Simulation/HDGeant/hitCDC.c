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
#define THRESH_KEV	1.

binTree_t* centralDCTree = 0;
static int strawCount = 0;
static int pointCount = 0;
static int stripCount = 0;


/* register hits during tracking (from gustep) */

void hitCentralDC (float xin[4], float xout[4],
                   float pin[5], float pout[5], float dEsum,
                   int track, int stack)
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
      s_CdcStrawHits_t* hits;
#if CATHODE_STRIPS_IN_CDC
      s_CdcCathodeStrips_t* chits;
      int cell = getcell_();
#endif
      int layer = getlayer_();
      int ring = getring_();
      int sector = getsector_();
      float dradius = sqrt(xlocal[0]*xlocal[0] + xlocal[1]*xlocal[1]);
      float tdrift = t + dradius/DRIFT_SPEED;
      if (layer == 0)		/* in a straw */
      {
         int mark = (ring<<20) + sector;
         void** twig = getTwig(&centralDCTree, mark);
         if (*twig == 0)
         {
            s_CentralDC_t* cdc = *twig = make_s_CentralDC();
            s_CdcStraws_t* straws = make_s_CdcStraws(1);
            straws->mult = 1;
            straws->in[0].ring = ring;
            straws->in[0].straw = sector;
            straws->in[0].cdcStrawHits = hits = make_s_CdcStrawHits(MAX_HITS);
            cdc->cdcStraws = straws;
            strawCount++;
         }
         else
         {
            s_CentralDC_t* cdc = (s_CentralDC_t*) *twig;
            hits = cdc->cdcStraws->in[0].cdcStrawHits;
         }

         for (nhit = 0; nhit < hits->mult; nhit++)
         {
            if (fabs(hits->in[nhit].t - tdrift) < TWO_HIT_RESOL)
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
            fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
         }

      /* post the hit to the truth tree, once per primary track per ring */

         mark = (1<<30) + (ring<<20) + track;
         twig = getTwig(&centralDCTree, mark);
         if (*twig == 0)
         {
            s_CentralDC_t* cdc = *twig = make_s_CentralDC();
            s_CdcTruthPoints_t* points = make_s_CdcTruthPoints(1);
            points->in[0].primary = (stack == 0);
            points->in[0].track = track;
            points->in[0].z = x[2];
            points->in[0].r = sqrt(x[0]*x[0] + x[1]*x[1]);
            points->in[0].phi = atan2(x[1],x[0]);
            points->in[0].dradius = dradius;
            points->in[0].dEdx = dEdx;
            points->mult = 1;
            cdc->cdcTruthPoints = points;
            pointCount++;
         }
      }

#if CATHODE_STRIPS_IN_CDC
      else			/* in one of the z-readout (strip) layers */
      {
         int nchit;
         int mark = (layer << 28) + (sector << 20) + (cell << 8);
         void** twig = getTwig(&centralDCTree, mark);
         if (*twig == 0)
         {
            s_CentralDC_t* cdc = *twig = make_s_CentralDC();
            s_CdcCathodeStrips strips = make_s_CdcCathodeStrips(1);
            strips->mult = 1;
            strips->in[0].layer = layer;
            strips->in[0].sector = sector;
            strips->in[0].strip = cell;
            strips->in[0].cdcCathodeHits = chits
                                       = make_s_CdcCathodeHits(MAX_HITS);
            cdc->cdcCathodeStrips = strips;
            stripCount++;
         }
         else
         {
            s_CentralDC_t* cdc = (s_CentralDC_t*) *twig;
            chits = cdc->cdcCathodeStrips->in[0].cdcCathodeHits;
         }

         for (nchit = 0; nchit < chits->mult; nchit++)
         {
            if (fabs(chits->in[nchit].t - tdrift) < TWO_HIT_RESOL)
            {
               break;
            }
         }
         if (nchit < chits->mult)		/* merge with former hit */
         {
            chits->in[nchit].t =
                    (chits->in[nchit].t * chits->in[nchit].dE + t * dEsum)
                  / (chits->in[nchit].dE += dEsum);
         }
         else if (nchit < MAX_HITS)		/* create new hit */
         {
            chits->in[nchit].t = t;
            chits->in[nchit].dE = dEsum;
            chits->mult++;
         }
         else
         {
            fprintf(stderr,"HDGeant error in hitCentralDC: ");
            fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
         }
      }
#endif
   }
}

/* entry points from fortran */

void hitcentraldc_(float* xin, float* xout,
                   float* pin, float* pout, float* dEsum,
                   int* track, int* stack)
{
   hitCentralDC(xin,xout,pin,pout,*dEsum,*track,*stack);
}


/* pick and package the hits for shipping */

s_CentralDC_t* pickCentralDC ()
{
   s_CentralDC_t* box;
   s_CentralDC_t* item;

   if ((strawCount == 0) && (stripCount == 0) && (pointCount == 0))
   {
      return HDDM_NULL;
   }

   box = make_s_CentralDC();
   box->cdcStraws = make_s_CdcStraws(strawCount);
   box->cdcTruthPoints = make_s_CdcTruthPoints(pointCount);
#if CATHODE_STRIPS_IN_CDC
   box->cdcCathodeStrips = make_s_CdcCathodeStrips(stripCount);
#endif
   while (item = (s_CentralDC_t*) pickTwig(&centralDCTree))
   {
      s_CdcStraws_t* straws = item->cdcStraws;
#if CATHODE_STRIPS_IN_CDC
      s_CdcCathodetrips_t* strips = item->cdcCathodeStrips;
#endif
      s_CdcTruthPoints_t* points = item->cdcTruthPoints;
      if (straws != HDDM_NULL)
      {
         int m = box->cdcStraws->mult;
         int mok=0;

         /* compress out the hits below threshold */
         s_CdcStrawHits_t* hits = straws->in[0].cdcStrawHits;
         if (hits != HDDM_NULL)
         {
            int i,iok;
            for (iok=i=0; i < hits->mult; i++)
            {
               if (hits->in[i].dE >= THRESH_KEV/1e6)
               {
                  if (iok < i)
                  {
                     hits->in[iok] = hits->in[i];
                  }
                  ++iok;
                  ++mok;
               }
            }
            hits->mult = iok;

            if (mok)
            {
               box->cdcStraws->in[m] = straws->in[0];
               box->cdcStraws->mult++;
            }
            else
            {
               FREE(hits);
            }
         }
         FREE(straws);
      }
#if CATHODE_STRIPS_IN_CDC
      else if (strips != HDDM_NULL)
      {
         int m = box->cdcCathodeStrips->mult;
         int mok=0;

         /* compress out the hits below threshold */
         s_CdcCathodeHits* hits = strips->in[0].cdcCathodeHits;
         if (hits != HDDM_NULL)
         {
            int i,iok;
            for (iok=i=0; i < hits->mult; i++)
            {
               if (hits->in[i].dE >= THRESH_KEV/1e6)
               {
                  if (iok < i)
                  {
                     hits->in[iok] = hits->in[i];
                  }
                  ++iok;
                  ++mok;
               }
            }
            hits->mult = mok;

            if (mok)
            {
               box->cdcCathodeStrips->in[m] = strips->in[0];
               box->cdcCathodeStrips->mult++;
            }
            else
            {
               FREE(hits);
            }
         }
         FREE(item->cdcCathodeStrips);
      }
#endif
      else if (points != HDDM_NULL)
      {
         int m = box->cdcTruthPoints->mult++;
         box->cdcTruthPoints->in[m] = points->in[0];
         FREE(points);
      }
      FREE(item);
   }
   strawCount = stripCount = pointCount = 0;

   return box;
}
