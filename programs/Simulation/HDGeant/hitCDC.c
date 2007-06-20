/*
 * hitCDC - registers hits for Central Drift Chamber
 *
 *	This is a part of the hits package for the
 *	HDGeant simulation program for Hall D.
 *
 *	version 1.0 	-Richard Jones July 16, 2001
 *
 * changes: Wed Jun 20 13:19:56 EDT 2007 B. Zihlmann 
 *          add ipart to the function call hitCentralDC
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <hddm_s.h>
#include <geant3.h>
#include <bintree.h>

// Drift speed 2.2cm/us is appropriate for a 90/10 Argon/Methane mixture
#define DRIFT_SPEED	.0055
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
                   int track, int stack, int history, int ipart )
{
   float x[3], t;
   float dx[3], dr;
   float dEdx;
   float xlocal[3];
   float xinlocal[3];
   float xoutlocal[3];
   float dradius,drin,drout;

   x[0] = (xin[0] + xout[0])/2;
   x[1] = (xin[1] + xout[1])/2;
   x[2] = (xin[2] + xout[2])/2;
   t    = (xin[3] + xout[3])/2 * 1e9;
   dx[0] = xin[0] - xout[0];
   dx[1] = xin[1] - xout[1];
   dx[2] = xin[2] - xout[2];
   transformCoord(xin,"global",xinlocal,"local");
   transformCoord(xout,"global",xoutlocal,"local");
   xlocal[0] = (xinlocal[0] + xoutlocal[0])/2;
   xlocal[1] = (xinlocal[1] + xoutlocal[1])/2;
   xlocal[2] = (xinlocal[2] + xoutlocal[2])/2;

   dradius = sqrt(xlocal[0]*xlocal[0] + xlocal[1]*xlocal[1]);
   drin = sqrt(xinlocal[0]*xinlocal[0] + xinlocal[1]*xinlocal[1]);
   drout = sqrt(xoutlocal[0]*xoutlocal[0] + xoutlocal[1]*xoutlocal[1]);

/* If the particle hits the wire or exits from the end of the tube,
 * then the midpoint from entrance to exit will not actually correspond
 * to the DOCA.  This is important since many tracks that exit through
 * the CDC endplate will exit through the end of the straw tube. In either
 * of these two cases, there are two possibilities:
 *  1.) The DOCA is at the exit point itself
 *  2.) The DOCA is at the mid-point between the entrance and where the
 *      exit would have been if the tube were infinitely long.
 *
 * First, we determine if the particle exited the end of the tube.
 * If so, then calculate the DOCA (dradius) for both possibilities
 * and keep the one that is smallest.
 */
   if (drin < 0.01)		/* particle is leaving the wire */
   {
      dradius = drin;
      t = xin[3]*1e9;
   }
   else if (drout < dradius)    /* particle is entering the wire */
   {
      dradius = drout;
      t = xout[3]*1e9;
   }
   else if (fabs(drout-drin) > 0.01)
   {
   /* Particle exited from end of the straw.
    *
    * Calculate exit point from an infinite straw 
    * We do this by defining the direction of the
    * track and finding the amount we need to extend 
    * in that direction in order to be at the tube
    * radius (determined by the entrance point).
    *
    * xout_local = xin_local + alpha*trackdir
    *
    * where xout_local, xin_local, and trackdir are all
    * vectors. "alpha" is a scaler multiplier to be
    * be solved for. The values of xin_local and trackdir
    * are detemined from xin and xout, while xout_local
    * is to be calculated once alpha is determined.
    *
    * We solve for alpha by setting the transverse
    * distance of xout_local to drin, the radius
    * of the tube. This leads to an equation quadratic
    * in alpha.  */

      float alpha;
      float A,B,C;
      float trackdir[3];
      float xoutlocal_i[3], xout_i[3];
      float docaout;
      transformCoord(dx,"global",trackdir,"local");
      A = trackdir[0]*trackdir[0] + trackdir[1]*trackdir[1];
      B = 2.0*(trackdir[0]*xinlocal[0] + trackdir[1]*xinlocal[1]);
      C = drin*drin + xinlocal[0]*xinlocal[0] + xinlocal[1]*xinlocal[1];
      
      alpha = (-B + sqrt(B*B - 4.0*A*C))/(2.0*A);
      xoutlocal_i[0] = xinlocal[0] + alpha*trackdir[0];
      xoutlocal_i[1] = xinlocal[1] + alpha*trackdir[1];
      xoutlocal_i[2] = xinlocal[2] + alpha*trackdir[2];
      docaout = sqrt(xoutlocal_i[0]*xoutlocal_i[0]
                     + xoutlocal_i[1]*xoutlocal_i[1]);
      transformCoord(xoutlocal_i,"local",xout_i,"global");
      
    /* Check which is the smallest DOCA and copy that into
     * dradius so it will be used for the drift time below.  */
      if (drout < dradius)
      {
         dradius = drout;
         x[0]=xoutlocal[0];
         x[1]=xoutlocal[1];
         x[2]=xoutlocal[2];
         t = xout[3]*1e9;
      }
      if (docaout < dradius)
      {
         dradius = docaout;
         x[0]=xoutlocal_i[0];
         x[1]=xoutlocal_i[1];
         x[2]=xoutlocal_i[2];
         t = xout[3]*1e9;
      }
   }

   /* Calculate dE/dx */

   dr = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
   if (dr > 1e-3)
   {
      dEdx = dEsum/dr;
   }
   else
   {
      dEdx = 0;
   }

   /* post the hit to the truth tree */

   if (history == 0)
   {
      int mark = (1<<30) + pointCount;
      void** twig = getTwig(&centralDCTree, mark);
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

   /* post the hit to the hits tree, mark sector as hit */

   if (dEsum > 0)
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
                   int* track, int* stack, int* history, int* ipart)
{
   hitCentralDC(xin,xout,pin,pout,*dEsum,*track,*stack,*history, *ipart);
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
      int straw;
#if CATHODE_STRIPS_IN_CDC
      s_CdcCathodetrips_t* strips = item->cdcCathodeStrips;
      int strip;
#endif
      s_CdcTruthPoints_t* points = item->cdcTruthPoints;
      int point;
      for (straw=0; straw < straws->mult; ++straw)
      {
         int m = box->cdcStraws->mult;

         s_CdcStrawHits_t* hits = straws->in[straw].cdcStrawHits;

         /* compress out the hits below threshold */
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
            }
         }

         if (iok)
         {
            hits->mult = iok;
            box->cdcStraws->in[m] = straws->in[straw];
            box->cdcStraws->mult++;
         }
         else if (hits != HDDM_NULL)
         {
            FREE(hits);
         }
      }
      if (straws != HDDM_NULL)
      {
         FREE(straws);
      }
#if CATHODE_STRIPS_IN_CDC
      for (strip=0; strip < strips->mult; ++strip)
      {
         int m = box->cdcCathodeStrips->mult;

         s_CdcCathodeHits* hits = strips->in[strip].cdcCathodeHits;

         /* compress out the hits below threshold */
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
            }
         }

         if (iok)
         {
            hits->mult = iok;
            box->cdcCathodeStrips->in[m] = strips->in[strip];
            box->cdcCathodeStrips->mult++;
         }
         else if (hits != HDDM_NULL)
         {
            FREE(hits);
         }
      }
      if (strips != HDDM_NULL)
      {
         FREE(strips);
      }
#endif
      for (point=0; point < points->mult; ++point)
      {
         int m = box->cdcTruthPoints->mult++;
         box->cdcTruthPoints->in[m] = points->in[point];
      }
      if (points != HDDM_NULL)
      {
         FREE(points);
      }
      FREE(item);
   }

   strawCount = stripCount = pointCount = 0;

#if CATHODE_STRIPS_IN_CDC
   if ((box->cdcCathodeStrips != HDDM_NULL) &&
       (box->cdcCathodeStrips->mult == 0))
   {
      FREE(box->cdcCathodeStrips);
      box->cdcCathodeStrips = HDDM_NULL;
   }
#endif
   if ((box->cdcStraws != HDDM_NULL) &&
       (box->cdcStraws->mult == 0))
   {
      FREE(box->cdcStraws);
      box->cdcStraws = HDDM_NULL;
   }
   if ((box->cdcTruthPoints != HDDM_NULL) &&
       (box->cdcTruthPoints->mult == 0))
   {
      FREE(box->cdcTruthPoints);
      box->cdcTruthPoints = HDDM_NULL;
   }
   if ((box->cdcStraws->mult == 0) &&
       (box->cdcTruthPoints->mult == 0))
   {
#if CATHODE_STRIPS_IN_CDC
    if (box->cdcCathodeStrips->mult == 0) {
#endif
      FREE(box);
      box = HDDM_NULL;
#if CATHODE_STRIPS_IN_CDC
    }
#endif
   }
   return box;
}
