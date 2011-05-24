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

#include "calibDB.h"
extern s_HDDM_t* thisInputEvent;

// Drift speed 2.2cm/us is appropriate for a 90/10 Argon/Methane mixture
static float DRIFT_SPEED     =   0.0055;
static float TWO_HIT_RESOL   =  25.;
static int   MAX_HITS 	     = 100;
static float THRESH_KEV	     =   1.;
static float STRAW_RADIUS    =   0.8;

binTree_t* centralDCTree = 0;
static int strawCount = 0;
static int pointCount = 0;
static int stripCount = 0;
static int initialized = 0;

/* void GetDOCA(int ipart, float x[3], float p[5], float doca[3]);  disabled 6/24/2009 */

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
  float doca[3]; 
  float trackdir[3];
  float alpha;
  
  if (!initialized) {
    mystr_t strings[50];
    float values[50];
    int nvalues = 50;
    int status = GetConstants("CDC/cdc_parms", &nvalues, values, strings);
    
    if (!status) {
      int ncounter = 0;
      int i;
      for ( i=0;i<(int)nvalues;i++){
	//printf("%d %s \n",i,strings[i].str);
	if (!strcmp(strings[i].str,"CDC_DRIFT_SPEED")) {
	  DRIFT_SPEED  = values[i];
	  ncounter++;
	}
	if (!strcmp(strings[i].str,"CDC_TWO_HIT_RESOL")) {
	  TWO_HIT_RESOL  = values[i];
	  ncounter++;
	}
	if (!strcmp(strings[i].str,"CDC_MAX_HITS")) {
	  MAX_HITS  = (int)values[i];
	  ncounter++;
	}
	if (!strcmp(strings[i].str,"CDC_THRESH_KEV")) {
	  THRESH_KEV  = values[i];
	  ncounter++;
	}
	if (!strcmp(strings[i].str,"CDC_STRAW_RADIUS")) {
	  THRESH_KEV  = values[i];
	  ncounter++;
	}
      }
      if (ncounter==5){
	printf("CDC: ALL parameters loaded from Data Base\n");
      } else if (ncounter<5){
	printf("CDC: NOT ALL necessary parameters found in Data Base %d out of 5\n",ncounter);
      } else {
	printf("CDC: SOME parameters found more than once in Data Base\n");
      } 	
    }
    initialized = 1;
  }
    
  x[0] = (xin[0] + xout[0])/2;
  x[1] = (xin[1] + xout[1])/2;
  x[2] = (xin[2] + xout[2])/2;
  t    = (xin[3] + xout[3])/2 * 1e9;
  dx[0] = xin[0] - xout[0];
  dx[1] = xin[1] - xout[1];
  dx[2] = xin[2] - xout[2];
  transformCoord(xin,"global",xinlocal,"local");
  transformCoord(xout,"global",xoutlocal,"local");
  /*
  xlocal[0] = (xinlocal[0] + xoutlocal[0])/2;
  xlocal[1] = (xinlocal[1] + xoutlocal[1])/2;
  xlocal[2] = (xinlocal[2] + xoutlocal[2])/2;
  */
  
  /* For particles that range out inside the active volume, the
   * "out" time seems to be set to something enormously high.
   * This screws up the hit. Check for this case here by looking
   * at xout[3] and making sure it is less than 1 second. If it's
   * not, then just use xin[3] for "t".
   */
  if(xout[3] > 1.0) t = xin[3] * 1e9;
	 
  drin = sqrt(xinlocal[0]*xinlocal[0] + xinlocal[1]*xinlocal[1]);
  drout = sqrt(xoutlocal[0]*xoutlocal[0] + xoutlocal[1]*xoutlocal[1]);
  
  trackdir[0] =-xinlocal[0] + xoutlocal[0];
  trackdir[1] =-xinlocal[1] + xoutlocal[1];
  trackdir[2] =-xinlocal[2] + xoutlocal[2];
  alpha=-(xinlocal[0]*trackdir[0]+xinlocal[1]*trackdir[1])
    /(trackdir[0]*trackdir[0]+trackdir[1]*trackdir[1]);
  xlocal[0]=xinlocal[0]+trackdir[0]*alpha;  
  xlocal[1]=xinlocal[1]+trackdir[1]*alpha;
  xlocal[2]=xinlocal[2]+trackdir[2]*alpha;
  
  // Deal with tracks exiting the ends of the straws
  if (fabs(xlocal[2])>=75.0){
    int ring = getring_();
    float sign=(xoutlocal[2]>0)?1.:-1.;
    if (ring<=4 || (ring>=13 && ring<=16) || ring>=25){
      alpha=(sign*75.0-xinlocal[2])/trackdir[2];
      xlocal[0]=xinlocal[0]+trackdir[0]*alpha;  
      xlocal[1]=xinlocal[1]+trackdir[1]*alpha;
      xlocal[2]=sign*75.0;
    }
    else if (fabs(xlocal[2])>=75.415){
      alpha=(sign*75.415-xinlocal[2])/trackdir[2]; 
      xlocal[0]=xinlocal[0]+trackdir[0]*alpha;  
      xlocal[1]=xinlocal[1]+trackdir[1]*alpha;
      xlocal[2]=sign*75.415;
    }
  } 

  /* This will get called when the particle actually passes through
   * the wire volume itself. For these cases, we should set the 
   * location of the hit to be the point on the wire itself. Do
   * determine if this is what is happening, we check drout to
   * see if it is very close to the wire and drin to see if it is
   * close to the tube.
   *
   * For the other case, when drin is close to the wire, we assume
   * it is because it is emerging from the wire volume and
   * automatically ignore those hits by returning immediately.
   */
  if(drin < 0.0050)return; /* entering straw within 50 microns of wire. ignore */
  if(drin>(STRAW_RADIUS-0.0200) && drout<0.0050){
    /* We entered within 200 microns of the straw tube and left
     * within 50 microns of the wire. Assume the track passed through
     * the wire volume.
     */
    
    x[0] = xin[0];
    x[1] = xin[1];
    x[2] = xin[2];
    t = xin[3] * 1e9;
    xlocal[0] = xinlocal[0];
    xlocal[1] = xinlocal[1];
    xlocal[2] = xinlocal[2];
    
    /* For dx, we will just assume it is twice the distance from
     * the straw to wire.
     */
    dx[0] *= 2.0;
    dx[1] *= 2.0;
    dx[2] *= 2.0;
  }
  
  /* Distance of hit from center of wire */
  dradius = sqrt(xlocal[0]*xlocal[0] + xlocal[1]*xlocal[1]);

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
        int a = thisInputEvent->physicsEvents->in[0].reactions->in[0].vertices->in[0].products->mult;
         points->in[0].primary = (stack <= a);
         points->in[0].track = track;
         points->in[0].t = t;
         points->in[0].z = x[2];
         points->in[0].r = sqrt(x[0]*x[0] + x[1]*x[1]);
         points->in[0].phi = atan2(x[1],x[0]);
         points->in[0].dradius = dradius;
         points->in[0].px = pin[0]*pin[4];
         points->in[0].py = pin[1]*pin[4];
         points->in[0].pz = pin[2]*pin[4];
         points->in[0].dEdx = dEdx;
         points->in[0].ptype = ipart;
         points->mult = 1;
         cdc->cdcTruthPoints = points;
         pointCount++;
      }
   }

   /* post the hit to the hits tree, mark sector as hit */

   if (dEsum > 0)
   {
      int nhit;
      s_CdcStrawTruthHits_t* hits;
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
            straws->in[0].cdcStrawTruthHits = hits = make_s_CdcStrawTruthHits(MAX_HITS);
            cdc->cdcStraws = straws;
            strawCount++;
         }
         else
         {
            s_CentralDC_t* cdc = (s_CentralDC_t*) *twig;
            hits = cdc->cdcStraws->in[0].cdcStrawTruthHits;
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
				/* keep the earlier hit and discard the later one */
				/* Feb. 11, 2008 D. L. */
				if(hits->in[nhit].t>tdrift){
					hits->in[nhit].t = tdrift;
					hits->in[nhit].dE = dEsum;
					hits->in[nhit].d = dradius;
					hits->in[nhit].itrack = track;
					hits->in[nhit].ptype = ipart;
				}
			
           /* hits->in[nhit].t =
                    (hits->in[nhit].t * hits->in[nhit].dE + tdrift * dEsum)
                  / (hits->in[nhit].dE += dEsum);
				*/
         }
         else if (nhit < MAX_HITS)		/* create new hit */
         {
            hits->in[nhit].t = tdrift;
            hits->in[nhit].dE = dEsum;
	    hits->in[nhit].d = dradius;
	    hits->in[nhit].itrack = track;
	    hits->in[nhit].ptype = ipart;
	    
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

         s_CdcStrawTruthHits_t* hits = straws->in[straw].cdcStrawTruthHits;

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

