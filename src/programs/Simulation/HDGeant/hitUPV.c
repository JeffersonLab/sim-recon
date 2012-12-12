/*
 * hitUPV - registers hits for UPV - ao
 *
 *
 *	This is a part of the hits package for the
 *	HDGeant simulation program for Hall D.
 *
 *
 * changes: Wed Jun 20 13:19:56 EDT 2007 B. Zihlmann 
 *          add ipart to the function hitUpstreamEMveto
 *
 * Programmer's Notes:
 * -------------------
 * 1) In applying the attenuation to light propagating down to both ends
 *    of the modules, there has to be some point where the attenuation
 *    factor is 1.  I chose it to be the midplane, so that in the middle
 *    of the paddle both ends see the unattenuated E values.  Closer to
 *    either end, that end has a larger E value and the opposite end a
 *    lower E value than the actual deposition.
 * 2) In applying the propagation delay to light propagating down to the
 *    ends of the modules, there has to be some point where the timing
 *    offset is 0.  I chose it to be the midplane, so that for hits in
 *    the middle of the paddle the t values measure time-of-flight from
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
extern s_HDDM_t* thisInputEvent;

#define ATTEN_LENGTH	150.
#define C_EFFECTIVE	19.   /* This assumes a single linear fiber path */
#define THRESH_MEV      5.
#define TWO_HIT_RESOL	50.
#define MAX_HITS 	100

binTree_t* upstreamEMvetoTree = 0;
static int paddleCount = 0;
static int rowCount = 0;
static int showerCount = 0;


/* register hits during tracking (from gustep) */

void hitUpstreamEMveto (float xin[4], float xout[4],
			float pin[5], float pout[5], float dEsum,
			int track, int stack, int history, int ipart)
{
  float x[3], t;
  float xlocal[3];
  float xupv[3];
  float zeroHat[] = {0,0,0};
  int nhit;
  s_UpvLeftHits_t* leftHits;
  s_UpvRightHits_t* rightHits;

  x[0] = (xin[0] + xout[0])/2;
  x[1] = (xin[1] + xout[1])/2;
  x[2] = (xin[2] + xout[2])/2;
  t    = (xin[3] + xout[3])/2 * 1e9;
  transformCoord(x,"global",xlocal,"UPV");
  transformCoord(zeroHat,"local",xupv,"UPV");
  
  int layer = getlayer_wrapper_();
  int row = getrow_wrapper_();
  int column = getcolumn_wrapper_();
  /*
    'column' is not used in the current code. It distinguishes long 
    paddles (column=0) from short paddles to the left(column=1) or 
    right(column=2) of the beam hole. However, we assume that a pair 
    of short paddles is connected with a lightguide which has the light 
    propagation properties of a scintillator. In other words, a left and 
    right pair of short paddles form a long paddle which just happens 
    to have an insensitive to hits area in the middle but otherwise is identical
    to a normal long paddle. If we later change our minds and start treating 3 
    types of paddles differently, then 'column' is available for use.
  */

  float dxleft = xlocal[0];
  float dxright = -xlocal[0];
  float tleft  = t + dxleft/C_EFFECTIVE;
  float tright = t + dxright/C_EFFECTIVE;
  float dEleft  = dEsum * exp(-dxleft/ATTEN_LENGTH);
  float dEright = dEsum * exp(-dxright/ATTEN_LENGTH);
  float ycenter = (fabs(xupv[1]) < 1e-4) ? 0 : xupv[1];
  float zcenter = (fabs(xupv[2]) < 1e-4) ? 0 : xupv[2];

  /* post the hit to the truth tree */

  if ((history == 0) && (pin[3] > THRESH_MEV/1e3))
  {
    int mark = (1<<30) + showerCount;
    void** twig = getTwig(&upstreamEMvetoTree, mark);
    if (*twig == 0) {
      s_UpstreamEMveto_t* upv = *twig = make_s_UpstreamEMveto();
      s_UpvTruthShowers_t* showers = make_s_UpvTruthShowers(1);
        int a = thisInputEvent->physicsEvents->in[0].reactions->in[0].vertices->in[0].products->mult;
      showers->in[0].primary = (stack <= a);
      showers->in[0].track = track;
      showers->in[0].x = xin[0];
      showers->in[0].y = xin[1];
      showers->in[0].z = xin[2];
      showers->in[0].t = xin[3]*1e9;
      showers->in[0].px = pin[0]*pin[4];
      showers->in[0].py = pin[1]*pin[4];
      showers->in[0].pz = pin[2]*pin[4];
      showers->in[0].E = pin[3];
      showers->in[0].ptype = ipart;
      showers->mult = 1;
      upv->upvTruthShowers = showers;
      showerCount++;
    }
  }

  /* post the hit to the hits tree, mark upvPaddle as hit */

  if (dEsum > 0)
  {
    int mark = (layer<<16) + row;
    void** twig = getTwig(&upstreamEMvetoTree, mark);
    if (*twig == 0)
    {
      s_UpstreamEMveto_t* upv = *twig = make_s_UpstreamEMveto();
      s_UpvPaddles_t* paddles = make_s_UpvPaddles(1);
      paddles->mult = 1;
      paddles->in[0].row = row;
      paddles->in[0].layer = layer;
      leftHits = HDDM_NULL;
      rightHits = HDDM_NULL;
      if (column == 0 || column == 1)
      {
         paddles->in[0].upvLeftHits = leftHits
                                    = make_s_UpvLeftHits(MAX_HITS);
      }
      if (column == 0 || column == 1)
      {
         paddles->in[0].upvRightHits = rightHits
                                     = make_s_UpvRightHits(MAX_HITS);
      }
      upv->upvPaddles = paddles;
      paddleCount++;
    }
    else
    {
      s_UpstreamEMveto_t* upv = *twig;
      leftHits = upv->upvPaddles->in[0].upvLeftHits;
      rightHits = upv->upvPaddles->in[0].upvRightHits;
    }

    if (leftHits != HDDM_NULL)
    {
      for (nhit = 0; nhit < leftHits->mult; nhit++)
      {
        if (fabs(leftHits->in[nhit].t - tleft) < TWO_HIT_RESOL)
        {
          break;
        }
      }

      if (nhit < leftHits->mult)		/* merge with former hit */
      {
        leftHits->in[nhit].t =
	  (leftHits->in[nhit].t * leftHits->in[nhit].E + tleft * dEleft)
	  / (leftHits->in[nhit].E += dEleft);
      }
      else if (nhit < MAX_HITS)			/* create new hit */
      {
        leftHits->in[nhit].t =
	  (leftHits->in[nhit].t * leftHits->in[nhit].E + tleft * dEleft)
	  / (leftHits->in[nhit].E += dEleft);
        leftHits->mult++;
      }
      else
      {
        fprintf(stderr,"HDGeant error in hitUpstreamEMveto: ");
        fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
      }
    }

    if (rightHits != HDDM_NULL)
    {
      for (nhit = 0; nhit < rightHits->mult; nhit++)
      {
        if (fabs(rightHits->in[nhit].t - tright) < TWO_HIT_RESOL)
        {
	  break;
        }
      }
        
      if (nhit < rightHits->mult) 		/* merge with former hit */
      {
        rightHits->in[nhit].t =
	  (rightHits->in[nhit].t * rightHits->in[nhit].E + tright * dEright)
	  / (rightHits->in[nhit].E += dEright);
      }
      else if (nhit < MAX_HITS) 		/* create new hit */
      {
        rightHits->in[nhit].t = 
	  (rightHits->in[nhit].t * rightHits->in[nhit].E +  tright * dEright)
	  / (rightHits->in[nhit].E += dEright);
        rightHits->mult++;
      }
      else
      {
        fprintf(stderr,"HDGeant error in hitUpstreamEMveto: ");
        fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
      }
    }
  }
}

/* entry point from fortran */

void hitupstreamemveto_(float* xin, float* xout,
			float* pin, float* pout, float* dEsum,
			int* track, int* stack, int* history, int* ipart)
{
  hitUpstreamEMveto(xin,xout,pin,pout,*dEsum,*track,*stack,*history,*ipart);
}




/* pick and package the hits for shipping */

s_UpstreamEMveto_t* pickUpstreamEMveto ()
{
   s_UpstreamEMveto_t* box;
   s_UpstreamEMveto_t* item;

   if ((paddleCount == 0) && (rowCount == 0) && (showerCount == 0))
      return HDDM_NULL;

   box = make_s_UpstreamEMveto();
   box->upvPaddles = make_s_UpvPaddles(paddleCount);
   box->upvTruthShowers = make_s_UpvTruthShowers(showerCount);
   while (item = (s_UpstreamEMveto_t*) pickTwig(&upstreamEMvetoTree))
   {
      s_UpvPaddles_t* paddles = item->upvPaddles;
      int paddle;
      s_UpvTruthShowers_t* showers = item->upvTruthShowers;
      int shower;
    
      for (paddle=0; paddle < paddles->mult; ++paddle)
      {
         int m = box->upvPaddles->mult;
         int mok = 0;

         s_UpvLeftHits_t* leftHits = paddles->in[paddle].upvLeftHits;
         s_UpvRightHits_t* rightHits = paddles->in[paddle].upvRightHits;

       /* compress out the hits below threshold */
         int i,iok;
         for (iok=i=0; i < leftHits->mult; i++)
         {
            if (leftHits->in[i].E >= THRESH_MEV/1e3)
            {
               if (iok < i)
               {
                  leftHits->in[iok] = leftHits->in[i];
               }
               ++iok;
               ++mok;
            }
         }
         if (iok)
         {
            leftHits->mult = iok;
         }
         else if (leftHits != HDDM_NULL)
         {
             paddles->in[paddle].upvLeftHits = HDDM_NULL;
             FREE(leftHits);
         }

         for (iok=i=0; i < rightHits->mult; i++)
         {
            if (rightHits->in[i].E >= THRESH_MEV/1e3)
            {
              if (iok < i)
              {
                rightHits->in[iok] = rightHits->in[i];
              }
              ++iok;
              ++mok;
            }
         }
         if (iok)
         {
            rightHits->mult = iok;
         }
         else if (rightHits != HDDM_NULL)
         {
            paddles->in[0].upvRightHits = HDDM_NULL;
            FREE(rightHits);
         }

         if (mok)
         {
            box->upvPaddles->in[m] = paddles->in[paddle];
            box->upvPaddles->mult++;
         }
      }
      if (paddles != HDDM_NULL)
      {
         FREE(paddles);
      }

      for (shower=0; shower < showers->mult; ++shower)
      {
         int m = box->upvTruthShowers->mult++;
         box->upvTruthShowers->in[m] = showers->in[shower];
      }
      if (showers != HDDM_NULL)
      {
         FREE(showers);
      }
      FREE(item);
   }

   paddleCount = showerCount = 0;

   if ((box->upvPaddles != HDDM_NULL) &&
       (box->upvPaddles->mult == 0))
   {
      FREE(box->upvPaddles);
      box->upvPaddles = HDDM_NULL;
   }
   if ((box->upvTruthShowers != HDDM_NULL) &&
       (box->upvTruthShowers->mult == 0))
   {
      FREE(box->upvTruthShowers);
      box->upvTruthShowers = HDDM_NULL;
   }
   if ((box->upvPaddles->mult == 0) &&
       (box->upvTruthShowers->mult == 0))
   {
      FREE(box);
      box = HDDM_NULL;
   }
   return box;
}
