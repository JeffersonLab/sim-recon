/*
 * hitUPV - registers hits for UPV - ao
 *
 *
 *	This is a part of the hits package for the
 *	HDGeant simulation program for Hall D.
 *
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

#define ATTEN_LENGTH	150.
#define C_EFFECTIVE	19.   /* Index of refraction for UPV is 1.58 */
#define THRESH_ROW_MEV	5.
#define THRESH_PADDLE_MEV 0.5
#define TWO_HIT_RESOL	50.
#define MAX_HITS 	100

binTree_t* upstreamEMvetoTree = 0;
static int paddleCount = 0;
static int rowCount = 0;
static int showerCount = 0;


/* register hits during tracking (from gustep) */

void hitUpstreamEMveto (float xin[4], float xout[4],
			float pin[5], float pout[5], float dEsum,
			int track, int stack)
{
  float x[3], t;
  float xlocal[3];
  float xupv[3];
  float zeroHat[] = {0,0,0};
  int nshot;
  s_Showers_t* leftshots;
  s_Showers_t* rightshots;

  if (dEsum == 0) return;              /* only seen if it deposits energy */

  x[0] = (xin[0] + xout[0])/2;
  x[1] = (xin[1] + xout[1])/2;
  x[2] = (xin[2] + xout[2])/2;
  t    = (xin[3] + xout[3])/2 * 1e9;
  transformCoord(x,"global",xlocal,"UPV");
  transformCoord(zeroHat,"local",xupv,"UPV");
  
  int layer = getlayer_();
  int row = getrow_();
  int column = getcolumn_();
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

  /* post the hit to the hits tree, mark upvPaddle as hit */
  /* upvPaddles are used in MC studies only. They are not read separately in the data */
  {
    int mark = (layer << 16) + row;
    void** twig = getTwig(&upstreamEMvetoTree, mark);
    if (*twig == 0) {
      s_UpstreamEMveto_t* upv = *twig = make_s_UpstreamEMveto();
      upv->upvPaddles = make_s_UpvPaddles(1);
      upv->upvPaddles->mult = 1;
      upv->upvPaddles->in[0].y = ycenter;
      upv->upvPaddles->in[0].z = zcenter;
      upv->upvPaddles->in[0].upvLeft = make_s_UpvLeft();
      upv->upvPaddles->in[0].upvRight = make_s_UpvRight();
      upv->upvPaddles->in[0].upvLeft->showers = leftshots = make_s_Showers(MAX_HITS);
      upv->upvPaddles->in[0].upvRight->showers = rightshots = make_s_Showers(MAX_HITS);
      paddleCount++;
    }
    else {
      s_UpstreamEMveto_t* upv = *twig;
      leftshots = upv->upvPaddles->in[0].upvLeft->showers;
      rightshots = upv->upvPaddles->in[0].upvRight->showers;
    }

    for (nshot = 0; nshot < leftshots->mult; nshot++)
      if (fabs(leftshots->in[nshot].t - tleft) < TWO_HIT_RESOL)
	break;

    if (nshot < leftshots->mult) {		/* merge with former hit */
      leftshots->in[nshot].t =
	(leftshots->in[nshot].t * leftshots->in[nshot].E + tleft * dEleft)
	/ (leftshots->in[nshot].E += dEleft);
    }
    else if (nshot < MAX_HITS) {		/* create new hit */
      leftshots->in[nshot].t =
	(leftshots->in[nshot].t * leftshots->in[nshot].E + tleft * dEleft)
	/ (leftshots->in[nshot].E += dEleft);
      leftshots->mult++;
    }
    else {
      fprintf(stderr,"HDGeant error in hitUpstreamEMveto: ");
      fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
    }

    for (nshot = 0; nshot < rightshots->mult; nshot++)
      if (fabs(rightshots->in[nshot].t - tright) < TWO_HIT_RESOL)
	break;
        
    if (nshot < rightshots->mult) {		/* merge with former hit */
      rightshots->in[nshot].t =
	(rightshots->in[nshot].t * rightshots->in[nshot].E + tright * dEright)
	/ (rightshots->in[nshot].E += dEright);
    }
    else if (nshot < MAX_HITS) {		/* create new hit */
      rightshots->in[nshot].t = 
	(rightshots->in[nshot].t * rightshots->in[nshot].E +  tright * dEright)
	/ (rightshots->in[nshot].E += dEright);
      rightshots->mult++;
    }
    else {
      fprintf(stderr,"HDGeant error in hitUpstreamEMveto: ");
      fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
    }
  }

  /* post the hit to the hits tree, mark upvRow as hit */
  {
    int mark = row;
    void** twig = getTwig(&upstreamEMvetoTree, mark);
    if (*twig == 0) {
      s_UpstreamEMveto_t* upv = *twig = make_s_UpstreamEMveto();
      upv->upvRows = make_s_UpvRows(1);
      upv->upvRows->mult = 1;
      upv->upvRows->in[0].y = ycenter;
      upv->upvRows->in[0].upvLeft = make_s_UpvLeft();
      upv->upvRows->in[0].upvRight = make_s_UpvRight();
      upv->upvRows->in[0].upvLeft->showers = leftshots = make_s_Showers(MAX_HITS);
      upv->upvRows->in[0].upvRight->showers = rightshots = make_s_Showers(MAX_HITS);
      rowCount++;
    }
    else {
      s_UpstreamEMveto_t* upv = *twig;
      leftshots = upv->upvRows->in[0].upvLeft->showers;
      rightshots = upv->upvRows->in[0].upvRight->showers;
    }

    for (nshot = 0; nshot < leftshots->mult; nshot++)
      if (fabs(leftshots->in[nshot].t - tleft) < TWO_HIT_RESOL)
	break;

    if (nshot < leftshots->mult) {		/* merge with former hit */
      leftshots->in[nshot].t =
	(leftshots->in[nshot].t * leftshots->in[nshot].E + tleft * dEleft)
	/ (leftshots->in[nshot].E += dEleft);
    }
    else if (nshot < MAX_HITS) {		/* create new hit */
      leftshots->in[nshot].t =
	(leftshots->in[nshot].t * leftshots->in[nshot].E + tleft * dEleft)
	/ (leftshots->in[nshot].E += dEleft);
      leftshots->mult++;
    }
    else {
      fprintf(stderr,"HDGeant error in hitUpstreamEMveto: ");
      fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
    }

    for (nshot = 0; nshot < rightshots->mult; nshot++)
      if (fabs(rightshots->in[nshot].t - tright) < TWO_HIT_RESOL)
	break;

    if (nshot < rightshots->mult) {		/* merge with former hit */
      rightshots->in[nshot].t =
	(rightshots->in[nshot].t * rightshots->in[nshot].E + tright * dEright)
	/ (rightshots->in[nshot].E += dEright);
    }
    else if (nshot < MAX_HITS) {		/* create new hit */
      rightshots->in[nshot].t = 
	(rightshots->in[nshot].t * rightshots->in[nshot].E +  tright * dEright)
	/ (rightshots->in[nshot].E += dEright);
      rightshots->mult++;
    }
    else {
      fprintf(stderr,"HDGeant error in hitUpstreamEMveto: ");
      fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
    }
  }

  /* post the hit to the truth tree, once per primary track */
  {
    int mark = (track << 24);
    void** twig = getTwig(&upstreamEMvetoTree, mark);
    if (*twig == 0) {
      s_UpstreamEMveto_t* upv = *twig = make_s_UpstreamEMveto();
      s_UpvShowers_t* showers = make_s_UpvShowers(1);
      upv->upvShowers = showers;
      showers->in[0].primary = (stack == 0);
      showers->in[0].track = track;
      showers->in[0].x = x[0];
      showers->in[0].y = x[1];
      showers->in[0].z = x[2];
      showers->in[0].t = t;
      showers->in[0].E = dEsum;
      showers->mult = 1;
      showerCount++;
    }
    else {
      s_UpvShowers_t* showers = ((s_UpstreamEMveto_t*) *twig)->upvShowers;
      showers->in[0].x = (showers->in[0].x * showers->in[0].E + x[0]*dEsum)
	/ (showers->in[0].E + dEsum);
      showers->in[0].y = (showers->in[0].y * showers->in[0].E + x[1]*dEsum)
	/ (showers->in[0].E + dEsum);
      showers->in[0].z = (showers->in[0].z * showers->in[0].E + x[2]*dEsum)
	/ (showers->in[0].E + dEsum);
      showers->in[0].t = (showers->in[0].t * showers->in[0].E + t*dEsum)
	/ (showers->in[0].E += dEsum);
    }
  }
}

/* entry point from fortran */

void hitupstreamemveto_(float* xin, float* xout,
			float* pin, float* pout, float* dEsum,
			int* track, int* stack)
{
  hitUpstreamEMveto(xin,xout,pin,pout,*dEsum,*track,*stack);
}




/* pick and package the hits for shipping */

s_UpstreamEMveto_t* pickUpstreamEMveto ()
{
  s_UpstreamEMveto_t* box;
  s_UpstreamEMveto_t* item;

  if ((paddleCount == 0) && (rowCount == 0) && (showerCount == 0))
    return 0;

  box = make_s_UpstreamEMveto();
  box->upvRows = make_s_UpvRows(paddleCount);
  box->upvPaddles = make_s_UpvPaddles(paddleCount);
  box->upvShowers = make_s_UpvShowers(showerCount);
  while (item = (s_UpstreamEMveto_t*) pickTwig(&upstreamEMvetoTree)) {
    if (item->upvPaddles) {
      float Eleft  = item->upvPaddles->in[0].upvLeft->showers->in[0].E;
      float Eright = item->upvPaddles->in[0].upvRight->showers->in[0].E;
      if (Eleft + Eright > THRESH_PADDLE_MEV/1000) {
	int m = box->upvPaddles->mult++;
	box->upvPaddles->in[m] = item->upvPaddles->in[0];
      }
      else {
	FREE(item->upvPaddles->in[0].upvLeft->showers);
	FREE(item->upvPaddles->in[0].upvRight->showers);
	FREE(item->upvPaddles->in[0].upvLeft);
	FREE(item->upvPaddles->in[0].upvRight);
      }
      FREE(item->upvPaddles);
    }
    else if (item->upvRows) {
      float Eleft  = item->upvRows->in[0].upvLeft->showers->in[0].E;
      float Eright = item->upvRows->in[0].upvRight->showers->in[0].E;
      if (Eleft + Eright > THRESH_ROW_MEV/1000) {
	int m = box->upvRows->mult++;
	box->upvRows->in[m] = item->upvRows->in[0];
      }
      else {
	FREE(item->upvRows->in[0].upvLeft->showers);
	FREE(item->upvRows->in[0].upvRight->showers);
	FREE(item->upvRows->in[0].upvLeft);
	FREE(item->upvRows->in[0].upvRight);
      }
      FREE(item->upvRows);
    }
    else if (item->upvShowers) {
      int m = box->upvShowers->mult++;
      box->upvShowers->in[m] = item->upvShowers->in[0];
      FREE(item->upvShowers);
    }
    FREE(item);
  }
  paddleCount = rowCount = showerCount = 0;
  return box;
}
