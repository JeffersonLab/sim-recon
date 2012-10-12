/*
 * hitRICH.c
 *
 *  Created on: Oct 11, 2012
 *      Author: yqiang
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <hddm_s.h>
#include <geant3.h>
#include <bintree.h>

extern s_HDDM_t* thisInputEvent;

binTree_t* richTree = 0;
static int richCount = 0;
static int richpointCount = 0;

/* register truth points during tracking (from gustep) */
void hitRich(float xin[4], float xout[4], float pin[5], float pout[5],
		float dEsum, int track, int stack, int history, int ipart) {
	float x[3], t;

	x[0] = (xin[0] + xout[0]) / 2;
	x[1] = (xin[1] + xout[1]) / 2;
	x[2] = (xin[2] + xout[2]) / 2;
	t = (xin[3] + xout[3]) / 2 * 1e9;

	// post to truth tree
	if ((history == 0) && (dEsum > 0)) {
		int mark = (1 << 25) + richpointCount;
		void** twig = getTwig(&richTree, mark);
		if (*twig == 0) {
			s_RICH_t* rich = *twig = make_s_RICH();
			s_RichTruthPoints_t* points = make_s_RichTruthPoints(1);
			rich->richTruthPoints = points;
			int a =
					thisInputEvent->physicsEvents->in[0].reactions->in[0].vertices->in[0].products->mult;
			points->in[0].primary = (stack <= a);
			points->in[0].track = track;
			points->in[0].x = xin[0];
			points->in[0].y = xin[1];
			points->in[0].z = xin[2];
			points->in[0].t = xin[3] * 1e9;
			points->in[0].px = pin[4] * pin[0];
			points->in[0].py = pin[4] * pin[1];
			points->in[0].pz = pin[4] * pin[2];
			points->in[0].E = pin[3];
			points->in[0].ptype = ipart;
			points->mult = 1;
			richpointCount++;
		}
	}

	// post rich hit
	if (dEsum < 0) {
		int mark = (1 << 20) + richCount;
		void** twig = getTwig(&richTree, mark);
		if (*twig == 0) {
			s_RICH_t* rich = *twig = make_s_RICH();
			s_RichHits_t* richHits = make_s_RichHits(1);
			rich->richHits = richHits;
			int a =
					thisInputEvent->physicsEvents->in[0].reactions->in[0].vertices->in[0].products->mult;
			richHits->in[0].x = xin[0];
			richHits->in[0].y = xin[1];
			richHits->in[0].z = xin[2];
			richHits->in[0].t = xin[3] * 1e9;
			richHits->mult = 1;
			richCount++;
		}
	}

}

/* entry points from fortran */
void hitrich_(float* xin, float* xout, float* pin, float* pout,
		float* dEsum, int* track, int* stack, int* history, int* ipart) {
	hitRich(xin, xout, pin, pout, *dEsum, *track, *stack, *history, *ipart);
}

/* pick and package the hits for shipping */

s_RICH_t* pickRich() {
	s_RICH_t* box;
	s_RICH_t* item;

	if ((richCount == 0) && (richpointCount == 0)) {
		return HDDM_NULL ;
	}

	box = make_s_RICH();
	// create RICH hits
	box->richHits = make_s_RichHits(richCount);
	box->richTruthPoints = make_s_RichTruthPoints(richpointCount);

	while (item = pickTwig(&richTree)) {

		// pack RICH hits
		s_RichHits_t* richhits = item->richHits;
		int richhit;
		for (richhit = 0; richhit < richhits->mult; ++richhit) {
			int m = box->richHits->mult++;
			box->richHits->in[m] = richhits->in[richhit];
		}
		if (richhits != HDDM_NULL) {
			FREE(richhits);
		}
		// pack RICH Truth points
		s_RichTruthPoints_t* richpoints = item->richTruthPoints;
		int richpoint;
		for (richpoint = 0; richpoint < richpoints->mult; ++richpoint) {
			int m = box->richTruthPoints->mult++;
			box->richTruthPoints->in[m] = richpoints->in[richpoint];
		}
		if (richpoints != HDDM_NULL) {
			FREE(richpoints);
		}
		FREE(item);
	}

	// clear RICH hits and truth
	richCount = richpointCount = 0;
	if ((box->richHits != HDDM_NULL ) && (box->richHits->mult == 0)) {
		FREE(box->richHits);
		box->richHits = HDDM_NULL;
	}
	if ((box->richTruthPoints != HDDM_NULL )
			&& (box->richTruthPoints->mult == 0)) {
		FREE(box->richTruthPoints);
		box->richTruthPoints = HDDM_NULL;
	}
	if ((box->richHits->mult == 0) && (box->richTruthPoints->mult == 0)) {
		FREE(box);
		box = HDDM_NULL;
	}
	return box;
}
