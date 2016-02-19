/*
 * hitDIRC.c
 *
 *  Created on: Oct 11, 2012
 *      Author: yqiang
 *  Modified on June 22, 2015: changed RICH -> DIRC and remove CERE
 *      Author: jrsteven
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <HDDM/hddm_s.h>
#include <geant3.h>
#include <bintree.h>
#include <gid_map.h>

extern s_HDDM_t* thisInputEvent;

binTree_t* dircTree = 0;
static int dircTruthHitCount = 0;
static int dircTruthPointCount = 0;

/* register truth points during tracking (from gustep) */
void hitDIRC(float xin[4], float xout[4], float pin[5], float pout[5],
		float dEsum, int track, int stack, int history, int ipart) {
	//float x[3], t;

	//x[0] = (xin[0] + xout[0]) / 2;
	//x[1] = (xin[1] + xout[1]) / 2;
	//x[2] = (xin[2] + xout[2]) / 2;
	//t = (xin[3] + xout[3]) / 2 * 1e9;

	//printf("%f\n",dEsum);

	// post dirc truth point
	if ((history == 0) && (dEsum > 0)) {
		int mark = (1 << 25) + dircTruthPointCount;
		void** twig = getTwig(&dircTree, mark);
		if (*twig == 0) {
			s_DIRC_t* dirc = *twig = make_s_DIRC();
			s_DircTruthPoints_t* truthPoints = make_s_DircTruthPoints(1);
			dirc->dircTruthPoints = truthPoints;
			int a = thisInputEvent->physicsEvents->in[0].reactions->in[0].vertices->in[0].products->mult;
			truthPoints->in[0].primary = (stack <= a);
			truthPoints->in[0].track = track;
			truthPoints->in[0].x = xin[0];
			truthPoints->in[0].y = xin[1];
			truthPoints->in[0].z = xin[2];
			truthPoints->in[0].t = xin[3] * 1e9;
			truthPoints->in[0].px = pin[4] * pin[0];
			truthPoints->in[0].py = pin[4] * pin[1];
			truthPoints->in[0].pz = pin[4] * pin[2];
			truthPoints->in[0].E = pin[3];
			truthPoints->in[0].ptype = ipart;
			truthPoints->in[0].trackID = make_s_TrackID();
			truthPoints->in[0].trackID->itrack = gidGetId(track);
			truthPoints->mult = 1;
			dircTruthPointCount++;
		}
	}

	// post dirc truth hit
	if (dEsum < 0) {
	        //printf("found truth hit\n");
		int mark = (1 << 20) + dircTruthHitCount;
		void** twig = getTwig(&dircTree, mark);
		if (*twig == 0) {
		        s_DIRC_t* dirc = *twig = make_s_DIRC();
			s_DircTruthHits_t* truthHits = make_s_DircTruthHits(1);
			dirc->dircTruthHits = truthHits;
			truthHits->in[0].x = xin[0];
			truthHits->in[0].y = xin[1];
			truthHits->in[0].z = xin[2];
			truthHits->in[0].t = xin[3] * 1e9;
			truthHits->in[0].px = pin[4] * pin[0];
			truthHits->in[0].py = pin[4] * pin[1];
			truthHits->in[0].pz = pin[4] * pin[2];
			truthHits->in[0].E = pin[3];
			truthHits->in[0].track = track;
			truthHits->mult = 1;
			dircTruthHitCount++;
		}
	}

}

/* entry points from fortran */
void hitdirc_(float* xin, float* xout, float* pin, float* pout,
		float* dEsum, int* track, int* stack, int* history, int* ipart) {
	hitDIRC(xin, xout, pin, pout, *dEsum, *track, *stack, *history, *ipart);
}

/* pick and package the hits for shipping */

s_DIRC_t* pickDirc() {
	s_DIRC_t* box;
	s_DIRC_t* item;

	if ((dircTruthHitCount == 0) && (dircTruthPointCount == 0)) {
		return HDDM_NULL ;
	}

	box = make_s_DIRC();
	// create DIRC hits
	box->dircTruthHits = make_s_DircTruthHits(dircTruthHitCount);
	box->dircTruthPoints = make_s_DircTruthPoints(dircTruthPointCount);

	while ((item = pickTwig(&dircTree))) {

		// pack DIRC hits
		s_DircTruthHits_t* dircTruthHits = item->dircTruthHits;
		int dircTruthHit;
		for (dircTruthHit = 0; dircTruthHit < dircTruthHits->mult; ++dircTruthHit) {
			int m = box->dircTruthHits->mult++;
			box->dircTruthHits->in[m] = dircTruthHits->in[dircTruthHit];
		}
		if (dircTruthHits != HDDM_NULL) {
			FREE(dircTruthHits);
		}
		// pack DIRC Truth points
		s_DircTruthPoints_t* dircTruthPoints = item->dircTruthPoints;
		int dircTruthPoint;
		for (dircTruthPoint = 0; dircTruthPoint < dircTruthPoints->mult; ++dircTruthPoint) {
			int m = box->dircTruthPoints->mult++;
			box->dircTruthPoints->in[m] = dircTruthPoints->in[dircTruthPoint];
		}
		if (dircTruthPoints != HDDM_NULL) {
			FREE(dircTruthPoints);
		}
		FREE(item);
	}

	// clear DIRC hits and truth
	dircTruthHitCount = dircTruthPointCount = 0;
	if ((box->dircTruthHits != HDDM_NULL ) && (box->dircTruthHits->mult == 0)) {
		FREE(box->dircTruthHits);
		box->dircTruthHits = HDDM_NULL;
	}
	if ((box->dircTruthPoints != HDDM_NULL ) && (box->dircTruthPoints->mult == 0)) {
		FREE(box->dircTruthPoints);
		box->dircTruthPoints = HDDM_NULL;
	}
	if ((box->dircTruthHits->mult == 0) && (box->dircTruthPoints->mult == 0)) {
		FREE(box);
		box = HDDM_NULL;
	}
	return box;
}
