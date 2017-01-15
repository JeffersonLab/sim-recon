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
static int dircCount = 0;
static int dircpointCount = 0;

/* register truth points during tracking (from gustep) */
void hitDIRC(float xin[4], float xout[4], float pin[5], float pout[5],
		float dEsum, int track, int stack, int history, int ipart) {

    int itrack = (stack == 0)? gidGetId(track) : -1;

	// post to truth tree
	if ((history == 0) && (dEsum > 0)) {
		int mark = (1 << 25) + dircpointCount;
		void** twig = getTwig(&dircTree, mark);
		if (*twig == 0) {
			s_DIRC_t* dirc = *twig = make_s_DIRC();
			s_DircTruthPoints_t* points = make_s_DircTruthPoints(1);
			dirc->dircTruthPoints = points;
			int a = thisInputEvent->physicsEvents->in[0].reactions->in[0].vertices->in[0].products->mult;
			points->in[0].primary = (track <= a && stack == 0);
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
			points->in[0].trackID = make_s_TrackID();
			points->in[0].trackID->itrack = itrack;
			points->mult = 1;
			dircpointCount++;
		}
	}

	// post dirc hit
	if (dEsum < 0) {
		int mark = (1 << 20) + dircCount;
		void** twig = getTwig(&dircTree, mark);
		if (*twig == 0) {
			s_DIRC_t* dirc = *twig = make_s_DIRC();
			s_DircTruthHits_t* dircHits = make_s_DircTruthHits(1);
			dirc->dircTruthHits = dircHits;
			dircHits->in[0].x = xin[0];
			dircHits->in[0].y = xin[1];
			dircHits->in[0].z = xin[2];
			dircHits->in[0].t = xin[3] * 1e9;
			dircHits->in[0].E = pin[3];
			dircHits->mult = 1;
			dircCount++;
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

	if ((dircCount == 0) && (dircpointCount == 0)) {
		return HDDM_NULL ;
	}

	box = make_s_DIRC();
	// create DIRC hits
	box->dircTruthHits = make_s_DircTruthHits(dircCount);
	box->dircTruthPoints = make_s_DircTruthPoints(dircpointCount);

	while ((item = pickTwig(&dircTree))) {

		// pack DIRC hits
		s_DircTruthHits_t* dirchits = item->dircTruthHits;
		int dirchit;
		for (dirchit = 0; dirchit < dirchits->mult; ++dirchit) {
			int m = box->dircTruthHits->mult++;
			box->dircTruthHits->in[m] = dirchits->in[dirchit];
		}
		if (dirchits != HDDM_NULL) {
			FREE(dirchits);
		}
		// pack DIRC Truth points
		s_DircTruthPoints_t* dircpoints = item->dircTruthPoints;
        int last_track = -1;
        double last_t = 1e9;
		int dircpoint;
		for (dircpoint = 0; dircpoint < dircpoints->mult; ++dircpoint) {
            if (dircpoints->in[dircpoint].trackID->itrack > 0 &&
               (dircpoints->in[dircpoint].track != last_track ||
                fabs(dircpoints->in[dircpoint].t - last_t) > 0.1))
            {
               int m = box->dircTruthPoints->mult++;
               box->dircTruthPoints->in[m] = dircpoints->in[dircpoint];
               last_track = dircpoints->in[dircpoint].track;
               last_t = dircpoints->in[dircpoint].t;
            }
        }
		if (dircpoints != HDDM_NULL) {
			FREE(dircpoints);
		}
		FREE(item);
	}

	// clear DIRC hits and truth
	dircCount = dircpointCount = 0;
	if ((box->dircTruthHits != HDDM_NULL ) && (box->dircTruthHits->mult == 0)) {
		FREE(box->dircTruthHits);
		box->dircTruthHits = HDDM_NULL;
	}
	if ((box->dircTruthPoints != HDDM_NULL )
			&& (box->dircTruthPoints->mult == 0)) {
		FREE(box->dircTruthPoints);
		box->dircTruthPoints = HDDM_NULL;
	}
	if ((box->dircTruthHits->mult == 0) && (box->dircTruthPoints->mult == 0)) {
		FREE(box);
		box = HDDM_NULL;
	}
	return box;
}
