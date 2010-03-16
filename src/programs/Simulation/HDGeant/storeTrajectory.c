

#include <stdio.h>

#include "hddm_s.h"

extern s_HDDM_t* thisOutputEvent;

unsigned int Npoints;
unsigned int Maxpoints;
int last_track_num;
int last_stack_num;
int track_id;
s_McTrajectoryPoint_t *traj_points = NULL;

/*---------------------
// cleartrajectories_
//--------------------*/
void cleartrajectories_(void)
{
	Npoints = 0;
	
	/* First time through, allocate buffer for trajectory points */
	if(traj_points==NULL){
		Maxpoints = 100000;
		traj_points = (s_McTrajectoryPoint_t*)malloc(Maxpoints*sizeof(s_McTrajectoryPoint_t));
	}
	
	last_track_num = -1;
	last_stack_num = -1;
	track_id = 0;
}

/*---------------------
// addtrajectorypoint_
//--------------------*/
void addtrajectorypoint_(float *VECT, float *TOFG, float *DESTEP
									,float *GETOT, int *ITRA, int *ISTAK, int *IPART
									,float *RADL, float *STEP, int *NMEC, int* LMEC, int *storetraj)
{
	static int warned = 0;
	if(!warned){
		printf("\n\n    ** WARNING **\n\nYou have set the TRAJ flag to %d.\n", *storetraj);
		printf("This will significantly increase the output file size.\n\n\n");
		warned = 1;
	}
	
	/* We want to record a unique id for every particle in the event.
	/* The value in ITRA is the primary track's number from which this
	/* particle originated. The value in ISTAK is the stack poisition
	/* of the current particle (ISTAK==0 means it's the primary).
	/* Because the stack is reused during showers, there is no
	/* combination of ITRA and ISTAK that is guaranteed to be unique
	/* for a given particle in the event. Therefore, we need to keep
	/* track of this ourselves. Use the track_id global to do this by
	/* watching for changes to ITRA or increases in ISTAK.
	*/
	/* The variable "point_on_this_track" is used to decide whether or
	/* not to add this as a new track point, or to overwrite the previous
	/* track point if we're only keeping birth/death info.
	*/
	static int point_on_this_track = 0;
	if(last_stack_num<*ISTAK || last_track_num!=*ITRA){
		track_id++;
		point_on_this_track=0;
	}else{
		point_on_this_track++;
	}
	last_track_num = *ITRA;
	last_stack_num = *ISTAK;

	/* storetraj = 0  don't store trajectory info
	/* storetraj = 1  store birth and death points of primary tracks
	/* storetraj = 2  store birth and death points of all particles
	/* storetraj = 3  store full trajectory of primary tracks
	/* storetraj = 4  store full trajectory of primary tracks and birth/death points of secondaries
	/* storetraj = 5  store full trajectory for all particles
	*/
	
	int is_primary = (*ISTAK==0);
	int store_full_traj = 0;
	switch(*storetraj){
		case 0: return;
		case 1:
			if(!is_primary)return;
			break;
		case 2:
			break;
		case 3:
			if(!is_primary)return;
			store_full_traj = 1;
			break;
		case 4:
			if(is_primary)store_full_traj=1;
			break;
		case 5:
			store_full_traj = 1;
			break;
		default:
			return;
	}

	/* If buffer is full, print warning and return */
	static int Nwarns = 0;
	if(Npoints>=Maxpoints){
		if(Nwarns<10){
			fprintf(stderr,"%s:%d Too many trajectory points to store! Dropping some.\n",__FILE__,__LINE__);
			if(++Nwarns == 10)fprintf(stderr,"******** LAST WARNING!! *********\n");
		}
		return;
	}

	/* If we're only storing birth and death points, then backup and
	/* overwrite the last trajectory point on this track. If it is a
	/* different track, then don't overwrite the previous one.
	*/
	if(store_full_traj==0 && point_on_this_track>1 && Npoints>0)Npoints--;
	
	/* Finally, fill in the new trajectory point info */
	s_McTrajectoryPoint_t *p = &traj_points[Npoints++];
	
	p->E = *GETOT;
	p->dE = *DESTEP;
	p->part = *IPART;
	p->x = VECT[0];
	p->y = VECT[1];
	p->z = VECT[2];
	p->px = VECT[6]*VECT[3];
	p->py = VECT[6]*VECT[4];
	p->pz = VECT[6]*VECT[5];
	p->t = *TOFG;
	p->primary_track = *ITRA;
	p->track = track_id;
	p->radlen = *RADL;
	p->step = *STEP;
	p->mech = LMEC[*NMEC-1];
}

/*---------------------
// pickMCTrajectory
//--------------------*/
s_McTrajectory_t* pickMCTrajectory(void)
{
	unsigned int i;
		
	if(Npoints==0)return HDDM_NULL;
	
	s_McTrajectory_t* McTrajectory = make_s_McTrajectory();
	s_McTrajectoryPoints_t *points = make_s_McTrajectoryPoints(Npoints);
	McTrajectory->mcTrajectoryPoints = points;
	
	/* copying all of this isn't really the most efficient way to do this,
	  but if we're actually writing all of this stuff out, then it doesn't
	  need to be!
	*/
	points->mult = 0;
	for(i=0; i<Npoints; i++){
		points->in[i] = traj_points[i];
		points->mult++;
	}
	
	return McTrajectory;
}



