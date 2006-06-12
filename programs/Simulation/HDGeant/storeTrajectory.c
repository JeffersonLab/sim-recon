

#include <stdio.h>

#include "hddm_s.h"

extern s_HDDM_t* thisOutputEvent;

unsigned int Npoints;
unsigned int Maxpoints;
s_McTrajectoryPoint_t *traj_points = NULL;

/*---------------------
// cleartrajectories_
//--------------------*/
void cleartrajectories_(void)
{
	Npoints = 0;
	
	/* First time through, allocate buffer for trajectory points */
	if(traj_points==NULL){
		Maxpoints = 5000;
		traj_points = (s_McTrajectoryPoint_t*)malloc(Maxpoints*sizeof(s_McTrajectoryPoint_t));
	}
}

/*---------------------
// addtrajectorypoint_
//--------------------*/
void addtrajectorypoint_(float *VECT, float *TOFG, float *DESTEP
									,float *GETOT, int *ITRA, int *ISTAK, int *IPART)
{
	/* Only keep points from primary track */
	if(*ISTAK != 0)return;
	
	/* If buffer is full, print warning and return */
	static int Nwarns = 0;
	if(Npoints>=Maxpoints){
		if(Nwarns<10){
			fprintf(stderr,"%s:%d Too many trajectory points to store! Dropping some.\n");
			if(++Nwarns == 10)fprintf(stderr,"******** LAST WARNING!! *********\n");
		}
		return;
	}
	
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
	p->track = *ITRA;
}

/*---------------------
// pickMCTrajectory
//--------------------*/
s_McTrajectory_t* pickMCTrajectory(void)
{
	unsigned int i;
	
	if(Npoints==0)return HDDM_NULL;
	printf("%s:%d Npoints=%d\n",__FILE__,__LINE__,Npoints);
	
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



