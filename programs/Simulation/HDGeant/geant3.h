void gsvert_(float v[3], int* ntbeam, int* nttarg,
             float ubuf[], int* nubuf, int* nvtx);

void gskine_(float p[3], Particle_t* kind, int* nvtx, 
             float ubuf[], int* nubuf, int* ntrk);

void grndm_(float v[], int* len);

/* convenience interface function for gmtod_ and gdtom_ */

#define transformCoord(xin,sin,xout,sout) \
   transformcoord_(xin,sin,xout,sout,strlen(sin),strlen(sout))
