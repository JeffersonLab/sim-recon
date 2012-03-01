void gsvert_(float vert[3], int* ntbeam, int* nttarg,
             float ubuf[], int* nubuf, int* nvtx);

void gfvert_(int* nvtx, float vert[3], int* ntbeam, int* nttarg,
             float* tofg, float ubuf[], int* nubuf);

void gskine_(float plab[3], Particle_t* ipart, int* nv, 
             float ubuf[], int* nubuf, int* nt);

void gfkine_(int* itra, float vert[3], float pvert[3], Particle_t* ipart,
             int* nvert, float ubuf[], int* nubuf);

void grndm_(float v[], int* len);

/* convenience interface function for gmtod_ and gdtom_ */

#define transformCoord(xin,sin,xout,sout) \
   transformcoord_(xin,sin,xout,sout,strlen(sin),strlen(sout))


/* Type declarations to avoid "implicit function declaration" errors */
void transformcoord_(float* xin, const char* sin, float* xout, const char* sout, int, int);
int getsector_(void);
int getlayer_(void);
int getmodule_(void);
int getrow_(void);
int getcolumn_(void);
int getplane_(void);
int getring_(void);
