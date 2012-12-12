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


void hddsgeant3_wrapper_(void);
void md5geom_wrapper_(char *md5);
float guplsh_wrapper_(int *medi0, int *medi1);
void gufld_wrapper_(float *r, float *B);
void getoptical_wrapper_(int *imat, float *E, float *refl, float *abs1, float *rind, float *plsh, float *eff);

int getcolumn_wrapper_(void);
int getlayer_wrapper_(void);
int getmap_wrapper_(void);
int getmodule_wrapper_(void);
int getpackage_wrapper_(void);
int getplane_wrapper_(void);
int getring_wrapper_(void);
int getrow_wrapper_(void);
int getsector_wrapper_(void);


