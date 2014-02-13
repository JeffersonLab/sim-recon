

#ifndef __CERN_H__
#define __CERN_H__

#define real    float
#define integer int
#define logical int


/* GEANT Particle types */
enum geant_particles{
   ptype_none,
   ptype_gamma,
   ptype_positron,
   ptype_electron,
   ptype_neutrino,
   ptype_muon_plus,
   ptype_muon_minus,
   ptype_pion_zero,
   ptype_pion_plus,
   ptype_pion_minus,
   ptype_kaon_zero_long,
   ptype_kaon_plus,
   ptype_kaon_minus,
   ptype_neutron,
   ptype_proton,
   ptype_antiproton,
   ptype_kaon_zero_short,
   ptype_eta,

   // The following are particle types not defined in GEANT
   // This is a bit dangerous since GEANT has other particle
   // types defined which will use these values
   ptype_rho,
   ptype_omega,
   ptype_eta_prime
};

#define ELECTRON_MASS     0.00051100
#define MUON_MASS         0.10566
#define PI_CHARGED_MASS   0.13957
#define PI_ZERO_MASS      0.13498
#define KAON_CHARGED_MASS 0.49368
#define KAON_ZERO_MASS    0.49767
#define ETA_MASS          0.54745
#define PROTON_MASS       0.93827
#define NEUTRON_MASS      0.93957 
#define OMEGA_MASS        0.78194
#define RHO_MASS          0.770
#define ETA_PRIME_MASS    0.95778
#define PHI_MASS          1.019413
#define DEUTERON_MASS     1.877
#define TRITON_MASS       2.1573
#define LAMBDA_MASS       1.115683
#define SIGMA_ZERO_MASS   1.192642


/* Common Block Declarations */
typedef struct {
    integer jdigi, jdraw, jhead, jhits, jkine, jmate, jpart, jrotm, jrung, 
	    jset, jstak, jgstat, jtmed, jtrack, jvertx, jvolum, jxyz, jgpar, 
	    jgpar2, jsklt;
} gclink_t;
extern gclink_t gclink_;

typedef struct {
    integer nrecrz, nrget, nrsave, lrget[20], lrsave[20];
} gcrz1_t;
extern gcrz1_t gcrz1_;

typedef struct {
    char rztags[32];
} gcrz2_t;
extern gcrz2_t gcrz2_;

typedef struct {
    integer idebug, idemin, idemax, itest, idrun, idevt, ieorun, ieotri, 
	    ievent, iswit[10], ifinit[20], nevent, nrndm[2];
} gcflag_t;
extern gcflag_t gcflag_;

typedef struct {
    logical batch, nolog;
} gcflax_t;
extern gcflax_t gcflax_;

typedef struct {
    integer ikine;
    real pkine[10];
    integer itra, istak, ivert, ipart, itrtyp, napart[5];
    real xamass, charge, tlife, vert[3], pvert[4];
    integer ipaold;
} gckine_t;
extern gckine_t gckine_;

typedef struct {
    integer ihset, ihdet, iset, idet, idtype, nvname, numbv[20];
} gcsets_t;
extern gcsets_t gcsets_;

typedef struct {
    real timint, timend;
    integer itime, igdate, igtime;
} gctime_t;
extern gctime_t gctime_;

typedef struct {
    real vect[7], getot, gekin, vout[7];
    integer nmec, lmec[30], namec[30], nstep, maxnst;
    real destep, destel, safety, sleng, step, snext, sfield, tofg, gekrat, 
	    upwght;
    integer ignext, inwvol, istop, igauto, iekbin, ilsol, imull, ingoto, 
	    nldown, nlevin, nlsav, istory;
} gctrak_t;
extern gctrak_t gctrak_;

typedef struct {
    real polar[3];
    integer namec1[30];
} gctpol_t;
extern gctpol_t gctpol_;

typedef struct {
    integer kcase, ngkine;
    real gkin[100][5], tofd[100];
    integer iflgk[100];
} gcking_t;
extern gcking_t gcking_;

typedef struct {
    integer ngphot;
    real xphot[800][11]	;
} gckin2_t;
extern gckin2_t gckin2_;

typedef struct {
    real gpos[100][3];
} gckin3_t;
extern gckin3_t gckin3_;

typedef struct {
    integer np;
    real tecm, amass[18];
    integer kgenev;
} genin_t;
extern genin_t genin_;

typedef struct {
    real pcm[18][5], wt;
} genout_t;
extern genout_t genout_;

typedef struct {
    integer nlevel;
    char names[60];
    integer number[15];
    integer lvolum[15];
    integer lindex[15];
    integer inform,nlevmx;
    integer nldev[15];
    integer linmx[15];
    real gtran[45];
    real grmat[150];
    real gonly[15];
    real glx[3];
} gcvolu_t;
extern gcvolu_t gcvolu_;

typedef struct {
	int ipair;
	real spair,slpair,zintpa,steppa;
	int icomp;
	real scomp,slcomp,zintco,stepco;
	int iphot;
	real sphot,slphot,zintph,stepph;
	int ipfis;
	real spfis,slpfis,zintpf,steppf;
	int idray;
	real sdray,sldray,zintdr,stepdr;
	int ianni;
	real sanni,slanni,zintan,stepan;
	int ibrem;
	real sbrem,slbrem,zintbr,stepbr;
	int ihadr;
	real shadr,slhadr,zintha,stepha;
	int imunu;
	real smunu,slmunu,zintmu,stepmu;
	int idcay;
	real sdcay,slife ,sumlif,dphys1;
	int iloss;
	real sloss,soloss,stloss,dphys2;
	int imuls;
	real smuls,somuls,stmuls,dphys3;
	int irayl;
	real srayl,slrayl,zintra,stepra;
} gcphys_t;
extern gcphys_t gcphys_;


/* Function declarations */
#ifdef __cplusplus
extern "C" {
#endif
void hbname(int i,char *b,void *v,char *f);
void hbnamc(int i,char *b,char *v,char *f);
void hbook1(int n,char *N,int b,real m,real M,real v);
void hbook2(int n,char *N,int nx,real xm,real xM,int ny,real ym,real yM,real v);
void hbprof(int n,char *N,int nx,float xm,float xM,float ym,float yM,char *chopt);
void hbookn(int id,char *t,int NVAR,char *rz,int nw,char *tags);
void hbnt(int id,char *t,char *o);
void hlabel(int ID, int NLAB, char *CLAB, char*CHOPT);
void hcopy(int id1,int id2,char *ti);
int  hropen(int lun,char *nam,char*fnam,char*stat,int stor,int istat);
void hlimit(int size);
void gzebra(int size);
void hrput(int id, char *file, char *opt);
void hf1(int id, real data, real weight);
void hf2(int id,real data1,real data2,real weight);
void hfill(int id,float data1,float data2,float weight);
void hdelet(int id);
void hrin( int id, int icycle , int iofset );
void hfithn(int id, char *chfun, char *chopt, int np, real par[],
        real step[], real pmin[], real pmax[], real sigpar[], real *chi2);
void hunpak(int histo,real contents[],char choice[], int num);
void hidopt( int id, char *chopt);
void hpak(int histo,real contents[]);
void hrget( int id, char *chfile, char *chopt);
void hldir(char dir[],char flag[]);
void hmdir(char dir[],char flag[]);
void hcdir(char dir[],char flag[]);
void hrout(int num, int icycle, char*opt);
void hrend(char*filename);
void hreset(int no, char* opt);
void hfn(int id, real data[]);
void hfnt(int id);
void hfntb(int id,char *chblok);
void hprnt(int id);
void grndmq(int iseed1,int iseed2,int iseq,char *chopt);
void gsxyz(void);
void gdxyz(int k);
void gsvolu(char *name,char *type, int media , real PAR[],int NPAR,int *IVOL);
void gspos(char *VOL ,int id, char *MOTH ,real x, real y, real z, int IROT,char *CHONLY);
void gsrotm(int id,real th1,real ph1,real th2,real ph2,real th3,real ph3);
void gsdetv(char *CHSET,char *CHDET,int IDTYPE,int NWHI,int NWDI,int *ISET,int *IDET);
void gstpar(int ITMED, char* CHPAR, float PARVAL);
void gdeca2(real XMO,real XM1,real XM2,real**PCM);
void gdeca3(real XMO,real XM1,real XM2,real XM3,real**PCM);
void gmate(void);
void gsmixt(int IMATE, char* NAMATE, float* A, float* Z,float DENS,int NLMAT, float* WMAT);
void gsmate(int IMATE, char* CHNAMA, float A, float Z, float DENS, float RADL, float ABSL, float* UBUF, int NWBUF);
void gfmate(int IMATE, char* CHNAMA, float* A, float* Z, float* DENS, float* RADL, float* ABSL, float* UBUF, int* NWBUF);
void gpmate(int IMATE);
void gstmed(int ITMED, char* NATMED, int NMAT, int ISVOL, int IFIELD, float FIELDM, float TMAXFD,float STEMAX, float DEEMAX, float EPSIL, float STMIN, float* UBUF, int NWBUF );
void gptmed(int ITMED);
void gsatt(char *volu,char *attr,int ival);
void gprint(char *volu,int ival);
void ixupdwi(int ival);
void gsxyz(void);
void gdxyz(int ival);
void gsking(int ival);
void ffkey(char *KEY, void *VAR, int NVAR, char *TYPE);
void mninit(int a,int b, int c);
#ifdef __cplusplus
}
#endif

#endif /* __CERN_H__ */
