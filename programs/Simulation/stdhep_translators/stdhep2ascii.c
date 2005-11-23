/************************************************
 * stdhep2ascii.c 
 * This program converts the StdHep format to the ascii 
 * format of genr8    See http://www-pat.fnal.gov/stdhep.html
 *
 * Paul Eugenio
 * Carnegie Mellon University
 * 24 Sept 98
 **********************************************/

#include<stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


#include<particleType.h>
#include<stdhep.h>
#include<stdlun.h>
#include<stdcnt.h>

/*
 * #include <stdlib.h> See it below.
 */

/******************* From stdhep.h **************************** 
* Basic COMMON block from STDHEP: the HEPEVT COMMON block 
* See product StDhep
*
*  note that to avoid alignment problems, structures and common blocks
*    should be in the order: double precision, real, integer.
***************************************************************
#define NMXHEP 4000
extern struct hepevt {
  int nevhep;             /* The event number *
  int nhep;               /* The number of entries in this event *
  int isthep[NMXHEP];     /* The Particle id *
  int idhep[NMXHEP];      /* The particle id *
  int jmohep[NMXHEP][2];    /* The position of the mother particle *
  int jdahep[NMXHEP][2];    /* Position of the first daughter... *
  double phep[NMXHEP][5];    /* 4-Momentum, mass *
  double vhep[NMXHEP][4];    /* Vertex information *
} hepevt_;
************************************************************/

#define MAX_PARTS 2000
/******************
 *  Local Structure
 ******************/
typedef struct { double x,y,z; } dvector3_t;
typedef struct { double t; dvector3_t space; } dvector4_t;

typedef struct {
  int statusCode; /* 1 is a final state particle */
  int pid; /* this use the PDG MC numbering scheme */
  double mass;
  dvector4_t p;
}mc_part_t;

typedef struct{
  int runNo;
  int eventNo;
  int nparts;
  mc_part_t part[MAX_PARTS];
}mc_evt_t;

/***************
 *  GLOBALS
 ***************/
int Debug=0;


/********************************
 * Prototypes for StdHep functions.
 *********************************/

int StdHepXdrReadInit(char *fileName,int ntries, int istream);
int StdHepXdrRead(int *ilbl,int istream);
int StdHepXdrEnd(int istream);


/********************************
 * Local Prototypes 
 *********************************/
int PrintUsage(char *processName);
int fill_mc_parts(mc_evt_t *mc_evt);
int write_mc_parts(FILE *fp, mc_evt_t *mc_evt);
int write_gamp_parts(FILE *fp, mc_evt_t *mc_evt,int beamE);
int getCharge(int PDGpid);
int  gampID(int id);

/****************
 * main()
 **************/
 /* 
  * MAX_PARTS sets the max numbers of particles in an event
  */

int main(int argc,char **argv)
{
  char *argptr;
  int i,ntries=0,ret,gotANevent,written=0;
  int ilbl=1,istream=0,nparts=0;/*I guessed at istream */ 
  int kludge_beam =-1;
  char *evtfile ="default.evt";
  char outputfile[40];
  FILE *outputfp=stdout;
  mc_evt_t mc_evt;



  if (argc == 1){
    PrintUsage(argv[0]);
    exit (0);
  }
   else { /* good start */
    for (i=1; i<argc; i++) {
      argptr = argv[i];
      if ((*argptr == '-') && (strlen(argptr) > 1)) {
        argptr++;
        switch (*argptr) {
        case 'd':
          Debug =2;
          break;
	  case 'n':
          nparts=atoi(++argptr);
          break;
	  case 'g':
          kludge_beam=atoi(++argptr);
	  fprintf(stderr,"Writing gamp format kludge beam of %d\n",kludge_beam);
          break;
        case 'N':
          ntries=atoi(++argptr);
          break;
	case 'o':
          sprintf(outputfile,"%s",++argptr);
	  if(!(outputfp = fopen(outputfile,"w"))){
            fprintf(stderr,"Fail to open output file!\n");
            exit(-1);
          }
          break;
	case 'i':
	  evtfile= ++argptr;
          break;
	case 'h':
          PrintUsage(argv[0]);
          exit(0);
          break;
	default:
          fprintf(stderr,"Unrecognized argument -%s\n\n",argptr);
          PrintUsage(argv[0]);    
          exit(-1);
          break;
        }
      }
    }

    /* main code */

    /*
     * Open and init a stdhep file.
     */
    
    if(ret=StdHepXdrReadInit(evtfile,ntries, istream)){
      /* some error must have occured */
      fprintf(stderr,"err:StdHepXdrReadInit ret=%d\n",ret);
      exit(-1);
    }
    
    /*
     * Fill the mc_part_t structure from the HEPEVENT structure.
     */
      
    for(i=0;i<ntries;i++){
      if(ret=StdHepXdrRead(&ilbl,istream)){/* some error must have occured */
	fprintf(stderr,"err:StdHepXdrRead ret=%d\n",ret);
	exit(-1); 
      }
      fill_mc_parts(&mc_evt);
      if(kludge_beam>0)
	write_gamp_parts(outputfp,&mc_evt,kludge_beam);
      else
	write_mc_parts(outputfp,&mc_evt);
      if(!(++written %100))
	fprintf(stderr,"McFast events Read: %d\r",written);  
    }

    fprintf(stderr,"\nTotal McFast events Read: %d\n",written);  
  
    /*
     * Close the StdHep file and exit.
     */

     StdHepXdrEnd(istream);
     exit(0);

   } /* end of else /* good start */
} /* end of main */

/*********************
 *
 * fill_mc_parts  
 * 
 *******************/


int fill_mc_parts(mc_evt_t *mc_evt){
  int i;
  int runNo=-9000;
  
  /* 
   * fill header info 
   */
  mc_evt->runNo=runNo;
  mc_evt->eventNo=hepevt_.nevhep;
  mc_evt->nparts=hepevt_.nhep;
  
  /* 
   * now loop over the particle in the event 
   */
  for(i=0;i<mc_evt->nparts;i++){
    mc_evt->part[i].statusCode=hepevt_.isthep[i];
    mc_evt->part[i].pid=hepevt_.idhep[i];
    mc_evt->part[i].p.space.x=hepevt_.phep[i][0];
    mc_evt->part[i].p.space.y=hepevt_.phep[i][1];
    mc_evt->part[i].p.space.z=hepevt_.phep[i][2];
    mc_evt->part[i].p.t=hepevt_.phep[i][3];
    mc_evt->part[i].mass=hepevt_.phep[i][4];
  }
  return 1;
}


/********************
 * write_gamp_parts
 *******************/
int write_gamp_parts(FILE *fp, mc_evt_t *mc_evt, int photon_beam){
  int i;
  /*
   * Write the number of particles (including the beam)
   */
  fprintf(fp,"%d \n", mc_evt->nparts+1);
  /*
   * Write particle info (id charge p.x p.y p.z p.t)
   */
  /* kludged beam */
  fprintf(fp,"1 0 0 0 %d %d \n", photon_beam,photon_beam);
  /* now write out nparts */
  for(i=0;i<mc_evt->nparts;i++){ 
    fprintf(fp,"%d %d %lf %lf %lf %lf \n",
	    gampID(mc_evt->part[i].pid),getCharge( mc_evt->part[i].pid),
	    mc_evt->part[i].p.space.x,
	    mc_evt->part[i].p.space.y,
	    mc_evt->part[i].p.space.z,
	    mc_evt->part[i].p.t);
  }
}

/********************
 * write_mc_parts
 *******************/
int write_mc_parts(FILE *fp, mc_evt_t *mc_evt){
  int i;
  
  
  /* 
   * write header info 
   */
  fprintf(fp,"%d  %d\n", mc_evt->runNo, mc_evt->eventNo);
  
  /* 
   * now loop over the particle in the event 
   */
  for(i=0;i<mc_evt->nparts;i++){
    fprintf(fp,"%d %d %lf \n",(i+1), 
	    gampID(mc_evt->part[i].pid), mc_evt->part[i].mass);
    fprintf(fp,"   %d %lf %lf %lf %lf\n",
	    getCharge( mc_evt->part[i].pid),
	    /* mc_evt->part[i].pid/abs(mc_evt->part[i].pid),*/
	    mc_evt->part[i].p.space.x,
	    mc_evt->part[i].p.space.y,
	    mc_evt->part[i].p.space.z,
	    mc_evt->part[i].p.t);
  }
  return 1;
}


/********************
 * getCharge
 *******************/
int getCharge(int PDGpid){
  int charge;
  /* 
   * See 1998 PDG Rev. MC particle numbering scheme 
   */
  typedef enum {
    gamma=22,
    pizero=111,
    eta=221,
    etaprime=331,
    Kshort=310,
    Klong=130,
    neutron=2112,
    Lambda=3122
  } PGDneutralPID_t;

  charge = PDGpid/abs(PDGpid);
  /* 
   * now check for neutrals 
   */
  switch((PGDneutralPID_t) PDGpid){ /* See PDG Rev. MC particle numbering scheme */
  case gamma:
  case pizero:
  case eta:
  case etaprime:
  case Kshort:
  case Klong:
  case neutron:
  case -neutron:
  case Lambda:
  case -Lambda:
    charge =0;
  default:
    break;
  }

  return charge;  

}
/********************
 * gampID(int pdgID)
 * The convention is that
 * followed in StdHep 4.02
 *******************/
int  gampID(int id){
  Particle_t p=Unknown;
  switch (id) {
  case 0:
    p=Unknown;
    break;
  case 22:
    p=Gamma;
    break;
  case -11:
    p=Positron;
    break;
  case 11:
    p=Electron;
    break;
  case 12:
    p=Neutrino;
    break;
  case -13:
    p=MuonPlus;
    break;
  case 13:
    p=MuonMinus;
    break;
  case 111:
    p=Pi0;
    break;
  case 211:
    p=PiPlus;
    break;
  case -211:
    p=PiMinus;
    break;
  case 130:
    p=KLong;
    break;
  case 321:
    p=KPlus;
    break;
  case -321:
    p=KMinus;
    break;
  case 2112:
    p=Neutron;
    break;
  case 2212:
    p=Proton;
    break;
  case -2212:
    p=AntiProton;
    break;
  case 310:
    p=KShort;
    break;
  case 221:
    p=Eta;
    break;
  case 3122:
    p=Lambda;
    break;
  case 3222:
    p=SigmaPlus;
    break;
  case 3212:
    p=Sigma0;
    break;
  case 3112:
    p=SigmaMinus;
    break;
  case 3322:
    p=Xi0;
  case 3312:
    p=XiMinus;
    break;
  case 3334:
    p=OmegaMinus;
    break;
  case -2112:
    p=AntiNeutron;
    break;
  case -3122:
    p=AntiLambda;
    break;
  case -3112:
    p=AntiSigmaMinus;
    break;
  case -3212:
    p=AntiSigma0;
    break;
  case -3222:
    p=AntiSigmaPlus;
    break;
  case -3322:
    p=AntiXi0;
    break;
  case -3312:
    p=AntiXiPlus;
    break;
  case -3334:
    p=AntiOmegaPlus;
    break;
  case 113:
    p=Rho0;
    break;  
  case 213:
    p=RhoPlus;
    break;
  case -213:
    p=RhoMinus;
    break;
  case 223:
    p=omega;
    break;
  case 331:
    p=EtaPrime;
    break;
  case 333:
    p=phiMeson;
    break;
  default:
    p=Unknown;
    break;
  }
  return((int)p);
}



/********************
 * PrintUsage
 *******************/
int PrintUsage(char *processName)
{
  
  fprintf(stderr,"%s usage: [switches]   \n",processName);
  fprintf(stderr,"\t-i<name> The input mcfast evt file.\n");
  fprintf(stderr,"\t-o<name> The output ascii file(default is stdout).\n");
  fprintf(stderr,"\t-N<#> The number mcfast events.\n");
  fprintf(stderr,"\t-n<#> The number particles per event.\n");
  fprintf(stderr,"\t-g<beamE> Write gamp output w/ kludged photon beam.\n");
  fprintf(stderr,"\t-h Print this help message\n\n");
}


