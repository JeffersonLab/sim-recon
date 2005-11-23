/************************************************
 * itape2stdhep.c 
 * This program converts the esr itape output of genr8
 * to StdHep format. See http://www-pat.fnal.gov/stdhep.html
 * Paul Eugenio
 * Carnegie Mellon University
 * 24 Sept 98
 **********************************************/

#include<stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/*
 * needed for itape */
#include <sys/types.h>
#include <malloc.h>
#include <sys/signal.h>
#include <sys/stat.h>
#include <time.h>

#include <disData.h>
#include <dataIO.h>
#include <ntypes.h>
#include<esr.h>

#include<particleType.h>
#include<stdhep.h>
#include<stdlun.h>
#include<stdcnt.h>

/*
 * #include <stdhep.h> See it below.
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


/******************
 *  Local Structure
 ******************/
typedef struct { double x,y,z; } Dvector3_t;
typedef struct { double t; Dvector3_t space; } Dvector4_t;


typedef struct {
  int statusCode; /* 1 is a final state particle */
  int pid; /* this use the PDG MC numbering scheme */
  double mass;
  Dvector4_t p;
}mc_part_t;

/***************
 *  GLOBALS
 ***************/
int Debug=0;


/********************************
 * Prototypes for StdHep functions.
 *********************************/

int StdHepXdrWiteInit(char *fileName,char *title, int ntries, int istream);
int StdHepXdrWrite(int ilbl,int istream);
int StdHepXdrEnd(int istream);


/********************************
 * Local Prototypes 
 *********************************/
int PrintUsage(char *processName);
int getFromAscii(FILE *fp,int nparts,mc_part_t *P);
int fill_hepevt(int nparts, mc_part_t *parts);
int  pdgID(Particle_t p);


/****************
 * main()
 **************/
 /* 
  * MAX_PARTS sets the max numbers of particles in an event
  */
#define MAX_PARTS 20
#define BUFSIZE 200000
int main(int argc,char **argv)
{
  char *argptr;
  int i,ntries=0,ret,gotANevent,written=0;
  int ilbl=1,istream=0,nparts=0;/*I guessed at istream */ 
  char *outputfile ="default.evt";
  char *title ="The HallD MCfast";
  char inputfile[40];
  FILE *inputfp=stdin;
  mc_part_t parts[MAX_PARTS];



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
        case 'N':
          ntries=atoi(++argptr);
          break;
	case 'o':
          outputfile = ++argptr;
          break;
	case 'i':
          sprintf(inputfile,"%s",++argptr);
          if(!(inputfp = fopen(inputfile,"r"))){
            fprintf(stderr,"Fail to open input file!\n");
            exit(-1);
          }
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
    
    if(ret=StdHepXdrWriteInit(outputfile,title,ntries, istream)){
      /* some error must have occured */
      fprintf(stderr,"err:StdHepXdrWiteInit ret=%d\n");
      exit(-1);
    }
    
    /*
     * Fill the HEPEVENT structure.
     */
      
    /* raed itape */
    gotANevent=getFromItape(inputfp,nparts,parts);
    
    while(gotANevent>0){  /* I have an event! */
      if( gotANevent==-2){
	fprintf(stderr,"Broken input file.\n");
	exit(-1);
      }
      /* fill the hepevt structure */
      fill_hepevt(nparts, parts);
      
      /* write to stdhep file */
      ret=StdHepXdrWrite(ilbl,istream);
      gotANevent=getFromItape(inputfp,nparts,parts);
      if(!(++written %100))
	fprintf(stderr,"McFast events Written: %d\r",written);
      if(written == ntries)
	gotANevent=0; /* Stop reading events */
    }
    if(written<ntries)
      fprintf(stderr,"Warning!! wanted %d events found only %d\n",
	      ntries,written);   
    /*
     * Close the StdHep file and exit.
     */

    StdHepXdrEnd(istream); 
    exit(0);

   } /* end of else /* good start */
} /* end of main */

/*********************
 *
 *   getEvent
 *
 *******************/

int getFromItape(FILE *fp,int nparts, mc_part_t *P)
{
  static int nread=0;
  static itape_header_t *buffer=NULL;
  esr_nparticle_t *esr;
  int i,ret;
  Particle_t genr8id;
  double mass;

  if(!buffer)
    buffer=(itape_header_t *)malloc(BUFSIZE);

  /* get the event header info */
  ret=data_read(fileno(fp),buffer,BUFSIZE);
  if(ret<0){
    fprintf(stderr,"err: read itape\n");
    exit(-1);
  }
  if(ret==0 || ret==2){ /* end of data */
    return 0;
  }
  /*
   * we have an event
   */
  esr=data_getGroup(buffer,GROUP_ESR_NPARTICLE);
  nread++;

    /* get the particle information */
  for(i=0;i<nparts;i++){ /* four particles */
    genr8id = esr->p[i].particleType;
    P[i].pid= pdgID(genr8id);
    P[i].statusCode=1;
    
    P[i].p.space.x = esr->p[i].p.space.x;
    P[i].p.space.y = esr->p[i].p.space.y; 
    P[i].p.space.z = esr->p[i].p.space.z; 
    P[i].p.t = esr->p[i].p.t;
 
    mass = sqrt((P[i].p.t)*(P[i].p.t) -( P[i].p.space.x *  P[i].p.space.x +
				    P[i].p.space.y *  P[i].p.space.y +
				    P[i].p.space.z *  P[i].p.space.z  ));
    P[i].mass = mass;
        
    
    if(Debug>2) 
      fprintf(stderr,"\tThe four momentum is %lf %lf %lf %lf \n",
	      P[i].p.space.x,P[i].p.space.y,P[i].p.space.z,P[i].p.t);
    if(Debug>2) 
      fprintf(stderr,"\tThe  mass is %lf\n",P[i].mass);


  }
   
  return 1;/* it looks like a good event */
  
}

/********************
 * fill_hepevt
 *******************/
int fill_hepevt(int nparts, mc_part_t *part){
  int i;
  static int nevent=0;

  nevent++;
  /* hepevt header info */
  hepevt_.nevhep=nevent;
  hepevt_.nhep=nparts;
  if(Debug>1)
    fprintf(stderr,"EventNo: %d\n",nevent);
  /* now loop over the particle in the event */
  
    
    for(i=0;i<nparts;i++){
      hepevt_.isthep[i]=part[i].statusCode;
       hepevt_.idhep[i]=part[i].pid;
       hepevt_.phep[i][0]=part[i].p.space.x;
       hepevt_.phep[i][1]=part[i].p.space.y;
       hepevt_.phep[i][2]=part[i].p.space.z;
       hepevt_.phep[i][3]=part[i].p.t;
       hepevt_.phep[i][4]=part[i].mass;
       hepevt_.vhep[i][0]=0.0; /* v.x */
       hepevt_.vhep[i][1]=0.0; /* v.y */
       hepevt_.vhep[i][2]=0.0; /* v.z */
       /* file the left overs */
       /* hepevt_.jmohep[i]  */
      if(Debug>1){
	fprintf(stderr,"StatusCode: %d  Pid: %d  Mass: %lf\n",
		part[i].statusCode, part[i].pid, part[i].mass);
	fprintf(stderr,"(p.x, p.y, p.x, p.t): (%lf, %lf, %lf, %lf)\n",
		part[i].p.space.x,
		part[i].p.space.y,
		part[i].p.space.z,
		part[i].p.t); 
      }
  }
  return 1;
}

/********************
 * pdgID(int genr8id)
 * The convention is that
 * followed in StdHep 4.02
 *******************/
int  pdgID(Particle_t p){
  int id=0;
  switch (p) {
  case Unknown:
    
    break;
  case Gamma:
    id=22;
    break;
  case Positron:
    id=-11;
    break;
  case Electron:
    id=11;
    break;
  case Neutrino:
    id=12;
    break;
  case MuonPlus:
    id=-13;
    break;
  case MuonMinus:
    id=13;
    break;
  case Pi0:
    id=111;
    break;
  case PiPlus:
    id=211;
    break;
  case PiMinus:
    id=-211;
    break;
  case KLong:
    id=130;
    break;
  case KPlus:
    id=321;
    break;
  case KMinus:
    id=-321;
    break;
  case Neutron:
    id=2112;
    break;
  case Proton:
    id=2212;
    break;
  case AntiProton:
    id=-2212;
    break;
  case KShort:
    id=310;
    break;
  case Eta:
    id=221;
    break;
  case Lambda:
    id=3122;
    break;
  case SigmaPlus:
    id=3222;
    break;
  case Sigma0:
    id=3212;
    break;
  case SigmaMinus:
    id=3112;
    break;
  case Xi0:
    id=3322;
  case XiMinus:
    id=3312;
    break;
  case OmegaMinus:
    id=3334;
    break;
  case AntiNeutron:
    id=-2112;
    break;
  case AntiLambda:
    id=-3122;
    break;
  case AntiSigmaMinus:
     id=-3112;
    break;
  case AntiSigma0:
    id=-3212;
    break;
  case AntiSigmaPlus:
    id=-3222;
    break;
  case AntiXi0:
    id=-3322;
    break;
  case AntiXiPlus:
    id=-3312;
    break;
  case AntiOmegaPlus:
    id=-3334;
    break;
  case Rho0:
    id=113;
    break;  
  case RhoPlus:
    id=213;
    break;
  case RhoMinus:
    id=-213;
    break;
  case omega:
    id=223;
    break;
  case EtaPrime:
    id=331;
    break;
  case phiMeson:
    id=333;
    break;
  default:
    id=0;
    break;
  }
  return(id);
}


/********************
 * PrintUsage
 *******************/
int PrintUsage(char *processName)
{
  
  fprintf(stderr,"%s usage: [switches]   \n",processName);
   fprintf(stderr,"\t-i<name> The input itape file(default is stdin).\n");
  fprintf(stderr,"\t-o<name> The output  file.\n");
  fprintf(stderr,"\t-N<#> The number output events. (N <= NasciiEvents)\n");
  fprintf(stderr,"\t-n<#> The number particles per event.\n");
  fprintf(stderr,"\t-h Print this help message\n\n");
}


