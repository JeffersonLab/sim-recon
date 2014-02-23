/************************************************
 * ascii2stdhep.c 
 * This program converts the ascii output of genr8
 * to StdHep format. See http://www-pat.fnal.gov/stdhep.html
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

#define TRUE 1 
#define FALSE 0

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
int getFromGampFile(FILE *fp,int *nparts,mc_part_t *P);
int fill_hepevt(int nparts, mc_part_t *parts);
int  pdgID(Particle_t p);


/****************
 * main()
 **************/
 /* 
  * MAX_PARTS sets the max numbers of particles in an event
  */
#define MAX_PARTS 20
int main(int argc,char **argv)
{
  char *argptr;
  int i,ntries=0,nskip=0,ret,gotANevent,written=0;
  int ilbl=1,istream=0,nparts=0;/*I guessed at istream */ 
  int read_gamp_file=FALSE;
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
	case 'g':
	  read_gamp_file=TRUE;
	  break;
	case 'S':
          nskip=atoi(++argptr);
	  fprintf(stderr,"Skipping the first %d events\n",nskip);
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
      
    /* Read asci event. */
    if(read_gamp_file)
      gotANevent=getFromGampFile(inputfp,&nparts,parts);
    else
      gotANevent=getFromAscii(inputfp,nparts,parts);
    while(gotANevent>0){  /* I have an event! */
      if(!(nskip-- > 0)){
	if( gotANevent==-2){
	  fprintf(stderr,"Broken input file.\n");
	  exit(-1);
	}
	/* fill the hepevt structure */
	fill_hepevt(nparts, parts);
	
	/* write to stdhep file */
	ret=StdHepXdrWrite(ilbl,istream);
	if(!(++written %100))
	  fprintf(stderr,"McFast events Written: %d\r",written);
      } else{
	if(!(nskip %100))
	  fprintf(stderr,"                          \rSkipping events: %d\r",nskip);
      }
      if(read_gamp_file)
	gotANevent=getFromGampFile(inputfp,&nparts,parts);
      else
	gotANevent=getFromAscii(inputfp,nparts,parts);
      
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

int getFromGampFile(FILE *fp,int *nparts, mc_part_t *P)
{
  char line[2056];
  char *token;
  static int nread=0;
  int i;
  Particle_t genr8id;
  double mass;

  /* get the event header info */
  
  if(fgets(line,sizeof(line),fp)!=NULL){
    if(Debug>1)
      fprintf(stderr,"header: %s\n",line);
    nread++;
    token=strtok(line," ");/* nparts + beam */
    *nparts = atoi(token);

    /* skip the beam */
    fgets(line,sizeof(line),fp);
    (*nparts)--;
    if(Debug>1)
      fprintf(stderr," beam: %s\n",line);
    /* get the particle information */
    for(i=0;i<*nparts;i++){ 
      if(fgets(line,sizeof(line),fp)!=NULL){
	if(Debug>1)
	  fprintf(stderr,"part[%d]: %s\n",i,line);

        token=strtok(line," ");
        genr8id = atoi(token); 
	P[i].pid= pdgID(genr8id);
	P[i].statusCode=1;

        token=strtok(NULL," "); /* the charge */
        token=strtok(NULL," ");
        P[i].p.space.x = atof(token);
        token=strtok(NULL," ");
        P[i].p.space.y = atof(token);
        token=strtok(NULL," ");
        P[i].p.space.z = atof(token);
        token=strtok(NULL," ");
        P[i].p.t = atof(token);
        mass = sqrt( P[i].p.t *  P[i].p.t - (P[i].p.space.x * P[i].p.space.x +
					     P[i].p.space.y * P[i].p.space.y +
					     P[i].p.space.z * P[i].p.space.z));
					     
	P[i].mass = mass;
	if(Debug>1) 
          fprintf(stderr,"\tThe four momentum is %lf %lf %lf %lf \n",
                  P[i].p.space.x,P[i].p.space.y,P[i].p.space.z,P[i].p.t);
        if(Debug>1) 
          fprintf(stderr,"\tThe  mass is %lf\n",P[i].mass);
	
      } else return -2;

      
    }/* end of nparts */
    if(Debug>1) 
       fprintf(stderr,"It's a good event.\n");
    return 1;/* it looks like a good event */
  } else return -1;/* there are no more events */  
}


/*********************
 *
 *   getEvent
 *
 *******************/

int getFromAscii(FILE *fp,int nparts, mc_part_t *P)
{
  char line[2056];
  char *token;
  static int nread=0;
  int i;
  Particle_t genr8id;
  double mass;

  /* get the event header info */
  
  if(fgets(line,sizeof(line),fp)!=NULL){
    nread++;
    token=strtok(line," ");
    /* runNo = atoi(token); */
    
    token=strtok(NULL," ");
    /* conf = atof(token); */
     
    
    /* get the particle information */
for(i=0;i<nparts;i++){ /* four particles */
      if(fgets(line,sizeof(line),fp)!=NULL){

        token=strtok(line," ");
	/* particleNo[i] = atoi(token); */

        token=strtok(NULL," ");
        genr8id = atoi(token); 
	P[i].pid= pdgID(genr8id);
	P[i].statusCode=1;
        token=strtok(NULL," ");
        mass = atof(token);
	P[i].mass = mass;


      } else return -2;

      /* get next line */
      if(fgets(line,sizeof(line),fp)!=NULL){

        token=strtok(line," ");
        /* the charge 
        if(Debug>2) 
          fprintf(stderr,"Particle %d with charge = %d \n",i,P[i].iso3); */
        token=strtok(NULL," ");
        P[i].p.space.x = atof(token);
        token=strtok(NULL," ");
        P[i].p.space.y = atof(token);
        token=strtok(NULL," ");
        P[i].p.space.z = atof(token);
        token=strtok(NULL," ");
        P[i].p.t = atof(token);
        
        if(Debug>2) 
          fprintf(stderr,"\tThe four momentum is %lf %lf %lf %lf \n",
                  P[i].p.space.x,P[i].p.space.y,P[i].p.space.z,P[i].p.t);
        if(Debug>2) 
          fprintf(stderr,"\tThe  mass is %lf\n",P[i].mass);

      } else return -2;/* bad input file */
    } /* end of for(i=0 ... four particles */
    return 1;/* it looks like a good event */
  } else return -1;/* there are no more events */
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
   fprintf(stderr,"\t-i<name> The input ascii file(default is stdin).\n");
  fprintf(stderr,"\t-o<name> The output  file.\n");
  fprintf(stderr,"\t-N<#> The number output events. (N <= NasciiEvents)\n");
  fprintf(stderr,"\t-S<#> Skip the first # number of input events. \n");
  fprintf(stderr,"\t-n<#> The number particles per event.\n");
  fprintf(stderr,"\t-g Read gamp event format\n");
  fprintf(stderr,"\t-h Print this help message\n\n");
}


