 /********************************************************
 * 
 *  Usage: genr8 <options> < input.gen
 *         Use "genr8 -h" for help with options.
 ********************************************************
 *                      * Generate t-channel        
 *  genr8.c           * monte carlo events.        
 *                      * 
 ********************************************************
 * 
 * created by:  Paul M Eugenio
 *              Carnegie Mellon University
 *              25-Mar-98
 ******************************************************** */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <malloc.h>
#include <sys/signal.h>
#include <sys/stat.h>
#include <time.h>

#include <disData.h>
#include <dataIO.h>

#include <genkin.h>
#include <particleType.h>
#include <ntypes.h>

#include <esr.h>

#define TRUE 1
#define FALSE 0
#define CONV (180.0/M_PI)
#define BUFSIZE 100000
#define T_CUT 10.0 
#define RECOIL 0.938

#define RESTFRAME -1
#define PARENTFRAME +1

#define PRODUCTION_PARTICLE 1




/***********************/
/* STRUCTURES          */
/***********************/



struct particleMC_t{
  int flag;
  int nchildren;
  int charge;
  double mass;
  double bookmass;
  double width;
  Particle_t particleID;
  vector4_t p;
  struct particleMC_t *parent, *child[2];
} ;


/******************************************************************
 * GLOBAL VARIABLES 
 *     NOTE: Please start the global name 
 *           with one capital letter!!!!!!
 *           Use all lower case for local 
 *           variable names.  Thank You 8^)
 **********************************************************************/
 
FILE *NewItape;
int Debug = 0;
int Nprinted =0;
int PrintProduction=0;
int PrintRecoil=0;
double MassHighBW;
int UseName=0;
int FIRST_EVENT=1;
int PrintFlag=10;
static itape_header_t *Buffer=NULL;
int WriteEsr=0;
int WriteAscii=0;
int runNo=9000;
/***********************/
/* Declarations         */
/***********************/
/*
 *  These functions are coded after main().
 */


double rawthresh(struct particleMC_t *Isobar);
int decay(struct particleMC_t *Isobar);
int boost2lab(struct particleMC_t *Isobar);
int boostFamily(vector4_t *beta,struct particleMC_t *Isobar);
int boost(vector4_t *beta,vector4_t *vec);
int printParticle(struct particleMC_t *Isobar);
int printParticle(struct particleMC_t *Isobar);
vector4_t polarMake4v(double p, double theta, double phi, double mass);
double randm(double low, double high);
int printProduction(FILE *fp,struct particleMC_t *Isobar);
int printFinal(FILE *fp,struct particleMC_t *Isobar);
int printp2ascii(FILE *fp,struct particleMC_t *Isobar);
int setMass(struct particleMC_t *Isobar);
int initMass(struct particleMC_t *Isobar);
int fillEsr(struct particleMC_t *Isobar,int *npart,esr_particle_t *esrp);
int WriteItape(struct particleMC_t *CM,vector4_t *beam);

/*
 ***********************
 *                     *
 *  PrintUsage()       *
 *                     *
 ***********************
 */

int PrintUsage(char *processName)
{
 
  fprintf(stderr,"%s usage: [-o<name>]   < infile \n",processName);
  fprintf(stderr,"\t-d debug flag\n");
  fprintf(stderr,"\t-n Use a particle name and not its ID number (ascii only) \n");
  fprintf(stderr,"\t-M<max> Process first max events\n");
  fprintf(stderr,"\t-l<lfevents> Determine the lorentz factor with this many number of events (default is 10000)\n");
  fprintf(stderr,"\t-r<runNo> default runNo is 9000. \n");
  /*fprintf(stderr,"\t-o<name> The output file \n");*/
  fprintf(stderr,"\t-P save flag= 11 & 01 events(default saves 11 & 10 events) \n");
  /* fprintf(stderr,"\t-R Save recoiling baryon information. \n"); */
  fprintf(stderr,"\t-I<filename> Save in itape format. \n");
  fprintf(stderr,"\t-A<filename> Save in ascii format. \n");
  fprintf(stderr,"\t-h Print this help message\n\n");

}


/*
 ***********************
 *                     *
 *  Main()             *
 *                     *
 ***********************
 */

main(int argc,char **argv)
{
  char *outputFile = "genr8.out";
  char *argptr,*token,line[2056];
  int i,npart=0,ngenerated=0,naccepted=0;
  int nv4,max=10,part,chld1=-1,chld2=-1,prnt=-1,lfevents=10000;
  FILE *fout=stdout;
  struct particleMC_t particle[20],beam,target,recoil,CM;
  struct particleMC_t *X,*Y;
  vector4_t beta,v4[2],initBeam4;
  double t,expt_max,expt,expt_min,sqrt_s,t_min=0;
  double CMenergy, t_max,slope=5.0;
  double X_momentum, X_threshold, X_energy,xmass,ymass;
  double costheta,theta,phi,lf,lfmax;
  time_t now;
  int isacomment=TRUE,haveChildren=TRUE;

  Y= &(particle[0]);
  Y->parent = &CM;
  X= &(particle[1]);
  X->parent = &CM;
  /* recoil.parent = &CM; */
  CM.child[0]= X;
  CM.child[1]= Y;
  /* CM.child[1]= &recoil; */

  NewItape=stdout;

  if (argc == 1){
    PrintUsage(argv[0]);
    exit (0);
  }
  else {
    
    /*
     *  Read command line options.
     */

    for (i=1; i<argc; i++) {
      argptr = argv[i];
      if ((*argptr == '-') && (strlen(argptr) > 1)) {
	argptr++;
	switch (*argptr) {
	case 'd':
	  Debug =1;
	  break;
	case 'h':
	  PrintUsage(argv[0]);
	  exit(0);
	  break;
	case 'n':
	  UseName =1;
	  break;
	case 'A':
	  WriteAscii=1;
	  fout = fopen(++argptr,"w");
	  fprintf(stderr,"Opening file %s for output. \n",argptr);
	  break;
	case 'R':
	  fprintf(stderr,"Printing recoil information.\n");
	  PrintRecoil=1;
	  break;
	case 'P':
	  fprintf(stderr,"Printing eta and pizeros and not gammas.\n");
	  PrintProduction=1;
	  break;
	case 'l':
	  lfevents = atoi(++argptr);
	  fprintf(stderr,"Using %d events to determine the lorentz factor: %d\n",lfevents);
	  break;
	case 'M':
	  max = atoi(++argptr);
	  fprintf(stderr,"Maximum number of events: %d\n",max);
	  break;
	case 'r':
	  runNo = atoi(++argptr);
	  fprintf(stderr,"Using runNo: %d\n",runNo);
	  break;
	case 'I':
	  WriteEsr=1;
	  fprintf(stderr,"Using itapes\n");
	  NewItape=fopen(++argptr,"w");
	  break;
	default:
	  fprintf(stderr,"Unrecognized argument -%s\n\n",argptr);
	  PrintUsage(argv[0]);	  
	  exit(-1);
	  break;
	  
	}
      }
    }
  }

  /*
   *  Seed the random number generator.
   */
  now=time(NULL);
  srand48(now); 

  /*
   * Now read the input.gen file 
   * from the stdin.
   *
   *
   * Any line starting with a "%"
   * is a comment line and is ignored.
   *
   *
   */


/* Fill particle information */

  isacomment=TRUE;
  while(isacomment==TRUE){
    fgets(line,sizeof(line),stdin);
    token=strtok(line," ");
    if(!(*token == '%'))
      isacomment=FALSE;
  }  /* get beam information */
  beam.p.space.x = atof(token);
  token=strtok(NULL," ");
  beam.p.space.y = atof(token);
  token=strtok(NULL," ");
  beam.p.space.z = atof(token);
  token=strtok(NULL," ");
  beam.mass =  atof(token);
  fprintf(stderr,"Reading: \tbeamp.x \tbeamp.y \tbeamp.z \tbeamMass\n");
  fprintf(stderr,"Found: \t\t%lf \t%lf \t%lf \t%lf \n",
	  beam.p.space.x, beam.p.space.y, beam.p.space.z, beam.mass);  


  isacomment=TRUE;
  while(isacomment==TRUE){
    fgets(line,sizeof(line),stdin);
    token=strtok(line," ");
    if(!(*token == '%'))
      isacomment=FALSE;
  }  /* get target information */
  target.p.space.x = atof(token);
  token=strtok(NULL," ");
  target.p.space.y = atof(token);
  token=strtok(NULL," ");
  target.p.space.z = atof(token);
  token=strtok(NULL," ");
  target.mass =  atof(token);
  fprintf(stderr,"Reading: \ttargetp.x \ttargetp.y \ttargetp.z \ttargetMass\n");
  fprintf(stderr,"Found: \t\t%lf \t%lf \t%lf \t%lf \n",
	  target.p.space.x, target.p.space.y, target.p.space.z, target.mass);  


  isacomment=TRUE;
  while(isacomment==TRUE){
    fgets(line,sizeof(line),stdin);
    token=strtok(line," ");
    if(!(*token == '%'))
      isacomment=FALSE;
  }  
  /* get the t-channel slope */
  slope=atof(token);
  fprintf(stderr,"Reading: t-channelSlope\n");
  fprintf(stderr,"Found: \t%lf \n",slope);

  isacomment=TRUE;
  while(isacomment==TRUE){
    fgets(line,sizeof(line),stdin);
    token=strtok(line," ");
    if(!(*token == '%'))
      isacomment=FALSE;
  }  /* get the number of particles to read in below */
  npart = atoi(token);
  fprintf(stderr,"Reading: number of particles need to describe the decay\n");
  fprintf(stderr,"Found: \t%d \n",npart);  

  /*
   * read all particles needed
   * to decsribe an isobar decay
   * of the resonance (X)
   */

  fprintf(stderr,"Reading: \tpart# \tchld1# \tchld2# \tprnt# \tId \tnchld \tmass \t\twidth \t\tchrg \tflag  \n"); 

  for(i=0;i<npart;i++){ /* read particle information*/
    haveChildren=TRUE;
    isacomment=TRUE;
    while(isacomment==TRUE){
      fgets(line,sizeof(line),stdin);
      token=strtok(line," ");
      if(!(*token == '%'))
	isacomment=FALSE;
    } 
      if(!(*token == '*')) /* "*" means it's unknown at this time */
	part = atoi(token);
      token=strtok(NULL," ");
      if(!(*token == '*')){  /* set pointer to child */ 
	chld1 = atoi(token);
	particle[part].child[0] = &(particle[chld1]);
      }
      else
	haveChildren = FALSE;
      token=strtok(NULL," ");     
      if(!(*token == '*')) {/* set pointer to child */ 
	chld2 = atoi(token);
	particle[part].child[1] = &(particle[chld2]);
      }
      else
	haveChildren = FALSE;
      token=strtok(NULL," ");
      if(!(*token == '*')) {/* set pointer to parent */
	prnt = atoi(token);
	particle[part].parent = &(particle[prnt]);
      }
      token=strtok(NULL," ");
      if(!(*token == '*'))
	/* particle[part].particleID = part; */
	particle[part].particleID = atoi(token);/* see particleType.h */
      token=strtok(NULL," ");
      if(!(*token == '*')){
	particle[part].nchildren = atoi(token);
	if(haveChildren==FALSE && particle[part].nchildren>0 ){
	  fprintf(stderr,
		  "If a particle has children then it must point to them!\n");
	  exit(-1);
	}
      }
      token=strtok(NULL," ");
      if(!(*token == '*'))
	particle[part].bookmass = atof(token);
      else{ /* get a list of the particle  that need a mass generated  */
	fprintf(stderr,"Every Particle needs a mass\n");
	exit(-1);
      }
      token=strtok(NULL," ");
      if(!(*token == '*')){
	particle[part].width = atof(token);
	
      }
      else {/* for a fixed mass use a zero width  */
	fprintf(stderr,"Every Particle needs a width\n");
	exit(-1);
      }
      token=strtok(NULL," ");
      if(!(*token == '*'))
	particle[part].charge = atoi(token);
      token=strtok(NULL," ");
      if(!(*token == '*'))
	particle[part].flag = atoi(token);
      /* flag  00 = isobar or resonace
       * flag  01 = production particle that decays i.e. eta, pizero ..
       * flag  11 = production particle that does not decay i.e. piplus,...
       * flag  10 = final state particle not in production i.e. gamma
       */
      fprintf(stderr,
	      "Found: \t\t%d \t%d \t%d \t%d \t%d \t%d \t%lf \t%lf \t%d \t%d\n",
	      part,chld1,chld2,prnt,particle[part].particleID,particle[part].nchildren,particle[part].bookmass,particle[part].width,particle[part].charge,particle[part].flag);
    }
  isacomment=TRUE;
  while(isacomment==TRUE){
    fgets(line,sizeof(line),stdin);
    token=strtok(line," ");
    if(!(*token == '%'))
      isacomment=FALSE;
  } 
  if(!(*token == '!')){
    fprintf(stderr,"Failed to find EOI---- Check Input File\n");
    exit(-1);
  }

  checkFamily(X);
  fprintf(stderr,"Found EOI----  Input File appears Fine.\n");

    /* We are now done reading the input information */

  /*
   * The beam and target are in the lab frame.
   * Put them in the overall center of momentum (CM) frame 
   * and calculate |t| & recoil angles.
   */

    if(X->nchildren == 0) X->mass = X->bookmass;
     if(Y->nchildren == 0) Y->mass = Y->bookmass;
      
    target.p.t = energy(target.mass,&(target.p.space));        
    beam.p.t = energy(beam.mass,&(beam.p.space));
    initBeam4.t= beam.p.t;  initBeam4.space.x= beam.p.space.x;
    initBeam4.space.y= beam.p.space.y; initBeam4.space.z= beam.p.space.z;
    sqrt_s = sqrt( SQ(beam.mass) +SQ(target.mass) + 2.0*beam.p.t * target.p.t);
    /*MassHighBW = sqrt_s - recoil.mass; */
    MassHighBW = sqrt_s; /* see do loop below */

  v4[0]= beam.p;
  v4[1]= target.p;
  nv4=2;
  
  CM.mass = sqrt_s;
  CM.p = Sum4vec(v4,nv4);
  beta = get_beta(&(CM.p),RESTFRAME);
  boost(&beta,&(beam.p));
  boost(&beta,&(target.p));
  
  initMass(X);
  initMass(Y);
  CMenergy = beam.p.t + target.p.t;
 
  while(naccepted <max){

    
  
    /* 
     * Generate the resonance 
     * in the CM frame, and
     * fill the four vectors
     * for both X and the recoil.
     */
   
    do{ 
    if(!(X->width<0)) 
      do{/*use BreitWigner--phasespace distribution */
	initMass(X);
	setMass(X);

	/*
	 * set the children mass to the book mass or
	 * distribute the by a Breit-Wigner. If the isobar
	 * mass is unknown it's mass remains unknow at this time.
	 */
	
	for(i=0;i<X->nchildren;i++)
	  setChildrenMass(X->child[i]); 

      }while((X->mass > MassHighBW) ||  ( X->nchildren==0 ? FALSE :
	    (X->mass <  ( (X->child[0])->mass +  (X->child[1])->mass))) );
    else{/* there's an error.. */
	fprintf(stderr,"Cannot use a negative width!\n");
	exit(-1); 
      }
    /*
     * Now do loop it for Y
     */
    if(!(Y->width<0)) 
      do{/*use BreitWigner--phasespace distribution */
	initMass(Y);
	setMass(Y);

	/*
	 * set the children mass to the book mass or
	 * distribute the by a Breit-Wigner. If the isobar
	 * mass is unknown it's mass remains unknow at this time.
	 */
	
	for(i=0;i<Y->nchildren;i++)
	  setChildrenMass(Y->child[i]); 

      }while((Y->mass > MassHighBW) || ( Y->nchildren==0 ? FALSE :
	    (Y->mass <  ( (Y->child[0])->mass +  (Y->child[1])->mass)) ) );
    else{/* there's an error.. */
	fprintf(stderr,"Cannot use a negative width!\n");
	exit(-1); 
      }
    }while(sqrt_s < X->mass + Y->mass);
    /*
      xmass=rawthresh(X);
      ymass=rawthresh(Y);
      *
      * fprintf(stderr," xmass= %lf ymass= %lf  X->mass= %lf Y->mass= %lf\n",
      * xmass,ymass, X->mass , Y->mass);
      */
    xmass = X->mass;
    ymass = Y->mass;

    X_momentum = CMmomentum( CMenergy, X->mass, Y->mass);
    X_energy =  sqrt( (X->mass)*(X->mass) + X_momentum*X_momentum);

    if(Y->nchildren ==0){
      
      t_min = -( SQ( (SQ(beam.mass) -SQ(xmass) -SQ(target.mass) +SQ(ymass))/(2.0*sqrt_s))
		 -SQ(v3mag(&(beam.p.space)) - X_momentum ));
      t_max = -( SQ( (SQ(beam.mass) -SQ(xmass) -SQ(target.mass) +SQ(ymass))/(2.0*sqrt_s))
		 -SQ(v3mag(&(beam.p.space)) + X_momentum ));
      /* 
       *fprintf(stderr,
       "beam.mass= %lf xmass= %lf target.mass=%lf ymass= %lf sqrt_s= %lf beam.p= %lf X->p= %lf  X_momentum= %lf\n",
       beam.mass,xmass,target.mass,ymass,sqrt_s,
       v3mag(&(beam.p.space)),v3mag(&(X->p.space)), X_momentum);
       */

      /*fprintf(stderr,"t_min: %lf t_max: %lf\n", t_min,t_max);
       */ 
    } else{ /* it's some baryon pseudo t process */
       
    t_min=0.4;
    t_max=10.0;
    t_min = -( SQ( (SQ(beam.mass) -SQ(xmass) -SQ(target.mass) +SQ(ymass))/(2.0*sqrt_s))
	       -SQ(v3mag(&(beam.p.space)) - X_momentum ));
    t_max = -( SQ( (SQ(beam.mass) -SQ(xmass) -SQ(target.mass) +SQ(ymass))/(2.0*sqrt_s))
	       -SQ(v3mag(&(beam.p.space)) + X_momentum ));
    
    }     
    expt_max = exp(-slope * t_max);
    expt_min = exp(-slope * t_min);

    do{
     
      expt = randm(expt_max,expt_min);
     
      t= -log(expt)/slope;
      costheta = ( beam.p.t * X_energy -
		   0.5*(t + (beam.mass)*(beam.mass) + (X->mass)*(X->mass))
		   )/( v3mag(&(beam.p.space))*X_momentum ) ;
      
    }while(fabs(costheta)>1.0  );
    
   
    
    theta = acos(costheta);
    phi = randm(-1*M_PI,M_PI);
    
    X->p = polarMake4v(X_momentum,theta,phi,X->mass);
    Y->p=polarMake4v(X_momentum,(M_PI-theta),(M_PI+phi),Y->mass);

    /*
     *   Now decay X -> children -> grandchildren -> and so forth
     *
     * Note: all particles are generated in their parent's rest frame.
     */



    /*
    if(Debug) 
      fprintf(stderr,"before decay\n");
    if(Debug)   
      printFamily(X);
      */

    decay(X);
    decay(Y);
    if(Debug) {
      fprintf(stderr,"X after decay\n");
      printFamily(X);
      fprintf(stderr,"Y after decay\n");
       printFamily(Y);
    }
   

    /*
     * Compute Lorentz Factor (used for phasespace weighting)
     */
    lf=v3mag(&(X->p.space));
    lorentzFactor(&lf,X);
    lorentzFactor(&lf,Y);
    if(lfevents-->0){
      lfmax = lf >lfmax ? lf : lfmax; /* find the largest value */
      if( (lfevents % 10) == 0 )
	fprintf(stderr,"Calculating Lorentz Factor: %d \r",lfevents);
    }
    else{
     
      /*
       * Now generate the events weighted by phasespace
       * (the maximum Lorentz factor).
       *
       * Since each particle is in its parent's rest frame,
       * it must be boosted through each parent's -> parent's-> ...
       * rest frame to the lab frame.
       */
      /*
      fprintf(stderr,"expt_min: %lf \t expt_max: %lf\n",expt_min,expt_max);
      fprintf(stderr,"t_min: %lf \t t: %lf \t t_max: %lf\n",t_min,t,t_max);
      */
      ngenerated++;
      /* fprintf(stderr,"ngen: %d   \tlf:  %lf \tlfmax: %lf\n",ngenerated,lf,lfmax);*/
      if(lf > randm(0.0,lfmax) ){ /* phasespace distribution */
	
	naccepted++;
	boost2lab(X);
	boost2lab(Y);
	initBeam4;  /* use lab frame beam */
	/*
	 * We have a complete event.  Now save it!
	 */
	
	if(Debug) {
	  fprintf(stderr,"X after boost2lab\n");
	  printFamily(X);
	  fprintf(stderr,"Y after boost2lab\n");
	  printFamily(Y);
	}
	
	Nprinted =0;
	/* event header information 
	fprintf(fout,"RunNo %d EventNo %d\n",runNo,naccepted);*/
	fprintf(fout,"%d %d\n",runNo,naccepted);

	/*
	 * Print out the production
	 * or the final state particles.
	 */
	if(WriteEsr)
	  WriteItape(&CM,&initBeam4); 

	if(PrintProduction){
	  if(X->nchildren==0)
	    printp2ascii(fout,X);
	  else
	    printProduction(fout,X);
	  if(Y->nchildren==0)
	    printp2ascii(fout,Y);
	  else
	    printProduction(fout,Y);
	}
	else{
	  if(X->nchildren==0)
	    printp2ascii(fout,X);
	  else
	    printFinal(fout,X);
	  if(Y->nchildren==0  && Y->flag/10 == 1 )
	    printp2ascii(fout,Y);
	  else
	    printFinal(fout,Y);
	}
	

      }
      if(!(ngenerated % 100))
	  fprintf(stderr,"Events generated: %d Events accepted: %d \r",
		  ngenerated,naccepted);
      if(Debug) fprintf(stderr,"End of event\n");
    } /* end of else{ */
  }/* end of while */
  
  fprintf(stderr,
	  "Max Lorentz Factor:%lf Events generated:%d Events acepted:%d\n\n",
	  lfmax,ngenerated,naccepted);
  /*
   * Close the output file.
   */
  fflush(fout);
  fclose(fout);
  data_flush(fileno(NewItape));
  fclose(NewItape);

  
}/* end of main */


/********************
 *
 *  checkFamily()
 * Testing the input file.
 *
 * If I have children
 * then I should be
 * my children's
 * parent. 
 *******************/
int checkFamily(struct particleMC_t *Isobar)
{
  int i;
  
  for(i=0;i<Isobar->nchildren; i++){
    if(Isobar->nchildren != 2){
      fprintf(stderr,"Error in input file: Sorry, only 0 or 2 children are allowed.\n");
      exit(-1);
    }
    if(Isobar != Isobar->child[i]->parent){
      fprintf(stderr,"Error in input file: Parent to children mismatch\n");
      exit(-1);
    }
     checkFamily(Isobar->child[i]);
  }

}



/********************
 *
 * setChildrenMass()
 *
 * Sets all children masses
 * to bookmass or to a 
 * Breit-Wigner mass
 *********************/

int setChildrenMass(struct particleMC_t *Isobar)
{
  int i;

  initMass(Isobar);
  setMass(Isobar);
  for(i=0;i < Isobar->nchildren;i++){
    setChildrenMass(Isobar->child[i]);
  }
  if(Isobar->nchildren !=0) 
    if((Isobar->mass) < ((Isobar->child[0])->mass)+((Isobar->child[1])->mass))
      setChildrenMass(Isobar);
  
  
}




/********************
 *
 * setMass()
 *
 * Sets the particle mass
 * using a 
 * Breit-Wigner distribution.
 *********************/

int setMass(struct particleMC_t *Isobar)
{
  double n,height,thresH2,lowtail,hightail,hcut,lcut;

    if(Isobar->width > 0){
 
      lowtail = rawthresh(Isobar); 
      
      
      
      thresH2 = 
	Isobar == (Isobar->parent)->child[0] ? /* is the 1st child me? */
	rawthresh((Isobar->parent)->child[1]) : /* if true */
	rawthresh((Isobar->parent)->child[0]) ;/* if false */
      /*
      hightail = ((Isobar->parent)->mass);
      */
      hightail = (Isobar->parent->mass) - thresH2 ;
      
      
      /* cut off the tails */
      hcut= Isobar->mass + 4.0*Isobar->width ;
      lcut= Isobar->mass - 4.0*Isobar->width ;
      if(hightail> hcut)
	hightail=hcut;
      if(lowtail < lcut)
	lowtail= lcut;


      do{
	n=randm(0.0,0.9999);
	Isobar->mass = randm( lowtail , hightail);
	
	height= SQ((Isobar->bookmass)*(Isobar->width))/
	  ( SQ(SQ(Isobar->bookmass) - SQ(Isobar->mass)) + 
	    SQ((Isobar->bookmass) * (Isobar->width) ));
      }while(n > height );
      /*
	fprintf(stderr,"bookmass is %lf: low: %lf high: %lf bwmass: %lf\n",
	Isobar->bookmass, lowtail,hightail,Isobar->mass);    
      */
    }
}
/********************
 *
 * initMass()
 *
 * Sets the particle mass
 * to bookmass or UNKNOWN 
 * 
 *********************/

int initMass(struct particleMC_t *Isobar)
{
  int i;

  for(i=0;i<Isobar->nchildren;i++){
    if(Isobar->child[i]->width == 0.0)
      Isobar->child[i]->mass = Isobar->child[i]->bookmass;
    else
      Isobar->child[i]->mass = -1.0; /* UNKNOWN */
    initMass(Isobar->child[i]);
  } 
}
/********************
 *
 * rawthresh()
 *
 * Calculates the mass
 * threshold for the 
 * isobar.
 *******************/

double rawthresh(struct particleMC_t *Isobar)
{
  int i;
  double rmassThresh=0.0;

  if(Isobar->nchildren){
    for(i=0; i < Isobar->nchildren;i++){
      if (Isobar->child[i]->mass <0)/* it is not known now */
	rmassThresh += rawthresh(Isobar->child[i]);
      else /* it is known now */
	rmassThresh += Isobar->child[i]->mass;
    }
  }else{
    rmassThresh=Isobar->mass;
    if(Isobar->mass <0){ /* error */
      fprintf(stderr,"Error!: Isobar->mass <0 for Isobar with no children\nExit\n");
      exit(-1);
    }
  }
  return rmassThresh;
}


/**********************************
 *
 * decay(Isobar)
 *
 * Decay the isobar into its children 
 * and then repeat to decay each 
 * child isobar.
 *************************************/

int decay(struct particleMC_t *Isobar)
{
  int i,j,k;
  double breakup_p,theta,phi;

 
  
  if(Isobar->nchildren>0)
    {
      breakup_p = CMmomentum(Isobar->mass,
			     Isobar->child[0]->mass,
			     Isobar->child[1]->mass);
      theta = acos(randm(-0.9999, 0.9999));
      phi = randm(-1*M_PI,M_PI);
      Isobar->child[0]->p = 
	polarMake4v(breakup_p,theta,phi,Isobar->child[0]->mass);
      Isobar->child[1]->p = 
	polarMake4v(breakup_p,(M_PI - theta),(M_PI + phi),Isobar->child[1]->mass);
      for(i=0;i<Isobar->nchildren;i++)
	decay(Isobar->child[i]);
    }
}


/********************************
 *
 * lorentzFactor()
 *
 * Returns the multiplication of
 * all of the break-up momenta.
 *
 *********************************/

int lorentzFactor(double *lf,struct particleMC_t *Isobar)
{
  int i;

  if(!(Isobar->nchildren == 0 )){
    *lf  *= v3mag(&(Isobar->child[0]->p.space));
    for(i=0;i<Isobar->nchildren;i++)
      lorentzFactor(lf,Isobar->child[i]);
  }
 
}


/********************************
 *
 * boost2lab()
 *
 * o Starting w/ final state particles
 *   boost to parent's frame.
 *
 * o Then boost parent & children to
 *   parent's parent's frame and repeat
 *   to the lab frame.
 *
 *********************************/
int boost2lab(struct particleMC_t *Isobar)
{
  int i;
  vector4_t beta;
  
  for(i=0;i<Isobar->nchildren;i++)
     boost2lab(Isobar->child[i]);
  
  beta = get_beta(&(Isobar->parent->p),PARENTFRAME);/* see kinematics.c */
  boostFamily(&beta,Isobar);
}

/********************************
 *
 * boostFamily()
 *
 * Boost particle and all children,
 * children's children, ...
 *
 *********************************/
int boostFamily(vector4_t *beta,struct particleMC_t *Isobar)
{
  int j;
  boost(beta,&(Isobar->p));
  for(j=0;j<Isobar->nchildren;j++)
    boostFamily(beta, Isobar->child[j]);
}

/********************************
 *
 * boost()
 *
 * Boost a four vector.
 *
 *********************************/
int boost(vector4_t *beta,vector4_t *vec)
{
  vector4_t temp;
  
  temp = lorentz(beta,vec);/* see kinematics.c */
  vec->t = temp.t;
  vec->space.x = temp.space.x;
  vec->space.y = temp.space.y;
  vec->space.z = temp.space.z;
}

/********************************
 *
 * printProduction()
 *
 * Print out production particles.
 *******************************/
int printProduction(FILE *fp,struct particleMC_t *Isobar)
{
  int i;
  
  for(i=0;i<Isobar->nchildren;i++){
    if((Isobar->child[i]->flag%10 ) == 1)
      printp2ascii(fp,Isobar->child[i]);
    printProduction(fp,Isobar->child[i]);
  }
}

/********************************
 *
 * printFinal()
 *
 * Print out final state particles
 *******************************/
int printFinal(FILE *fp,struct particleMC_t *Isobar)
{
  int i;
 
  for(i=0;i<Isobar->nchildren;i++){
    if((Isobar->child[i]->flag/10 ) == 1)
      printp2ascii(fp,Isobar->child[i]);
    printFinal(fp,Isobar->child[i]);
  }
}

/********************************
 *
 * printp2ascii()
 *
 *******************************/
int printp2ascii(FILE *fp,struct particleMC_t *Isobar)
{
  Nprinted++;
  if(UseName)
    fprintf(fp,"%d %s %lf\n",Nprinted,ParticleType(Isobar->particleID),Isobar->mass);
  else
    fprintf(fp,"%d %d %lf\n",Nprinted,Isobar->particleID,Isobar->mass);
  
  fprintf(fp,"   %d %lf %lf %lf %lf\n",Isobar->charge,
	  Isobar->p.space.x,
	  Isobar->p.space.y,
	  Isobar->p.space.z,
	  Isobar->p.t);

}

/********************************
 *
 * printFamily()
 *
 *******************************/
int printFamily(struct particleMC_t *Isobar)
{
  int j;
  printParticle(Isobar);
  for(j=0;j<Isobar->nchildren;j++)
    printFamily(Isobar->child[j]); 
}

/*************************
 *
 * WriteItape((struct particleMC_t *CM,vector4_t *beam)
 *
 *************************/
int WriteItape(struct particleMC_t *CM, vector4_t *beam) 
{
  int group=802; /* esrGroup */
  int nparticles,i;
  esr_nparticle_t *esr;
  esr_particle_t esrp[20];
  static int icount;

  if(FIRST_EVENT==TRUE){
    Buffer = (itape_header_t *) malloc(BUFSIZE);
    
    data_newItape(Buffer);
    FIRST_EVENT=FALSE;
  }

  nparticles=0;
  fillEsr(CM->child[0],&nparticles,esrp);
  fillEsr(CM->child[1],&nparticles,esrp);
 

  Buffer->runNo= runNo;
  Buffer->spillNo =0;
  Buffer->eventNo =icount++;
  Buffer->trigger =1;
  esr= data_addGroup(Buffer,BUFSIZE,802,sizeof(esr_nparticle_t)+ 
		     nparticles*sizeof(esr_particle_t));
  esr->nparticles=nparticles;
  for(i=0;i<nparticles;i++)
    esr->p[i] = esrp[i]; 
  esr->beam = *beam;
  /*  esr->miss=CM->child[1]->p;
   */
  data_saveGroups(Buffer,1,&group);
  data_clean(Buffer);
  data_write(fileno(NewItape),Buffer);
}

int fillEsr(struct particleMC_t *Isobar,int *npart,esr_particle_t *esrp)
{
  int i;
  
  if(Isobar->flag == 11 || Isobar->flag==PrintFlag)
    {
      (esrp + *npart)->p = Isobar->p;
      (esrp + *npart)->charge = Isobar->charge;
      (esrp + *npart)->particleType = Isobar->particleID;
      (*npart)++;
      if(*npart>19){
	fprintf(stderr,"Increase the number of esr particles in WriteItape()\n");
	exit(-1);
      }
    }
  for(i=0;i<Isobar->nchildren;i++)
    fillEsr(Isobar->child[i],npart,esrp);
}



/********************************
 *
 * printParticle()
 *
 *******************************/
int printParticle(struct particleMC_t *Isobar)
{
  fprintf(stderr,"Particle ID %s with %d children\n",
	  ParticleType(Isobar->particleID),Isobar->nchildren);
  fprintf(stderr,"four momentum (E,p): %lf %lf %lf %lf\n\n",
	  Isobar->p.t,
	  Isobar->p.space.x,
	  Isobar->p.space.y,
	  Isobar->p.space.z);
}

/******************************************************
 *
 * polarMake4v()
 *
 * make a four vector given (p,theta,phi) and it's mass
 ********************************************************/
vector4_t polarMake4v(double p, double theta, double phi, double mass)
{
  vector4_t temp;

  temp.t = sqrt( SQ(mass) + SQ(p));
  temp.space.z = p*cos(theta);
  temp.space.x =  p*sin(theta)*cos(phi);
  temp.space.y =  p*sin(theta)*sin(phi);
  
  return temp;
}

/********************************
 *
 * randm(double low, double high)
 *
 *******************************/
double randm(double low, double high)
{
  /* Seed the random number generator using:
   *  int now = time(NULL);
   *  srand48(now);
   */
  return ((high - low) * drand48() + low);
}



/*
 ***********************
 *                     *
 *  END OF FILE        *
 *                     *
 ***********************
 */




