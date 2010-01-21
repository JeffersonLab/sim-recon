/************************************************
 * stdhep2hddm.c 
 * This program converts the StdHep format to a
 * generic hddm format.
 *
 *  See http://www-pat.fnal.gov/stdhep.html
 *      http://zeus.phys.uconn.edu/halld/datamodel/doc
 *
 * Richard Jones
 * University of Connecticut
 * June 1, 2001
 **********************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "hddm_s.h"

#include <stdhep.h>
#include <stdlun.h>
#include <stdcnt.h>

int runNo=-9000;

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


/********************************
 * Prototypes for StdHep functions.
 *********************************/

int StdHepXdrReadInit(char *fileName,int ntries, int istream);
int StdHepXdrRead(int *ilbl,int istream);
int StdHepXdrEnd(int istream);


int  gampID(int id)
{
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


int fill_mc_part(int i, s_Vertices_t* vs, int v)
{
  Particle_t ptype;
  s_Origin_t* or = vs->in[v].origin;
  s_Products_t* ps = vs->in[v].products;
  if (or == &hddm_nullTarget)
  {
    or = make_s_Origin();
    vs->in[v].origin = or;
    or->vx = hepevt_.vhep[i][0];
    or->vy = hepevt_.vhep[i][1];
    or->vz = hepevt_.vhep[i][2];
  }
  if (ps == &hddm_nullTarget)
  {
    ps = make_s_Products(30);
    vs->in[v].products = ps;
    ps->mult = 0;
  }
  ptype = gampID(hepevt_.idhep[i]);
  ps->in[ps->mult].type = ptype;
  ps->in[ps->mult].momentum = make_s_Momentum();
  ps->in[ps->mult].momentum->px = hepevt_.phep[i][0];
  ps->in[ps->mult].momentum->py = hepevt_.phep[i][1];
  ps->in[ps->mult].momentum->pz = hepevt_.phep[i][2];
  ps->in[ps->mult].momentum->E  = hepevt_.phep[i][3];
  ps->in[ps->mult].properties = make_s_Properties();
  ps->in[ps->mult].properties->mass = hepevt_.phep[i][4];
  ps->in[ps->mult].properties->charge = ParticleCharge(ptype);
  ps->mult++;

  if (hepevt_.jdahep[i][0] != 0)
  {
    int j;
    for (j = 0; j < hepevt_.nhep; j++)
    {
      int iv = vs->mult;
      if (hepevt_.jmohep[j][0] == i)
      {
        fill_mc_part(j,vs,iv);
      }
    }
    vs->mult++;
  }
}


int fill_mc_parts(s_HDDM_t* mc_evt)
{
  int i;
  s_PhysicsEvents_t* pes = make_s_PhysicsEvents(1);
  s_Reactions_t* rs = make_s_Reactions(1);
  s_Vertices_t* vs = make_s_Vertices(10);
  mc_evt->physicsEvents = pes;
  pes->in[0].reactions = rs;
  rs->in[0].vertices = vs;
  pes->mult = 1;
  rs->mult = 1;
  vs->mult = 1;
  for (i = 0; i < hepevt_.nhep; i++)
  {
    if (hepevt_.jmohep[i][0] == 0)
    {
      fill_mc_part(i,vs,0);
    }
  }
  pes->in[0].runNo = runNo;
  pes->in[0].eventNo = hepevt_.nevhep;
  return 1;
}


int PrintUsage(char *processName)
{
  fprintf(stderr,"%s usage: [switches]   \n",processName);
  fprintf(stderr,"\t-i<name> input stdhep evt file (no default)\n");
  fprintf(stderr,"\t-o<name> output hddm file (default is stdhep.hddm)\n");
  fprintf(stderr,"\t-N<#> number stdhep events to process (default is 0)\n");
  fprintf(stderr,"\t-r<#> run number saved in events (default is -9000)\n");
  fprintf(stderr,"\t-h Print this help message\n\n");
}


int main(int argc,char **argv)
{
  char *argptr;
  int i, ntries=0, ret, written=0;
  char *evtfile = "default.evt";
  int istream=0, ilbl;
  s_HDDM_t *mc_evt;
  s_iostream_t *outputfp;
  char hddmfile[400];
  strcpy(hddmfile,"stdhep.hddm");

  if (argc == 1)
  {
    PrintUsage(argv[0]);
    exit (0);
  }
  else
  {
    for (i = 1; i < argc; i++) {
      argptr = argv[i];
      if ((*argptr == '-') && (strlen(argptr) > 1)) {
        argptr++;
        switch (*argptr) {
        case 'N':
          ntries=atoi(++argptr);
          break;
	case 'o':
          sprintf(hddmfile,"%s.hddm",++argptr);
          break;
	case 'i':
	  evtfile= ++argptr;
          break;
	case 'r':
	  runNo= atoi(++argptr);
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

/*
 * Open and init a stdhep file.
 */
    
    if (ret=StdHepXdrReadInit(evtfile,ntries,istream))
    {
      /* some error must have occured */
      fprintf(stderr,"err:StdHepXdrReadInit ret=%d\n",ret);
      exit(-1);
    }

    if (!(outputfp = init_s_HDDM(hddmfile)))
    {
      fprintf(stderr,"Fail to open output file!\n");
      exit(-1);
    }

    for (i=0;i<ntries;i++)
    {
      if (ret = StdHepXdrRead(&ilbl,istream))
      {/* some error must have occured */
	fprintf(stderr,"err:StdHepXdrRead ret=%d\n",ret);
	exit(-1); 
      }
      mc_evt = make_s_HDDM();
      fill_mc_parts(mc_evt);
      flush_s_HDDM(mc_evt,outputfp);
      if (!(++written %100))
      {
	fprintf(stderr,"stdhep events read: %d\r",written);  
      }
    }

    fprintf(stderr,"\nTotal stdhep events read: %d\n",written);  
  
    StdHepXdrEnd(istream);
    return 0;
  }
  return 9;
}
