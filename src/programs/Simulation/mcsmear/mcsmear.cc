// $Id: mcsmear.cc 19023 2015-07-14 20:23:27Z beattite $
//
// Created June 22, 2005  David Lawrence

#include <iostream>
#include <iomanip>
using namespace std;

#include <TF1.h>
#include <TFile.h>
#include <TH2.h>
#include <TH1.h>

#include <signal.h>
#include <time.h>

#include <DANA/DApplication.h>
#include "MyProcessor.h"
#include "JFactoryGenerator_ThreadCancelHandler.h"

#include "units.h"
#include "HDDM/hddm_s.hpp"

void Smear(hddm_s::HDDM *record);
void ParseCommandLineArguments(int narg, char* argv[]);
void Usage(void);

extern void SetSeeds(const char *vals);

char *INFILENAME = NULL;
char *OUTFILENAME = NULL;
int QUIT = 0;

using namespace jana;

// histogram
pthread_mutex_t root_mutex = PTHREAD_MUTEX_INITIALIZER;


//-----------
// main
//-----------
int main(int narg,char* argv[])
{
   ParseCommandLineArguments(narg, argv);

   // Create DApplication object
   DApplication dapp(narg, argv);
   dapp.AddFactoryGenerator(new JFactoryGenerator_ThreadCancelHandler());
   
   TFile *hfile = new TFile("smear.root","RECREATE","smearing histograms");

   MyProcessor myproc;   
   dapp.Run(&myproc);
   
   hfile->Write();
   hfile->Close();

   return 0;
}

//-----------
// ParseCommandLineArguments
//-----------
void ParseCommandLineArguments(int narg, char* argv[])
{

  for (int i=1; i<narg; i++) {
    char *ptr = argv[i];
    
    if (ptr[0] == '-') {
      switch(ptr[1]) {
      case 'h': Usage();                                   break;
      case 'o': OUTFILENAME = strdup(&ptr[2]);             break;
      case 'N': ADD_NOISE=true;                            break;
      case 's': SMEAR_HITS=false;                          break;
      case 'i': IGNORE_SEEDS=true;                         break;
      case 'r': SetSeeds(&ptr[2]);                         break;
      case 'd': DROP_TRUTH_HITS=true;                      break;
      }
    }
    else {
      INFILENAME = argv[i];
    }
  }

   if (!INFILENAME){
      cout << endl << "You must enter a filename!" << endl << endl;
      Usage();
   }
  
   
   // Generate output filename based on input filename
   if (OUTFILENAME == NULL) {
      char *ptr, *path_stripped, *pdup;
      path_stripped = ptr = pdup = strdup(INFILENAME);
      while((ptr = strstr(ptr, "/")))path_stripped = ++ptr;
      ptr = strstr(path_stripped, ".hddm");
      if(ptr)*ptr=0;
      char str[256];
      sprintf(str, "%s_%ssmeared.hddm", path_stripped, ADD_NOISE ? "n":"");
      OUTFILENAME = strdup(str);
      free(pdup);
   }
   
}


//-----------
// Usage
//-----------
void Usage(void)
{
   cout << endl << "Usage:" << endl;
   cout << "     mcsmear [options] file.hddm" << endl;
   cout << endl;
   cout << " Read the given, Geant-produced HDDM file as input and smear" << endl;
   cout << "the truth values for \"hit\" data before writing out to a" << endl;
   cout << "separate file. The truth values for the thrown particles are" << endl;
   cout << "not changed. Noise hits can also be added using the -n option (deprecated for the moment)." << endl;
   cout << "Note that all smearing is done using Gaussians." << endl;
   cout << endl;
   cout << "  options:" << endl;
   cout << "    -ofname  Write output to a file named \"fname\" (default auto-generate name)" << endl;
   cout << "    -s       Don't smear real hits (default is to smear)" << endl;
   cout << "    -i       Ignore random number seeds found in input HDDM file" << endl;
   cout << "    -r\"s1 s2 s3\" Set initial random number seeds" << endl;
//   cout << "    -u#      Sigma CDC anode drift time in ns (def:" << CDC_TDRIFT_SIGMA*1.0E9 << "ns)" << endl;
//   cout << "             (NOTE: this is only used if -y is also specified!)" << endl;
//   cout << "    -y       Do NOT apply drift distance dependence error to" << endl;
//   cout << "             CDC (default is to apply)" << endl;
//   cout << "    -Y       Apply constant sigma smearing for FDC drift time. "  << endl;
//   cout << "             Default is to use a drift-distance dependent parameterization."  << endl;
//   cout << "    -t#      CDC time window for background hits in ns (def:" << CDC_TIME_WINDOW*1.0E9 << "ns)" << endl;
//   cout << "    -U#      Sigma FDC anode drift time in ns (def:" << FDC_TDRIFT_SIGMA*1.0E9 << "ns)" << endl;
//   cout << "    -C#      Sigma FDC cathode strips in microns (def:" << FDC_TDRIFT_SIGMA << "ns)" << endl;
//   cout << "    -T#      FDC time window for background hits in ns (def:" << FDC_TIME_WINDOW*1.0E9 << "ns)" << endl;
//   cout << "    -e       hdgeant was run with LOSS=0 so scale the FDC cathode" << endl;
//   cout << "             pedestal noise (def:false)" << endl;
   cout << "    -d       Drop truth hits (default: keep truth hits)" << endl;
//   cout << "    -p#      FCAL photo-statistics smearing factor in GeV^3/2 (def:" << FCAL_PHOT_STAT_COEF << ")" << endl;
//   cout << "    -b#      FCAL single block threshold in MeV (def:" << FCAL_BLOCK_THRESHOLD/k_MeV << ")" << endl;
//   cout << "    -B       Don't process BCAL hits at all (def. process)" << endl;
 //  cout << "    -Vthresh BCAL ADC threshold (def. " << BCAL_ADC_THRESHOLD_MEV << " MeV)" << endl;
 //  cout << "    -Xsigma  BCAL fADC time resolution (def. " << BCAL_FADC_TIME_RESOLUTION << " ns)" << endl;
 //  cout << "    -G       Don't smear BCAL times (def. smear)" << endl;
 //  cout << "    -H       Don't add BCAL dark hits (def. add)" << endl;
 //  cout << "    -K       Don't apply BCAL sampling fluctuations (def. apply)" << endl;
 //  cout << "    -L       Don't apply BCAL sampling floor term (def. apply)" << endl;
 //  cout << "    -M       Don't apply BCAL Poisson statistics (def. apply)" << endl;
 //  cout << "    -f#      TOF sigma in psec (def: " <<  TOF_SIGMA/k_psec << ")" << endl;
   cout << "    -h       Print this usage statement." << endl;
   cout << endl;
//   cout << " Example:" << endl;
//   cout << endl;
//   cout << "     mcsmear -u3.5 -t500 hdgeant.hddm" << endl;
//   cout << endl;
//   cout << " This will produce a file named hdgeant_nsmeared.hddm that" << endl;
//   cout << " includes the hit information from the input file hdgeant.hddm" << endl;
//   cout << " but with the FDC and CDC hits smeared out. The CDC hits will" << endl;
//   cout << " have their drift times smeared via a gaussian with a 3.5ns width" << endl;
//   cout << " while the FDC will be smeared using the default values." << endl;
//   cout << " In addition, background hits will be added, the exact number of" << endl;
//   cout << " of which are determined by the time windows specified for the" << endl;
//   cout << " CDC and FDC. In this examplem the CDC time window was explicitly" << endl;
//   cout << " set to 500 ns." << endl;
//   cout << endl;

   exit(0);
}

