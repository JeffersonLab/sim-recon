// $Id: mcsmear.cc 19023 2015-07-14 20:23:27Z beattite $
//
// Created June 22, 2005  David Lawrence

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>

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
#include "mcsmear_config.h" 

#include "units.h"
#include "HDDM/hddm_s.hpp"

void Smear(hddm_s::HDDM *record);
void ParseCommandLineArguments(int narg, char* argv[], mcsmear_config_t *in_config);
void Usage(void);

extern void SetSeeds(const char *vals);

char *INFILENAME = NULL;
char *OUTFILENAME = NULL;
int QUIT = 0;

std::map<hddm_s::istream*,double> files2merge;
std::map<hddm_s::istream*,hddm_s::streamposition> start2merge;
std::map<hddm_s::istream*,int> skip2merge;

using namespace jana;

// for histogramming
//pthread_mutex_t root_mutex = PTHREAD_MUTEX_INITIALIZER;

// GLOBAL RANDOM NUMBER GENERATOR
// Note, the argument is zero to cause the seeds to
// be initialized using the UUID (see code for ROOT's
// TRandom2 constructor) No argument, or an argument 
// greater than zero will result in the same seeds 
// being set every time mcsmear is run.
DRandom2 gDRandom(0); // declared extern in DRandom2.h


//-----------
// main
//-----------
int main(int narg,char* argv[])
{
   mcsmear_config_t *config = new mcsmear_config_t();
   ParseCommandLineArguments(narg, argv, config);

   // Create DApplication object
   DApplication dapp(narg, argv);
   dapp.AddFactoryGenerator(new JFactoryGenerator_ThreadCancelHandler());

   TFile *hfile = new TFile("smear.root","RECREATE","smearing histograms");  // note: not used for anything right now

   MyProcessor myproc(config);   
   jerror_t error_code = dapp.Run(&myproc);

   hfile->Write();
   hfile->Close();

   if(error_code != NOERROR) 
       return static_cast<int>(error_code);
   else
       return 0;
}

//-----------
// ParseCommandLineArguments
//-----------
void ParseCommandLineArguments(int narg, char* argv[], mcsmear_config_t *config)
{

   for (int i=1; i<narg; i++) {
      char *ptr = argv[i];
    
      if (ptr[0] == '-') {
         switch(ptr[1]) {
          case 'h': Usage();                                     break;
          case 'o': OUTFILENAME = strdup(&ptr[2]);               break;
          case 'N': config->ADD_NOISE=true;                      break;
          case 's': config->SMEAR_HITS=false;                    break;
          case 'i': config->IGNORE_SEEDS=true;                   break;
          case 'r': config->SetSeeds(&ptr[2]);                   break;
          case 'd': config->DROP_TRUTH_HITS=true;                break;
          case 'D': config->DUMP_RCDB_CONFIG=true;               break;
          case 'e': config->APPLY_EFFICIENCY_CORRECTIONS=false;  break;
          case 'E': config->FCAL_ADD_LIGHTGUIDE_HITS=true;       break;

	      // BCAL parameters
          case 'G': config->BCAL_NO_T_SMEAR = true;              break;
          case 'H': config->BCAL_NO_DARK_PULSES = true;          break;
          case 'K': config->BCAL_NO_SAMPLING_FLUCTUATIONS = true; break;
          case 'L': config->BCAL_NO_SAMPLING_FLOOR_TERM = true;  break;
          case 'M': config->BCAL_NO_POISSON_STATISTICS = true;   break;
	  case 'S': config->BCAL_NO_FADC_SATURATION = true;      break;
         }
      }
      else {
         std::string filename(ptr);
         size_t colon = filename.find_first_of(":");
         if (colon != filename.npos) {
            double wgt = std::stod(filename.substr(colon + 1));
            size_t plus = filename.substr(colon + 1).find_first_of("+");
            size_t decimal = filename.substr(colon + 1, plus).find_first_of(".");
            if (decimal != filename.npos) // distinguish float from int
               wgt += 1e-10;
            int skip = 0;
            if (plus != filename.npos)
               skip = std::stoi(filename.substr(colon + plus + 1));
            std::ifstream fin(filename.substr(0, colon));
            hddm_s::istream stin(fin);
            hddm_s::HDDM record;
            stin >> record;
            std::ifstream *ifs = new std::ifstream(filename.substr(0, colon));
            hddm_s::istream *istr = new hddm_s::istream(*ifs);
            start2merge[istr] = stin.getPosition();
            files2merge[istr] = wgt;
            skip2merge[istr] = skip;
            std::fill(ptr, ptr + strlen(ptr), '-');
            continue;
         }
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
      sprintf(str, "%s_smeared.hddm", path_stripped);
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
   cout << "     mcsmear [options] file.hddm [noise1.hddm:<N1> [...] ]" << endl;
   cout << endl;
   cout << "Read the given, Geant-produced HDDM file as input and smear" << endl;
   cout << "the truth values for \"hit\" data before writing out to a" << endl;
   cout << "separate file. The truth values for the thrown particles are" << endl;
   cout << "not changed. Noise hits can also be added appending additional" << endl;
   cout << "input hddm files after the primary input file, denoted above" << endl;
   cout << "as noise1.hddm:<N1>. Each event in the primary input file will" << endl;
   cout << "be merged at hits level with <N1> events from the first listed" << endl;
   cout << "noise file, <N2> events from the second noise file, and so on" << endl;
   cout << "for as many noise files as are listed. If the pileup factor <N>" << endl;
   cout << "is a float (contains a decimal point) then the number of events" << endl;
   cout << "from the noise file that get merged into each event in the" << endl;
   cout << "primary input file is generated at random from a Poisson" << endl;
   cout << "distribution with a mean of <N>. When all of the input events" << endl;
   cout << "in any of the noise files are exhausted, the file is opened" << endl;
   cout << "again and reading of noise events restarts from the beginning" << endl;
   cout << "of the file. If you want to skip S events at the beginning of" << endl;
   cout << "the noise file at startup, append \"+S\" to the <N> argument." << endl;
   cout << "Note that all smearing is done using Gaussians." << endl;
   cout << endl;
   cout << "  options:" << endl;
   cout << "    -ofname  Write output to a file named \"fname\" (default auto-generate name)" << endl;
   cout << "    -s       Don't smear real hits (default is to smear)" << endl;
   cout << "    -i       Ignore random number seeds found in input HDDM file" << endl;
   cout << "    -r\"s1 s2 s3\" Set initial random number seeds" << endl;
   cout << "    -e       Don't apply channel dependent efficiency corrections" << endl;
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
   cout << "    -D       Dump configuration debug information" << endl;
   cout << "    -G       Don't smear BCAL times (def. smear)" << endl;
   cout << "    -H       Don't add BCAL dark hits (def. add)" << endl;
   cout << "    -K       Don't apply BCAL sampling fluctuations (def. apply)" << endl;
   cout << "    -L       Don't apply BCAL sampling floor term (def. apply)" << endl;
   cout << "    -M       Don't apply BCAL Poisson statistics (def. apply)" << endl;
   cout << "    -S       Don't apply BCAL fADC saturation (def. apply)" << endl;
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
