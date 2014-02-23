// standard C++ header files 
#include <iostream>
#include <fstream>

// ROOT header files
#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"

// HallD header files
#include "kinfit/KinFit.h"

void PrintOptions(const char*);

using namespace std;

int main(int argc, char *argv[]) {

  int entry,Nevents,c;
  char *progName = argv[0];
  char *outFile = "default.root";
  string inFile;
  int max = (int)1E9;
  int verbose = 0;
  string trackConv = "EASY";
  extern char* optarg;
  extern int optind;
  int isMC = 0;
  char name[256];
  double E, px, py, pz;
  string junk;

  float mm;

  TLorentzVector postKfitPhoton, preKfitPhoton, target, p, pip, pim, missm, pi0, kp;
  TLorentzVector preMM, postMM, twoPions, phiMeson;
  TLorentzVector ppim;
  TLorentzVector dummyP4;
  int sector[3];

  target.SetXYZM(0, 0, 0, 0.938272);

  string pName[3];
  pName[0] = "proton";
  pName[1] = "pi+";
  pName[2] = "pi-";

  // Check command line arguments
  if(argc == 1){
    PrintOptions(progName);
    return 0;
  }
  else{   
    while((c = getopt(argc,argv,"hm:o:vi:t:")) != -1){
      switch(c){
        case 'h':
          PrintOptions(progName);
          return 0;
          break;
        case 'm':
          max = atoi(optarg);
          cout << "Max Events: " << max << endl;
          break;
        case 'o': 
          outFile = optarg;
          cout << "Output file name: " << outFile << endl;
          break;
        case 'i': 
          inFile = optarg;
          cout << "Input file name: " << inFile << endl;
          break;
        case 't':
          trackConv = optarg;
          break;
        case 'v':
          verbose++;
          break;
        default:
          break;
      }
    }
  } 

  cout << "Tracking conversion: " << trackConv << endl;
  // Set up the ROOT environment
  TROOT simple("simple","");
  TFile rfile(outFile,"RECREATE");
  TTree *tree = new TTree("T","analysis");

  // Declare some histograms
  TH1F *hprob[4]; // Confidence level of fit
  TH1F *hpulls[4][10]; // Pulls from the fit.
  TH1F *hmm[4]; // Missing mass
  TH1F *h2pi[4]; // Invariant mass of 2 pions
  TH1F *h3mesons[4]; // Invariant mass of 2 pions and missing mass

  // Initialize some histograms
  for(int i=0;i<4;i++)
  {
    sprintf(name,"hmm%d",i);
    hmm[i] = new TH1F(name,"",100,-0.1,0.4);

    sprintf(name,"hprob%d",i);
    hprob[i] = new TH1F(name,"",100,0.0,1.0);

    sprintf(name,"h2pi%d",i);
    h2pi[i] = new TH1F(name,"",90,0.3,1.2);

    sprintf(name,"h3mesons%d",i);
    h3mesons[i] = new TH1F(name,"",80,0.5,1.3);

    for(int j=0;j<10;j++)
    {
      sprintf(name,"hpulls%d_%d",i,j);
      hpulls[i][j] = new TH1F(name,"",100,-5.0,5.0);
    }
  }

  // Delcare a few items that we will need for the 
  // kinematic fit.
  // Note that the covariance matrix starts out as a 1x1 matrix
  // which would be the error on the photon energy for now.
  TMatrixD cov_matrix(1,1);
  double EphotError= 0.0;
  vector<TLorentzVector> preKfitP4(3);
  vector<TLorentzVector> postKfitP4(3);
  vector<string> reconstruction(3);
  reconstruction[0] = "DC";
  reconstruction[1] = "DC";
  reconstruction[2] = "DC";

  // Open the input file
  cerr << "Processing file: " << inFile << endl;
  ifstream IN(inFile.c_str());

  entry = 0; // reset entry number for this file

  /*____________________Loop Over Events in Current File___________________*/
  Nevents = 0;
  while(IN >> junk && Nevents < max)
  {
    if(Nevents%1000==0) cerr << Nevents << "\r";
    Nevents++;

    // Declare a new kinematic fit and set up the 
    // covariance matrix for a new fit.
    KinFit *kfit = new KinFit();
    cov_matrix.ResizeTo(1,1);

    // Read in the events from the file.
    // Read in the photon information
    IN >> E >> px >> py >> pz;
    IN >> EphotError;
    preKfitPhoton.SetXYZT(px, py, E, E);
    // Start to set the covariance matrix
    cov_matrix(0,0) = EphotError;
    // Read in the final state particles
    for(int i=0;i<3;i++)
    {
      IN >> pName[i];
      IN >> E >> px >> py >> pz; // 4-momentum
      // Set the 4-vectors to be passed in to the kinematic fit
      dummyP4.SetXYZT(px, py, pz, E);
      preKfitP4[i] = dummyP4;
      // Our covariance matrix is block diagonal.
      // Resize the covariance matrix by adding another block of 3x3 for the three 
      // tracking parameters for each particle, and read in the errors. 
      cov_matrix.ResizeTo(cov_matrix.GetNrows()+3, cov_matrix.GetNcols()+3);
      for(int j=0;j<3;j++) // Read in the errors on the tracking parameters
      {
        for(int k=0;k<3;k++)
        {
          IN >> cov_matrix(3*i+1+j,3*i+1+k);
        }
      } 
    }

    // Grab the mising mass and two pi mass PRIOR to fitting.
    preMM = preKfitPhoton + target - preKfitP4[0] - preKfitP4[1] - preKfitP4[2];
    twoPions = preKfitP4[1] + preKfitP4[2];

    // Fill a histogram with the original missing mass squared.
    hmm[0]->Fill(preMM.M2());

    // Loop over three different hypothesis 
    // for this final state.
    for(int trial=1;trial<=3;trial++)
    {
      // Set up the kimatic fit
      kfit->SetPhotonEnergy(preKfitPhoton.E());
      kfit->SetP4(preKfitP4);
      kfit->SetReconstruction(reconstruction); // Tells how the track was reconstructed (DC, Calorimetry, etc)
      kfit->SetCovMatrix(cov_matrix);
      kfit->SetTargetMass(target.M());
      kfit->SetTrackingConversion(trackConv);
      // Could also call
      //kfit->SetEvent(preKfitPhoton.E(), preKfitP4, reconstruction, cov_matrix, target.M(), reconstruction);
      // Start the fits!
      if(trial==1) 
      {
        // Fit to nothing missing. 4-C fit
        kfit->Fit();
      }
      else if(trial==2) 
      {
        // Fit to a missing pi0. 1-C fit
        kfit->Fit("pi0");
      }
      else if(trial==3) 
      {
        // Fit to a missing K0 and constrain the pion
        // mass to be that of another K0. 2-C fit.
        vector<bool> constrain_measured(3); // Used to determine which measured
                                            // particles are in the extra mass
                                            // constraint.
        bool constrain_missing; // Is the missing particle used in the extra mass
                                // constraint.
        constrain_measured[0] = false; // proton
        constrain_measured[1] = true; // pi+
        constrain_measured[2] = true; // pi-
        constrain_missing     = false; // missing K in this case.
        kfit->Fit(0.49767, constrain_measured, constrain_missing, 0.49767);
      }

      // Now that the fit is over, let's grab some information 
      // from the fit.
      // First grab the wiggled photon energy.
      double fitPhotonE = kfit->FitPhotonEnergy();
      for(int i=0;i<3;i++)
      {
        // Grab the wiggled 4-vectors
        postKfitP4[i] = kfit->FitP4(i);
      }
      double prob = kfit->Prob(); // Confidence level.
      if(verbose && trial==1) cerr << "prob: " << prob << endl;

      // Calculate a few post-fit quantities
      postKfitPhoton.SetXYZT(0, 0, fitPhotonE, fitPhotonE);
      postMM = postKfitPhoton + target - postKfitP4[0] - postKfitP4[1] - postKfitP4[2];
      phiMeson = (postKfitP4[1] + postKfitP4[2]) + postMM;

      hprob[trial]->Fill(prob);

      if(prob>0.01)
      {
        hmm[trial]->Fill(preMM.M2()); // Fill with the prefit missing mass,
                                      // but only if it passes a 1% CL cut.
        h2pi[trial]->Fill(twoPions.M()); // Fill with the prefit two pi mass,
                                         // but only if it passes a 1% CL cut.
        h3mesons[trial]->Fill(phiMeson.M()); // Fill with the post fit mass of the 
                                             // two pi's and missing mass, but only
                                             // if if passes the 1% CL cut.
        for(int i=0;i<10;i++)
        {
          // Grab the pulls on each of the fit quantities, but only
          // if it passes the 1% CL cut.
          hpulls[trial][i]->Fill(kfit->GetPull(i));
        }
      }
    }
    tree->Fill();
    delete kfit;
  }
  cout << "Events Processed: " << Nevents << "." << endl;
  // Write out to outFile
  rfile.Write();
  tree->Print();  
  rfile.Close();

  return 0; 
}
//_____________________________________________________________________________

void PrintOptions(const char *name) {
  // Print options to the screen.
  cout << "Usage:" << endl;
  cout << name << " [option 1] [option 2] " << endl;
  cout << "\t-m[#] Max number of events to be processed." << endl;
  cout << "\t-i<inFile> Input file name." << endl;
  cout << "\t-o<outFile> Output file name." << endl;
  cout << "\t-t<trackingConversion> How to convert momentum to tracking parameters. (EASY, CLAS)" << endl;
  cout << "\t-v  Verbose (prints to screen events written to output TTree)" << endl;
  cout << "\t-h  Print this message." << endl;
  cout << endl;
}
