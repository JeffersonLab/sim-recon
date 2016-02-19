#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMatrixT.h"
#include "TMatrixTSym.h"
#include "TMatrixDSymEigen.h"
#include "TVectorT.h"

#include <iostream>    

using namespace std;

int main( int argc, char** argv ){

if(string(argv[1]) == "help" || string(argv[1]) == "-help" || string(argv[1]) == "--help" || string(argv[1]) == "-h") {
    cout << "This program expects the following argument(s) in order: " << endl;
    cout << "Argument 1: path to merged gain_factor_matrices.root" << endl;
    return 0;
}

if(argc != 2) {
    cout << "Error: you must specify the path to your gain_factor_matrices.root file!!!" << endl;
    return 0;
}

TFile* myfile = new TFile(argv[1]);
TH1F* h1D_mC = (TH1F*)myfile->Get("h2D_mC");

int n_channels = h1D_mC->GetNbinsX();
int n_channelsy = h1D_mC->GetNbinsY();

 if (n_channels != n_channelsy) {
   cout << "*** Non-square matrix***, nchannels=" << n_channels << " nchannelsy=" << n_channelsy << endl;
   return -1;
 }
 else {
   cout << "Found square matrix m_mc, nchannels=" << n_channels << endl;
 }

TMatrixDSym m_mC;

m_mC.ResizeTo(n_channels,n_channels);
    //Generate Matrices back from TH1F and TH2F's

cout << "Generating Matrix from imported TH1F..." << endl;
for(int i = 0  ; i<n_channels; ++i) {
     for( int j = 0 ; j<n_channels; ++j) {
         m_mC[i][j]=h1D_mC->GetBinContent(i+1,j+1);
     }
}


cout << "Getting Eigensystem (this will take a while...) \n";
TMatrixDSymEigen mC_eigen(m_mC);
TVectorD eval = mC_eigen.GetEigenValues();
TMatrixD evec = mC_eigen.GetEigenVectors();
cout << "Filling Histograms... \n";


double mymin=eval.Min();
double mymax=eval.Max();

  TH1F* h1_evals = new TH1F( "h1_evals", "Eigenvalues of C Matrix by Channel Number",n_channels, 0., n_channels );
  TH2F* h2_evecs = new TH2F( "h2_evecs DO NOT PLOT", "Eigenvectors (good luck drawing this)", n_channels, 0., n_channels, n_channels,0.,n_channels );
  TH1F* eig = new TH1F( "Eigenvalues, unordered", "Eigenvalues of C matrix" ,100, mymin, mymax/12. );



for(int i = 0  ; i<n_channels; ++i) {

   if(eval[i]>0.0 &&eval[i]<10000.0) h1_evals->SetBinContent(i+1,eval[i]);
    eig->Fill(eval[i]);
    for(int j = 0  ; j<n_channels; ++j) {
         h2_evecs->SetBinContent(i+1,j+1,evec[i][j]);
    }
}

  TFile* m_rootFile = new TFile( "C_Eigensystem.root", "RECREATE" );

  m_rootFile->cd();
  h1_evals->Write();
  h2_evecs->Write();
  eig->Write();
  m_rootFile->Close();


cout << "done" << endl;

  return 0;

}
