
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMatrixT.h"
#include "TMatrixTSym.h"
#include "TMatrixDSymEigen.h"
#include "TVectorT.h"

#include <string>
#include <fstream>
#include <sstream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "BCAL/DBCALGeometry.h"
#include "FCAL/DFCALGeometry.h"    

DFCALGeometry fcalgeom;

using namespace std;

pair<int,int> AbsNumtoXY(int channel);
vector<vector<string> > ParseTSV(const char* s);
TMatrixD GetNewGains(TMatrixD oldmatrix, string lastfile);


int main( int argc, char** argv ){
  int hits_cut = 10;

if(string(argv[1]) == "help" || string(argv[1]) == "-help" || string(argv[1]) == "--help" || string(argv[1]) == "-h") {
    cout << "This program expects the following argument(s) in order: " << endl;
    cout << "Argument 1: path to root file that contains C,D,L matrices as histos" << endl;
    cout << "Argument 2: path to root file containing eigenvalues and eigenvectors" << endl;
    cout << "Argument 3: lower limit cut on eigenvalues of C matrix" << endl;
    cout << "Argument 4: old gain factor text file (if not first iteration)" << endl;
    return 0;
}

else if(argc != 4 && argc != 5) {
    cout << "Error: unexpected number of arguments" << endl;
    return 0;
}

string lastfilename;
string nextfilename;
string nextroot;

vector<vector<string> > oldgains;

if(argc == 4) {
    cout << "Current iteration is: 1" << endl;
     nextroot = "GainFactorIter1.root";
     nextfilename = "GainFactorIter1.txt";
}

if(argc ==5){
    char* textfilename = argv[4];
    const char* it = &textfilename[14];
    int iteration = atoi(it)+1;


    if(iteration == 1) {
        cout << "Error: Something's wrong with GainFactorIter%i.txt file!!" << endl;
        return 0;
    }




    cout << "Current iteration is " << iteration << endl;
    char iter[50]; char nexty[50];
    sprintf(iter,"GainFactorIter%i.txt",iteration);
    sprintf(nexty,"GainFactorIter%i.root",iteration);
    nextfilename = iter; nextroot = nexty;
    cout << "Next file name: " << nextfilename << endl;
    lastfilename = argv[4];
    if(iteration == 10) {
        cout << "if it didnt work by now it never will :(" << endl;
        return 0;
    }
 }

TFile* eigfile = new TFile(argv[2]);
TFile* matrixfile  = new TFile(argv[1]);

TH1F* h1D_mD = (TH1F*)matrixfile->Get("h1D_mD");
TH1F* h1D_nhits = (TH1F*)matrixfile->Get("h1D_nhits");

TH2F* h2_evecs = (TH2F*)eigfile->Get("h2_evecs DO NOT PLOT");
TH1F* h1_evals = (TH1F*)eigfile->Get("h1_evals");

int n_channels = h1D_mD->GetNbinsX();

 int m_nElements = n_channels;

 cout << "Found historgram h1D_mD, dimension=" << n_channels << endl;


 // define new histograms

  TH1F* h_GainFactors = new TH1F( "h_GainFactors", "Gain Factor as Function of Channel Number" ,n_channels, 0., n_channels );
  TH1F* h_GainFactors_old = new TH1F( "h_GainFactors_old", "Old Gain Factor as Function of Channel Number" ,n_channels, 0., n_channels );
  TH1F* h_unordered_new = new TH1F( "h_unordered_new", "Gain Factors Unordered" ,100, 0., 2. );
  TH1F* h_unordered = new TH1F( "h_unordered", "Gain Factors Unordered" ,100, 0., 2. );
  TH1F* h_unordered_diff = new TH1F( "h_unordered_diff", "Old gain factor - New Gain Factor" ,100, -1., 1. );
  TH1F* h_old_unordered = new TH1F( "h_old_unordered", "Old Gain Factors Unordered" ,100, 0., 2. );
  TH1F* h_GainFactorsNew = new TH1F("h_GainFactorsNew", "Actual gain factors, ordered",n_channels,0,n_channels);
  TH2D* h_yvsx = new TH2D("h_yvsx", "Gain factors y vs x",60,0,60,60,0,60);


TMatrixD m_mD;
TMatrixD m_mL;
TMatrixD m_mLt;

m_mD.ResizeTo(m_nElements,1);
m_mL.ResizeTo(m_nElements,1);
m_mLt.ResizeTo(1,m_nElements);
m_mD.Zero();
m_mL.Zero();

TMatrixD evec;
evec.ResizeTo(n_channels,n_channels);
evec.Zero();

TVectorD eval;
eval.ResizeTo(n_channels);
eval.Zero();

TMatrixD eiginv;
eiginv.ResizeTo(n_channels,n_channels);
eiginv.Zero();

TMatrixD evecinv;
evecinv.ResizeTo(n_channels,n_channels);
evecinv.Zero();

double cutoff = atof(argv[3]);

cout << "Generating Matrices..." << endl;
for(int i = 0; i<n_channels; ++i) {
     m_mD[i][0]=h1D_mD->GetBinContent(i+1);
    if(h1_evals->GetBinContent(i+1) > cutoff) eiginv[i][i] = (1./h1_evals->GetBinContent(i+1));
    for(int j = 0  ; j<n_channels; ++j) {
       evec[i][j] = h2_evecs->GetBinContent(i+1,j+1);
    }
}

evecinv.Transpose(evec);

cout << "Doing some matrix multiplication now (this could take a while...)" << endl;
TMatrixD C_inv = evec * eiginv * evecinv;


//Now Generate gain factor vector // Comment out because gain factors c_k are not 1+eps_k
TMatrixD GainFactors = C_inv*m_mD;
/*for(int i = 0; i<n_channels; ++i) {
    GainFactors[i][0]+=1;
    }*/

if(argc == 5) {
  oldgains = ParseTSV(lastfilename.c_str());


    for(int i = 0; i<n_channels; ++i) {
        h_GainFactorsNew->SetBinContent(i+1,GainFactors[i][0]);
        h_unordered_new->Fill(GainFactors[i][0]);
    }



    GainFactors = GetNewGains(GainFactors, lastfilename);
}


ofstream next_file;
next_file.open(nextfilename.c_str());
int nskipped =0;


  for(int i=0; i<n_channels; i++){
     if(h1D_nhits->GetBinContent(i+1) < hits_cut) {
            nskipped++;
            GainFactors[i][0] = 1;
  }
  }


  for(int i=0; i<n_channels; i++){

    double my_gain = GainFactors[i][0];
    ostringstream sstream; sstream << my_gain;
    string smy_gain = sstream.str();
    ostringstream channel; channel<<i+1;

   // smy_gain = smy_gain+"\t"+channel.str();
	
    // next_file << smy_gain << endl;
    // add information to file
    int row = fcalgeom.row(i);
    int column = fcalgeom.column(i);
    DVector2 pos = fcalgeom.positionOnFace(i);
    double r = sqrt(pos.X()*pos.X() + pos.Y()*pos.Y());
    next_file  << smy_gain << "  channel=" << i << " row=" << row << " column=" << column << " x=" << pos.X() << " y=" << pos.Y() << " r=" << r << endl;
    if (my_gain > 0) h_yvsx->Fill(column,row,my_gain);    
    
  }

next_file.close();

//Now export to a 1D Histogram and save text file



  for(int i=0; i<n_channels; i++){

    if(argc ==5) {

        double oldgain = atof(oldgains[i][0].c_str());
        h_GainFactors_old->SetBinContent(i+1,oldgain);
        h_old_unordered->Fill(oldgain);
        h_unordered_diff->Fill(oldgain-GainFactors[i][0]);
    }

    h_GainFactors->SetBinContent(i+1,GainFactors[i][0]);
    h_unordered->Fill(GainFactors[i][0]);
  }


  TFile* m_rootFile = new TFile( nextroot.c_str(), "RECREATE" );
  m_rootFile->cd();
  h_GainFactors->Write();
  h_GainFactorsNew->Write();
  h_unordered_new->Write();
  h_GainFactors_old->Write();
  h_unordered->Write();
  h_old_unordered->Write();
  h_unordered_diff->Write();
  h_yvsx->Write();
  m_rootFile->Close();



cout << "Total Channels skipped: " << nskipped << endl;

cout << "done" << endl;

  return 0;

}



vector<vector<string> > ParseTSV(const char* s) {
typedef vector<vector<string> > Rows;
Rows parsed;
ifstream input(s);
char const row_delim = '\n';
char const field_delim = '\t';
for (string row; getline(input, row, row_delim); ) {
  parsed.push_back(Rows::value_type());
  istringstream ss(row);
  for (string field; getline(ss, field, field_delim); ) {
    parsed.back().push_back(field);
  }
}

return parsed;
}

TMatrixD GetNewGains(TMatrixD oldmatrix, string lastfile){

  int n_channels = 2800;

    TMatrixD newmatrix;
    newmatrix.ResizeTo(n_channels,1);

    vector<vector<string> > oldgains = ParseTSV(lastfile.c_str());

    for(int i = 0; i<n_channels; ++i) {
        double lastgain = atof(oldgains[i][0].c_str());
        newmatrix[i][0] = oldmatrix[i][0]*lastgain;
    }

    return newmatrix;
}

