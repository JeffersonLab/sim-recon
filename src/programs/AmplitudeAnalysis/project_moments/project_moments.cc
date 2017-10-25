#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <unistd.h>

#include "IUAmpTools/FitResults.h"

#include "wave.h"
#include "3j.h"

#include "TFile.h"

using namespace std;

int main( int argc, char* argv[] ){

  // set default parameters
  
  double lowMass = 0.28;
  double highMass = 2.00;
  int kNumBins = 86;
  string fitDir("");
  
  string outfileName("moments.root");
  bool print = false;

  // parse command line
  
  for (int i = 1; i < argc; i++){
    
    string arg(argv[i]);
    
    if (arg == "-f"){  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  fitDir = argv[++i]; }
    if (arg == "-o"){  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  outfileName = argv[++i]; }
    if (arg == "-p")  
      print = true;
    if (arg == "-h"){
      cout << endl << " Usage for: " << argv[0] << endl << endl;
      cout << "\t   -f <fit dir>\t : Fit Directory" << endl;
      cout << "(optional) -o <file>\t : Output file (default: moments.root)" << endl;
      cout << "(optional) -p\t\t : Print equations" << endl;
      exit(1);}
  }
  
  if (fitDir.size() == 0){
    cout << "No fit directory specified. Try -h for usage." << endl;
    exit(1);
  }
  
  TFile *outfile = new TFile(outfileName.c_str(), "recreate");
  if (!outfile->IsOpen()) exit(1);

  // Set waveset, has to be same order as in fit.cfg!

  vector<wave> negative;
  negative.push_back(wave("S0", 0, 0));
  negative.push_back(wave("P0", 1, 0));
  negative.push_back(wave("P-", 1, 1));
  negative.push_back(wave("D0", 2, 0));
  negative.push_back(wave("D-", 2, 1));

  vector<wave> positive;
  positive.push_back(wave("P+", 1, 1));
  positive.push_back(wave("D+", 2, 1));

  coherent_waves wsPos, wsNeg;
  wsPos.reflectivity = +1;
  wsPos.waves = positive;

  wsNeg.reflectivity = -1;
  wsNeg.waves = negative;

  waveset ws;
  ws.push_back(wsNeg);
  ws.push_back(wsPos);

  size_t lastIdx = 0;
  for (size_t i = 0; i < ws.size(); i++)
    for (size_t j = 0; j < ws[i].waves.size(); j++, lastIdx += 2)
      ws[i].waves[j].setIndex(lastIdx);

  // Find the non-zero moments for the given waveset.
  std::vector<std::pair<size_t, size_t> > vecMom = listOfMoments(ws);

  // Prepare the moment histograms and print moment equations
  map<std::pair<size_t,size_t>, TH1D*> hMoments;
  for (std::vector<std::pair<size_t, size_t> >::const_iterator it = vecMom.begin(); it != vecMom.end(); it++)
    {
      
      char name[999];
      char title[999];
      snprintf(name, 999, "hMoment%zd%zd", it->first, it->second);
      snprintf(title, 999, "Moment H(%zd,%zd)", it->first, it->second);
      hMoments[*it] = new TH1D(name, title, kNumBins, lowMass, highMass);

      if (print){
	cout << "H(" << it->first << ", " << it->second << ") = ";
	for (size_t iCoh = 0; iCoh < ws.size(); iCoh++)
	  {
	    const coherent_waves& waves = ws[iCoh];
	    int eps = waves.reflectivity;
	    
	    for (size_t i = 0; i < waves.waves.size(); i++)
	      {
		const wave& w1 = waves.waves[i];
		for (size_t j = 0; j <= i; j++)
		  {
		    const wave& w2 = waves.waves[j];
		    
		    double coeff = getCoefficient(eps, it->first, it->second, w1.getL(), w1.getM(), w2.getL(), w2.getM());
		    if (coeff == 0)
		      continue;
		    
		    if (j == i)
		      {
			if (coeff > 0)
			  cout << "+ " << coeff;
			else if (coeff < 0)
			  cout << "- " << -coeff;
			cout << " |" << w1.getName() << "|^2 ";
		      }
		    else
		      {
			if (coeff > 0)
			  cout << "+ " << 2*coeff;
			else if (coeff < 0)
			  cout << "- " << -2*coeff;
			cout << " Re(" << w1.getName() << "* " << w2.getName() << ") ";
		      }
		  }
	      }
	  }
	cout << endl;
      }
    }
  
  
  // descend into the directory that contains the bins
  chdir( fitDir.c_str() );
  
  for( int i = 0; i < kNumBins; ++i ){
    
    ostringstream dir;
    dir << "bin_" << i;
    chdir( dir.str().c_str() );
    
    ostringstream resultsFile;
    resultsFile << "bin_" << i << ".fit";
    
    FitResults results( resultsFile.str() );
    if( !results.valid() ){
      
      chdir( ".." );
      continue;
    }

    if (  2*ws.getNwaves() != results.parValueList().size() ){
      cout << "Different number of waves in fit result. Check waveset!" << endl;
      outfile->Close();
      exit(1);
    }
    
    for (std::vector<std::pair<size_t, size_t> >::const_iterator it = vecMom.begin(); it != vecMom.end(); it++)
      {
	hMoments[*it]->SetBinContent(i + 1, decomposeMoment(*it, ws, results.parValueList()));
	hMoments[*it]->SetBinError(i + 1, decomposeMomentError(*it, ws, results.parValueList(), results.errorMatrix()));
      }
    
    chdir( ".." );
  }
  
  
  for (std::vector<std::pair<size_t, size_t> >::const_iterator it = vecMom.begin(); it != vecMom.end(); it++)
    {
      hMoments[*it]->Write();
    }
  
  outfile->Close();
  cout << "Moment histograms written to " << outfileName << endl;
  
  return 0;
}
