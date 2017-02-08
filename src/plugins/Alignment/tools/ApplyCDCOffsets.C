// Script to get current values of the offsets from the alignment root file, and shift them according to the MILLEPEDE result.
#include <map>
#include <fstream>
#include <iostream>
#include <string>
#include "TProfile.h"
#include "TFile.h"

using namespace std;

void ApplyCDCOffsets(TString rootFile = "hd_root.root", TString pedeOutFile = "millepede.res"){
  
   TFile *file = TFile::Open(rootFile);
   TProfile *hCDCConstants;
   file->GetObject("AlignmentConstants/CDCAlignmentConstants",hCDCConstants);

   ifstream pedeResult;
   pedeResult.open(pedeOutFile.Data());

   map < int, double> resultMap;

   bool firstLine = true;
   string line;
   while (getline(pedeResult,line)){
      if (firstLine){
         firstLine = false;
         continue;
      }
      int index; double value;
      istringstream iss(line);
      iss >> index >> value;

      resultMap[index]=value;
   }

  // Global CDC offsets 
   ofstream outFile;
   outFile.open("cdc_global.txt");
   for (unsigned int i = 1; i<=6; i++){
      int histIndex = i;
      int pedeIndex = i;
      outFile << -resultMap[pedeIndex] + hCDCConstants->GetBinContent(histIndex) << " " ;
   }
   outFile.close();

   // Wire alignment
   outFile.open("cdc_wire_alignment.txt");
   for (unsigned int i = 1001; i<=15088; i++){
      int histIndex = i;
      int pedeIndex = i;
      outFile << -resultMap[pedeIndex] + hCDCConstants->GetBinContent(histIndex) << " " ;
      if (i%4 == 0) outFile << endl;
   }
   outFile.close();

}
