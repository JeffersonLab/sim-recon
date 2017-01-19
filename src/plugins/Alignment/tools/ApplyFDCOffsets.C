// Script to get current values of the offsets from the alignment root file, and shift them according to the MILLEPEDE result.
#include <map>
#include <fstream>
#include <iostream>
#include <string>
#include "TProfile.h"
#include "TFile.h"

using namespace std;

void ApplyFDCOffsets(TString rootFile = "hd_root.root", TString pedeOutFile = "millepede.res"){
  
   TFile *file = TFile::Open(rootFile);
   TProfile *hFDCConstants;
   file->GetObject("AlignmentConstants/FDCAlignmentConstants",hFDCConstants);

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

  // Cell offsets 
  // Indices are 101001, 101100 etc...
   ofstream outFile;
   outFile.open("cell_offsets.txt");
   for (unsigned int i = 1; i<=24; i++){
      int histIndex = i*1000;
      int pedeIndex = 100000 + i*1000;
      outFile << resultMap[pedeIndex + 1] + hFDCConstants->GetBinContent(histIndex+1) << " " ;
      outFile << resultMap[pedeIndex + 100] + hFDCConstants->GetBinContent(histIndex+100) << endl ;
   }
   outFile.close();

   // Cathode Offsets
   // Indices are 101, 102, 103, 104
   outFile.open("cathode_alignment.txt");
   for (unsigned int i = 1; i<=24; i++){
      int histIndex = i*1000;
      int pedeIndex = 100000 + i*1000;
      outFile << resultMap[pedeIndex + 103] + hFDCConstants->GetBinContent(histIndex+103) << " " ;//dPhiU
      outFile << resultMap[pedeIndex + 101] + hFDCConstants->GetBinContent(histIndex+101) << " " ;//dU
      outFile << resultMap[pedeIndex + 104] + hFDCConstants->GetBinContent(histIndex+104) << " " ;//dPhiV
      outFile << resultMap[pedeIndex + 100] + hFDCConstants->GetBinContent(histIndex+102) << endl ;//dV
   }
   outFile.close();

   // Strip Pitches
   // Indices are 200...209
   outFile.open("strip_pitches.txt");
   for (unsigned int i = 1; i<=24; i++){
      int histIndex = i*1000;
      int pedeIndex = 100000 + i*1000;
      outFile << resultMap[pedeIndex + 200] + hFDCConstants->GetBinContent(histIndex+200) << " " ;
      outFile << resultMap[pedeIndex + 201] + hFDCConstants->GetBinContent(histIndex+201) << " " ;
      outFile << resultMap[pedeIndex + 202] + hFDCConstants->GetBinContent(histIndex+202) << " " ;
      outFile << resultMap[pedeIndex + 203] + hFDCConstants->GetBinContent(histIndex+203) << " " ;
      outFile << resultMap[pedeIndex + 204] + hFDCConstants->GetBinContent(histIndex+204) << " " ;
      outFile << resultMap[pedeIndex + 205] + hFDCConstants->GetBinContent(histIndex+205) << " " ;
      outFile << resultMap[pedeIndex + 206] + hFDCConstants->GetBinContent(histIndex+206) << " " ;
      outFile << resultMap[pedeIndex + 207] + hFDCConstants->GetBinContent(histIndex+207) << " " ;
      outFile << resultMap[pedeIndex + 208] + hFDCConstants->GetBinContent(histIndex+208) << " " ;
      outFile << resultMap[pedeIndex + 209] + hFDCConstants->GetBinContent(histIndex+209) << endl;
   }
   outFile.close();

   // Wire Alignment
   // Indices are 2
   outFile.open("wire_alignment.txt");
   for (unsigned int i = 1; i<=24; i++){
      int histIndex = i*1000;
      int pedeIndex = 100000 + i*1000;
      outFile << resultMap[pedeIndex + 2] + hFDCConstants->GetBinContent(histIndex+2) << " 0.0 0.0" << endl;
   }
   outFile.close();

}
