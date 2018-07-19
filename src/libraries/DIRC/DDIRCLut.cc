// $Id$
//
//    File: DDIRCLut.cc
//

#include <cassert>
#include <math.h>
using namespace std;

#include "DDIRCLut.h"

//---------------------------------
// DDIRCLut    (Constructor)
//---------------------------------
DDIRCLut::DDIRCLut() 
{
	// retrieve from LUT from file
        const int luts = 48;
	
	TFile *fLut = new TFile("/group/halld/Users/jrsteven/2018-dirc/lut_all_flat.root");
        TTree *tLut=(TTree*) fLut->Get("lut_dirc_flat");
	//tLut->Print();
	
	vector<Double_t> *LutPixelAngleX[luts];
	vector<Double_t> *LutPixelAngleY[luts];
	vector<Double_t> *LutPixelAngleZ[luts];
	vector<Double_t> *LutPixelTime[luts];
	vector<Long64_t> *LutPixelPath[luts];

	// clear arrays to fill from TTree
	for(int l=0; l<luts; l++){
		LutPixelAngleX[l] = 0;
		LutPixelAngleY[l] = 0;
		LutPixelAngleZ[l] = 0;
		LutPixelTime[l] = 0;
		LutPixelPath[l] = 0;
	}

        for(int l=0; l<luts; l++){
		tLut->SetBranchAddress(Form("LUT_AngleX_%d",l),&LutPixelAngleX[l]); 
		tLut->SetBranchAddress(Form("LUT_AngleY_%d",l),&LutPixelAngleY[l]); 
		tLut->SetBranchAddress(Form("LUT_AngleZ_%d",l),&LutPixelAngleZ[l]); 
		tLut->SetBranchAddress(Form("LUT_Time_%d",l),&LutPixelTime[l]); 
		tLut->SetBranchAddress(Form("LUT_Path_%d",l),&LutPixelPath[l]); 
        }

        // fill nodes with LUT info for each bar/pixel combination
	for(int i=0; i<tLut->GetEntries(); i++) { // get pixels from TTree
		tLut->GetEntry(i);
		//if(i<6800 || i>7000) continue;
		//cout<<i<<endl;
		for(int l=0; l<luts; l++){ // loop over bars
			//if(l==0) cout<<"size="<<LutPixelAngleX[l]->size()<<endl;
			for(uint j=0; j<LutPixelAngleX[l]->size(); j++) { // loop over possible paths
				//if(l==31 && i==4432) 
				//	cout<<j<<" "<<LutPixelTime[l]->at(j)<<" "<<LutPixelPath[l]->at(j)<<endl;
				//cout<<"i="<<i<<" l="<<l<<" j="<<j<<endl;
				TVector3 angle(LutPixelAngleX[l]->at(j), LutPixelAngleY[l]->at(j), LutPixelAngleZ[l]->at(j));
				lutNodeAngle[l][i].push_back(angle);
				lutNodeTime[l][i].push_back(LutPixelTime[l]->at(j));
				lutNodePath[l][i].push_back(LutPixelPath[l]->at(j));
			}
		}
	}
	
	//delete LutPixelAngleX;
	//delete LutPixelAngleY;
	//delete LutPixelAngleZ;
	//delete LutPixelTime;
	//delete LutPixelPath;
}

uint DDIRCLut::GetLutPixelAngleSize(int bar, int pixel) const
{
	return lutNodeAngle[bar][pixel].size();
}
	
uint DDIRCLut::GetLutPixelTimeSize(int bar, int pixel) const
{
	return lutNodeTime[bar][pixel].size();
}
	
uint DDIRCLut::GetLutPixelPathSize(int bar, int pixel) const
{
	return lutNodePath[bar][pixel].size();
}

TVector3 DDIRCLut::GetLutPixelAngle(int bar, int pixel, int entry) const
{
	return lutNodeAngle[bar][pixel].at(entry);
}

Double_t DDIRCLut::GetLutPixelTime(int bar, int pixel, int entry) const
{
	return lutNodeTime[bar][pixel].at(entry);
}

Long64_t DDIRCLut::GetLutPixelPath(int bar, int pixel, int entry) const
{
	return lutNodePath[bar][pixel].at(entry);
}
