// $Id$
//
//    File: DDIRCLut.h
//

#ifndef _DDIRCLut_
#define _DDIRCLut_

#include <JANA/JFactory.h>
#include <JANA/JObject.h>
using namespace jana;

#include "TVector3.h"
#include "TFile.h"
#include "TTree.h"

class DDIRCLut : public JObject {

public:
	
	JOBJECT_PUBLIC(DDIRCLut);

	DDIRCLut();
	~DDIRCLut(){}

	uint GetLutPixelAngleSize(int bar, int pixel) const;
	uint GetLutPixelTimeSize(int bar, int pixel) const;
	uint GetLutPixelPathSize(int bar, int pixel) const;
	TVector3 GetLutPixelAngle(int bar, int pixel, int entry) const;
	Double_t GetLutPixelTime(int bar, int pixel, int entry) const;
	Long64_t GetLutPixelPath(int bar, int pixel, int entry) const;
	
private:

	vector<TVector3> lutNodeAngle[48][10864];
	vector<Double_t> lutNodeTime[48][10864];
	vector<Long64_t> lutNodePath[48][10864];
};

#endif // _DDIRCLut_

