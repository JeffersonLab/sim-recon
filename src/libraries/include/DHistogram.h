// $Id$
//
//    File: DHistogram.h
// Created: Mon Aug  1 08:56:58 EDT 2011
// Creator: davidl (on Linux ifarm1102 2.6.18-128.7.1.el5 x86_64)
//

#ifndef _DHistogram_
#define _DHistogram_


#include <iostream>
#include <cstddef> // for NULL
#include <cmath>

#include <TH1.h>

#include <JANA/jerror.h>

class DHistogram{
	public:
		inline DHistogram(int Nbins, float lowEdge, float highEdge);
		inline DHistogram(const DHistogram &hsrc);
		inline virtual ~DHistogram();
		
		inline void Add(const DHistogram *h);
		
		inline float Fill(float x, float weight=1.0);
		inline int FindBin(float x) const;
		inline int FindFirstBinAbove(float val, int start_bin=1) const;
		inline int FindFirstNonZeroBin(void) const;
		inline int FindLastNonZeroBin(void) const;
		
		inline int GetNbins(void) const;
		inline float GetBinLowEdge(int bin) const;
		inline float GetBinCenter(int bin) const;
		inline float GetBinContent(int bin) const;
		inline float GetBinWidth(void) const;
		inline float GetLowEdge(void) const;
		inline float GetHighEdge(void) const;
      inline float* GetContentPointer(void) const; ///< CAUTION!
		
		inline float Integral(void) const;
		inline TH1D* MakeTH1D(string name, string title) const;
		
		inline void Reset(void);
		inline void Scale(float scale);
		inline void SetBinContent(int bin, float val);
		
		inline DHistogram& operator=(const DHistogram &hsrc);
		
	protected:
	
	
	private:
		DHistogram(); // prevent default constructor usage
	
		int Nbins;
		float lowEdge;
		float highEdge;
		float binWidth;
		float half_binWidth;
		float bin_zero_center;

		float *content;
};


//---------------------------------
// DHistogram    (Constructor)
//---------------------------------
inline DHistogram::DHistogram(int Nbins, float lowEdge, float highEdge)
{
	
	this->Nbins = Nbins;
	this->lowEdge = lowEdge;
	this->highEdge = highEdge;
	binWidth = (highEdge - lowEdge)/(float)Nbins;
	half_binWidth = binWidth/2.0;
	bin_zero_center = lowEdge - half_binWidth;

	if(Nbins>0 && Nbins<100000000){
		content = new float[Nbins+1];
		
		for(int bin=1; bin<=Nbins; bin++){
			content[bin] = 0.0;
		}
	}else{
		content = NULL;
	}
}

//---------------------------------
// DHistogram    (Constructor)
//---------------------------------
inline DHistogram::DHistogram(const DHistogram &hsrc):Nbins(0),lowEdge(0.0),highEdge(0.0),binWidth(0.0),half_binWidth(0.0),content(NULL)
{
	*this = hsrc;
}

//---------------------------------
// ~DHistogram    (Destructor)
//---------------------------------
inline DHistogram::~DHistogram()
{
	if(content!=NULL) delete[] content;
}

//---------------------------------
// Add
//---------------------------------
inline void DHistogram::Add(const DHistogram *h)
{
	if(h->GetBinWidth()!=binWidth || h->GetNbins()!=Nbins){
		_DBG_<<"Histogram parameters don't match!!"<<endl;
		return;
	}
	
	for(int bin=1; bin<=Nbins; bin++){
		content[bin] += h->GetBinContent(bin);
	}
}

//---------------------------------
// FindBin
//---------------------------------
inline int DHistogram::FindBin(float x) const
{
	int bin = 1 + (int)floor((x-lowEdge)/binWidth);
	
	return bin;
}

//---------------------------------
// FindFirstBinAbove
//---------------------------------
inline int DHistogram::FindFirstBinAbove(float val, int start_bin) const
{
	for(int bin=start_bin; bin<=Nbins; bin++){
		if(content[bin] >= val)return bin;
	}
	
	return -1;
}

//---------------------------------
// FindFirstNonZeroBin
//---------------------------------
int DHistogram::FindFirstNonZeroBin(void) const
{
	for(int bin=1; bin<=Nbins; bin++){
		if(content[bin] != 0.0)return bin;
	}
	
	return -1;
}

//---------------------------------
// FindLastNonZeroBin
//---------------------------------
int DHistogram::FindLastNonZeroBin(void) const
{
	for(int bin=Nbins; bin>=1; bin--){
		if(content[bin] != 0.0)return bin;
	}
	
	return -1;
}

//---------------------------------
// Fill
//---------------------------------
inline float DHistogram::Fill(float x, float weight)
{
	int bin = 1 + (int)floor((x-lowEdge)/binWidth);
	if(bin<1 || bin>Nbins)return 0.0;

	content[bin] += weight;
	
	return content[bin];
}
		
//---------------------------------
// GetNbins
//---------------------------------
inline int DHistogram::GetNbins(void) const
{
	return Nbins;
}

//---------------------------------
// GetBinLowEdge
//---------------------------------
inline float DHistogram::GetBinLowEdge(int bin) const
{
	return lowEdge + (float)(bin-1)*binWidth;
}

//---------------------------------
// GetBinCenter
//---------------------------------
inline float DHistogram::GetBinCenter(int bin) const
{
	if(bin<1 || bin>Nbins)return 0.0;

	return bin_zero_center + binWidth*(float)bin;
}

//---------------------------------
// GetBinContent
//---------------------------------
inline float DHistogram::GetBinContent(int bin) const
{
	if(bin<1 || bin>Nbins)return 0.0;

	return content[bin];
}

//---------------------------------
// GetBinWidth
//---------------------------------
inline float DHistogram::GetBinWidth(void) const
{
	return binWidth;
}

//---------------------------------
// GetLowEdge
//---------------------------------
inline float DHistogram::GetLowEdge(void) const
{
	return lowEdge;
}

//---------------------------------
// GetHighEdge
//---------------------------------
inline float DHistogram::GetHighEdge(void) const
{
	return highEdge;
}

//---------------------------------
// GetContentPointer
//---------------------------------
inline float* DHistogram::GetContentPointer(void) const
{
   return content;
}

//---------------------------------
// Integral
//---------------------------------
inline float DHistogram::Integral(void) const
{
	float I=0.0;

	for(int bin=1; bin<=Nbins; bin++) I += content[bin];

	return I;
}

//---------------------------------
// MakeTH1D
//---------------------------------
inline TH1D* DHistogram::MakeTH1D(string name, string title) const
{
	TH1D *h = new TH1D(name.c_str(), title.c_str(), Nbins, lowEdge, highEdge);

	for(int bin=1; bin<=Nbins; bin++) h->SetBinContent(bin, content[bin]);

	return h;
}

//---------------------------------
// Reset
//---------------------------------
inline void DHistogram::Reset(void)
{
	for(int bin=1; bin<=Nbins; bin++)content[bin] = 0.0;
}

//---------------------------------
// Scale
//---------------------------------
inline void DHistogram::Scale(float scale)
{
	for(int bin=1; bin<=Nbins; bin++)content[bin] *= scale;
}

//---------------------------------
// SetBinContent
//---------------------------------
inline void DHistogram::SetBinContent(int bin, float val)
{
	if(bin<1 || bin>Nbins)return;
	
	content[bin] = val;
}

//---------------------------------
// operator=
//---------------------------------
inline DHistogram& DHistogram::operator=(const DHistogram &hsrc)
{
	// If histogram has different parameters, reallocate
	if( (hsrc.GetNbins()!=Nbins) || (hsrc.GetLowEdge()!=lowEdge) || (hsrc.GetHighEdge()!=highEdge)){
		this->Nbins = hsrc.GetNbins();
		this->lowEdge = hsrc.GetLowEdge();
		this->highEdge = hsrc.GetHighEdge();
		binWidth = (highEdge - lowEdge)/(float)Nbins;
		half_binWidth = binWidth/2.0;
		bin_zero_center = lowEdge - half_binWidth;
		
		if(content!=NULL) delete[] content;

		if((Nbins>0) && (Nbins<100000000) && (hsrc.content!=NULL)){
			content = new float[Nbins+1];
		}else{
			content = NULL;
		}
	}
	
	if(content!=NULL){
		
		float *ysrc  = hsrc.content;
		float *ydest = content;
		
		for(int bin=0; bin<=Nbins; bin++){
			*ydest++ = *ysrc++;
		}
	}

	return *this;
}


#endif // _DHistogram_

