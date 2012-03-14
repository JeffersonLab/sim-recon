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
		inline DHistogram(int Nbins, double lowEdge, double highEdge);
		inline DHistogram(DHistogram &hsrc);
		inline virtual ~DHistogram();
		
		inline void Add(const DHistogram *h);
		
		inline double Fill(double x, double weight=1.0);
		inline int FindBin(double x) const;
		inline int FindFirstBinAbove(double val, int start_bin=1) const;
		inline int FindFirstNonZeroBin(void) const;
		inline int FindLastNonZeroBin(void) const;
		
		inline int GetNbins(void) const;
		inline double GetBinLowEdge(int bin) const;
		inline double GetBinCenter(int bin) const;
		inline double GetBinContent(int bin) const;
		inline double GetBinWidth(void) const;
		inline double GetLowEdge(void) const;
		inline double GetHighEdge(void) const;
      inline double* GetContentPointer(void) const; ///< CAUTION!
      inline double* GetXValsPointer(void) const;   ///< CAUTION!
		
		inline double Integral(void) const;
		inline TH1D* MakeTH1D(string name, string title) const;
		
		inline void Reset(void);
		inline void Scale(double scale);
		inline void SetBinContent(int bin, double val);
		
		inline DHistogram& operator=(const DHistogram &hsrc);
		
	protected:
	
	
	private:
		DHistogram(); // prevent default constructor usage
	
		int Nbins;
		double lowEdge;
		double highEdge;
		double binWidth;
		double half_binWidth;

		double *content;
		double *xvals; // bin center
};


//---------------------------------
// DHistogram    (Constructor)
//---------------------------------
inline DHistogram::DHistogram(int Nbins, double lowEdge, double highEdge)
{
	
	this->Nbins = Nbins;
	this->lowEdge = lowEdge;
	this->highEdge = highEdge;
	binWidth = (highEdge - lowEdge)/(double)Nbins;
	half_binWidth = binWidth/2.0;

	if(Nbins>0 && Nbins<100000000){
		content = new double[Nbins+1];
		xvals = new double[Nbins+1];
		
		for(int bin=1; bin<=Nbins; bin++){
			content[bin] = 0.0;
			xvals[bin] = lowEdge + (double)(bin-1)*binWidth + half_binWidth;
		}
	}else{
		content = NULL;
		xvals = NULL;
	}
}

//---------------------------------
// DHistogram    (Constructor)
//---------------------------------
inline DHistogram::DHistogram(DHistogram &hsrc):Nbins(0),lowEdge(0.0),highEdge(0.0),binWidth(0.0),half_binWidth(0.0),content(NULL),xvals(NULL)
{
	*this = hsrc;
}

//---------------------------------
// ~DHistogram    (Destructor)
//---------------------------------
inline DHistogram::~DHistogram()
{
	if(content!=NULL) delete[] content;
	if(xvals!=NULL) delete[] xvals;
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
inline int DHistogram::FindBin(double x) const
{
	int bin = 1 + (int)floor((x-lowEdge)/binWidth);
	
	return bin;
}

//---------------------------------
// FindFirstBinAbove
//---------------------------------
inline int DHistogram::FindFirstBinAbove(double val, int start_bin) const
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
inline double DHistogram::Fill(double x, double weight)
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
inline double DHistogram::GetBinLowEdge(int bin) const
{
	return lowEdge + (double)(bin-1)*binWidth;
}

//---------------------------------
// GetBinCenter
//---------------------------------
inline double DHistogram::GetBinCenter(int bin) const
{
	if(bin<1 || bin>Nbins)return 0.0;

	return xvals[bin];
}

//---------------------------------
// GetBinContent
//---------------------------------
inline double DHistogram::GetBinContent(int bin) const
{
	if(bin<1 || bin>Nbins)return 0.0;

	return content[bin];
}

//---------------------------------
// GetBinWidth
//---------------------------------
inline double DHistogram::GetBinWidth(void) const
{
	return binWidth;
}

//---------------------------------
// GetLowEdge
//---------------------------------
inline double DHistogram::GetLowEdge(void) const
{
	return lowEdge;
}

//---------------------------------
// GetHighEdge
//---------------------------------
inline double DHistogram::GetHighEdge(void) const
{
	return highEdge;
}

//---------------------------------
// GetContentPointer
//---------------------------------
inline double* DHistogram::GetContentPointer(void) const
{
   return content;
}

//---------------------------------
// GetXValsPointer
//---------------------------------
inline double* DHistogram::GetXValsPointer(void) const
{
   return xvals;
}

//---------------------------------
// Integral
//---------------------------------
inline double DHistogram::Integral(void) const
{
	double I=0.0;

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
inline void DHistogram::Scale(double scale)
{
	for(int bin=1; bin<=Nbins; bin++)content[bin] *= scale;
}

//---------------------------------
// SetBinContent
//---------------------------------
inline void DHistogram::SetBinContent(int bin, double val)
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
		binWidth = (highEdge - lowEdge)/(double)Nbins;
		half_binWidth = binWidth/2.0;
		
		if(content!=NULL) delete[] content;
		if(xvals!=NULL) delete[] xvals;

		if((Nbins>0) && (Nbins<100000000) && (hsrc.content!=NULL) && (hsrc.xvals!=NULL)){
			content = new double[Nbins+1];
			xvals = new double[Nbins+1];
		}else{
			content = NULL;
			xvals = NULL;
		}
	}
	
	if(content!=NULL || xvals!=NULL){
		
		double *ysrc  = hsrc.content;
		double *xsrc  = hsrc.xvals;
		double *ydest = content;
		double *xdest = xvals;
		
		for(int bin=0; bin<=Nbins; bin++){
			*ydest++ = *ysrc++;
			*xdest++ = *xsrc++;
		}
	}

	return *this;
}


#endif // _DHistogram_

