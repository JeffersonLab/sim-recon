// Author: David Lawrence  Feb. 4, 2005
//
//
// DArcHit.cc

#include <iostream>
#include <iomanip>
using namespace std;

#include "DArcHit.h"


//--------------
// DArcHit
//--------------
DArcHit::DArcHit(float x, float y, float z)
{
	SetXYZ(x,y,z);
}

//--------------
// SetXYZ
//--------------
derror_t DArcHit::SetXYZ(float x, float y, float z)
{
	this->xhit = x;
	this->yhit = y;
	this->zhit = z;
	this->rhit = sqrt(x*x + y*y);
	this->phihit = atan2(y, x);
	this->used = 0;
	 
	if(y>x){
		m = -x/y;
		b = -1.0/(2.0*y)*(x*x + y*y);
		orientation = Y_OF_X;
	}else{
		m = -y/x;
		b = -1.0/(2.0*x)*(x*x + y*y);
		orientation = X_OF_Y;
	}
	
	return NOERROR;
}
 
//--------------
// Dist2ToLine
//--------------
float DArcHit::Dist2ToLine(float x, float y)
{
	/// x and y represent a point (the focus as identified
	/// by a peak on the density plot).
	/// This routine returns the distance from the x,y point to
	/// the line on which the circle centers must lie as defined by
	/// the hit.
	float x1,y1;
	if(orientation == Y_OF_X){		
		x1 = (x-m*(b-y))/(1.0+m*m);
		y1 = m*x1 + b;
	}else{
		y1 = (y-m*(b-x))/(1.0+m*m);
		x1 = m*y1 + b;
	}

	//float d2 = pow((double)(x-x1),2.0) + pow((double)(y-y1),2.0);
	float dx = x-x1;
	float dy = y-y1;
	
	return dx*dx + dy*dy;
}

//--------------
// DistToLine
//--------------
float DArcHit::DistToLine(float x, float y)
{
	return sqrt(Dist2ToLine(x,y));
}
//--------------
// Density
//--------------
float DArcHit::Density(float x, float y)
{
	float d2 = Dist2ToLine(x,y);

	//float d = sqrt(d2);
	//float frac = (5.0 - d)/5.0;
	float frac = (25.0 - d2)/25.0;
	return frac<0.0 ? 0.0:frac;

	//return exp(-d2*d2/1.0);
}
 
//--------------
// FillDensityHistogram
//--------------
derror_t DArcHit::FillArcDensityHistogram(TH2F *density)
{
	/// Fill the specified density histogram using the
	/// line-of-circle-centers determined by this hit.
	
	// Get useful info from histogram
	TAxis *xaxis = density->GetXaxis();
	TAxis *yaxis = density->GetYaxis();
	int NbinsX = xaxis->GetNbins();
	int NbinsY = yaxis->GetNbins();
	float xmin = xaxis->GetBinCenter(1);
	float xmax = xaxis->GetBinCenter(NbinsX);
	float ymin = yaxis->GetBinCenter(1);
	float ymax = yaxis->GetBinCenter(NbinsY);

	// Things will go much quicker if we just fill in the bins
	// surrounding the line. To find these though, we need to 
	// first determine the two points where the line crosses
	// the boundaries of the histogram.
	// Here, the point the line intersects with each of the
	// 4 limits (top, bottom, left, and right) is found. In
	// the end, only zero or two of these intersection points
	// will be contained in the limits of the histogram.
	// (we actually use the boundary at 1/2 a bin in from the
	// edges of the histo.)
	float x[4], y[4];
	if(orientation==Y_OF_X){
		// left
		x[0] = xmin;
		y[0] = m*x[0]+b;
		// right
		x[1] = xmax;
		y[1] = m*x[1]+b;
		//top
		y[2] = ymax;
		x[2] = (y[2] - b)/m;
		if(!finite(x[2]))x[2] = 2.0*xmax; // for m=0 case
		//bottom
		y[3] = ymin;
		x[3] = (y[3] - b)/m;
		if(!finite(x[3]))x[3] = 2.0*xmax; // for m=0 case
	}else{
		// top
		y[0] = ymax;
		x[0] = m*y[0]+b;
		// bottom
		y[1] = ymin;
		x[1] = m*y[1]+b;
		// left
		x[2] = xmin;
		y[2] = (x[2] - b)/m;
		if(!finite(y[2]))y[2] = 2.0*ymax; // for m=0 case
		// right
		x[3] = xmax;
		y[3] = (x[3] - b)/m;
		if(!finite(y[3]))y[3] = 2.0*ymax; // for m=0 case
	}

	// Count and remember the border crossings which occur inside the histo limits
	int goodpoints[4]; // index of x[] and y[] for points contained in histo
	int Ngoodpoints = 0;
	for(int i=0;i<4;i++){
		if(x[i]<xmin || x[i]>xmax)continue;
		if(y[i]<ymin || y[i]>ymax)continue;
		goodpoints[Ngoodpoints++] = i;
	}
	
	switch(Ngoodpoints){
		case 0:
			cout<<__FILE__<<":"<<__LINE__<<" Line does not enter histogram!"<<endl;
			return NOERROR;
			break;
		case 1:
		case 3:
			cout<<__FILE__<<":"<<__LINE__<<" Odd number of border crossings. This should NEVER happen!"<<endl;
			break;
		case 2:
			// This is what SHOULD happen.
			break;
	}
	
	// Now we basically want to loop over bins which span across the
	// line for every row or column. The number of bins across the
	// line to span is dependant on the "width" of the line (i.e.
	// what we use as the "density as a function of distance" function)
	// The key thing here is that the inner loop will loop over a fixed
	// and small number of bins so it needs to be in the direction
	// most perpendicular to the line (either x or y).
	float spanwidth = 4.0; // in cm
	if(orientation==Y_OF_X){
		// line has shallow slope.
		// histogram will be filled for spanbins bins on either
		// side of each bin the line cross through.
		int spanbins = 1+ (int)(spanwidth/yaxis->GetBinWidth(1));
		int xbin1 = xaxis->FindBin(x[goodpoints[0]]);
		int xbin2 = xaxis->FindBin(x[goodpoints[1]]);
		for(int xbin=xbin1 ; xbin!=xbin2 ; xbin+=xbin1<xbin2 ? 1:-1){
			float xcenter = xaxis->GetBinCenter(xbin);
			int ybincenter = yaxis->FindBin(m*xaxis->GetBinCenter(xbin) + b);
			for(int ybin=ybincenter-spanbins ; ybin<=ybincenter+spanbins ; ybin++){
				if(ybin<1 || ybin>NbinsY)continue;
				float ycenter = yaxis->GetBinCenter(ybin);
				density->Fill(xcenter, ycenter, Density(xcenter,ycenter));
			}
		}
		
	}else{
		// line has steep slope.
		// histogram will be filled for spanbins bins on either
		// side of each bin the line cross through.
		int spanbins = 1+ (int)(spanwidth/xaxis->GetBinWidth(1));
		int ybin1 = yaxis->FindBin(y[goodpoints[0]]);
		int ybin2 = yaxis->FindBin(y[goodpoints[1]]);
		for(int ybin=ybin1 ; ybin!=ybin2 ; ybin+=ybin1<ybin2 ? 1:-1){
			float ycenter = yaxis->GetBinCenter(ybin);
			int xbincenter = xaxis->FindBin(m*yaxis->GetBinCenter(ybin) + b);
			for(int xbin=xbincenter-spanbins ; xbin<=xbincenter+spanbins ; xbin++){
				if(xbin<1 || xbin>NbinsX)continue;
				float xcenter = xaxis->GetBinCenter(xbin);
				density->Fill(xcenter, ycenter, Density(xcenter,ycenter));
			}
		}
	}

	return NOERROR;
}
 
 
