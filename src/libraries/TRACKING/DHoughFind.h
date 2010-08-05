// $Id$
//
//    File: DHoughFind.h
// Created: Wed Oct 31 05:59:26 EDT 2007
// Creator: davidl (on Darwin Amelia.local 8.10.1 i386)
//

#ifndef _DHoughFind_
#define _DHoughFind_

#include <vector>
using std::vector;

#include <JANA/jerror.h>
#include <JANA/JObject.h>
using namespace jana;

#include <DVector2.h>

class DHoughFind:public JObject{
	public:
		DHoughFind();
		DHoughFind(double xmin, double xmax, double ymin, double ymax, unsigned int Nbinsx, unsigned int Nbinsy);
		virtual ~DHoughFind();
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DHoughFind";}
		
		inline double signof(double a){return a<0.0 ? -1.0:1.0;}
		
		void SetLimits(double xmin, double xmax, double ymin, double ymax, unsigned int Nbinsx, unsigned int Nbinsy);
		void ResetHist(void);
		DVector2 GetMaxBinLocation(void);
		double GetMaxBinContent(void);
		double GetSigmaX(void);
		double GetSigmaY(void);
		DVector2 Find(void);
		DVector2 Find(const vector<DVector2> &points);
		
		void Fill(double x, double sigmax, double y, double sigmay);
		static DVector2 GetMaxBinLocation(vector<const DHoughFind*> &houghs); // does not look at "this" object!
		
		void AddPoint(const DVector2 &point);
		void AddPoint(const double &x, const double &y);
		void AddPoints(const vector<DVector2> &points);
		unsigned int GetNPoints(void){return points.size();}
		void ClearPoints(void);
		void PrintHist(void);
		
		inline void FindIndexes(const DVector2 &pos, int &ix, int &iy);
		inline double FindBeta(double xlo, double ylo, double widthx, double widthy, DVector2 &pos, DVector2 &step);
		inline double FindBeta(const DVector2 &a, const DVector2 &b, const DVector2 &c, const DVector2 &d);

	protected:
		vector<DVector2> points;
		double xmin, xmax, ymin, ymax;
		unsigned int Nbinsx, Nbinsy;
		double bin_widthx, bin_widthy, bin_size;
		bool max_bin_valid;
		unsigned int imax_binx, imax_biny;
		double max_bin_content;
	
		double *hist;

		DVector2 a0;
		DVector2 xdir;
		DVector2 ydir;
		DVector2 stepdir;
		DVector2 start;
	
	private:

};


// The following functions are inlined for speed

//---------------------------------
// FindIndexes
//---------------------------------
inline void DHoughFind::FindIndexes(const DVector2 &pos, int &ix, int &iy)
{
	ix = (int)floor((pos.X()-xmin)/bin_widthx);
	iy = (int)floor((pos.Y()-ymin)/bin_widthy);
}

//---------------------------------
// FindBeta
//---------------------------------
inline double DHoughFind::FindBeta(double xlo, double ylo, double widthx, double widthy, DVector2 &pos, DVector2 &step)
{
	//DVector2 a0(xlo, ylo);
	//DVector2 xdir(1.0, 0.0);
	//DVector2 ydir(0.0, 1.0);
	//DVector2 stepdir = step/step.Mod();
	a0.Set(xlo, ylo);
	xdir.Set(1.0, 0.0);
	ydir.Set(0.0, 1.0);
	stepdir = step/step.Mod();

	//vector<double> beta(4);
	double beta[4];
	//start = a0;
	start.Set(xlo, ylo);
	beta[0] = FindBeta(start, xdir, pos, stepdir);
	//start = a0+widthx*xdir;
	start.Set(xlo+widthx, ylo);
	beta[1] = FindBeta(start, ydir, pos, stepdir);
	//start = a0+widthx*xdir+widthy*ydir;
	start.Set(xlo+widthx, ylo+widthy);
	beta[2] = FindBeta(start, -1.0*xdir, pos, stepdir);
	//start = a0+widthy*ydir;
	start.Set(xlo, ylo+widthy);
	beta[3] = FindBeta(start, -1.0*ydir, pos, stepdir);

	// Here, we want to choose the closest boundary of the bin which is not
	// the boundary the point "pos" lies on. The values of beta give the
	// distances to the lines which define the boundaries, but not clipped
	// to the bin itself.

	// Loop over beta values
	DVector2 bin_center(xlo+widthx/2.0, ylo+widthy/2.0);
	double min_dist_to_center = 1.0E6;
	double min_beta = 1.0E6;
	for(unsigned int i=0; i<4; i++){
		if(fabs(beta[i])<1.0E-6*bin_size)continue; // This is most likely due to our starting pos being on a boundary already. Ignore it.
		
		DVector2 delta_center = pos + beta[i]*stepdir - bin_center;
		if(delta_center.Mod() < min_dist_to_center){
			min_dist_to_center = delta_center.Mod();
			min_beta = beta[i];
		}
	}
	
	return min_beta;
}

//---------------------------------
// FindBeta
//---------------------------------
inline double DHoughFind::FindBeta(const DVector2 &a, const DVector2 &b, const DVector2 &c, const DVector2 &d)
{
	/// Given the 2-D vectors a,b,c, and d find the scaler value "beta" such that
	///
	/// a + alpha*b = c + beta*d
	///
	/// It vectors b and d should be unit vectors.
	/// The value of alpha is not returned, but can be calculated via:
	///
	/// alpha = beta*(d.b) - b.(a-c)
	///
	/// If the vectors b and d are parallel, then the return value will be inf
	/// or whatever division by zero returns on the system.

	double ax = a.X();
	double ay = a.Y();
	double bx = b.X();
	double by = b.Y();
	double cx = c.X();
	double cy = c.Y();
	double dx = d.X();
	double dy = d.Y();
	
	double k = bx*dx + by*dy;
	double alphax = ax - cx;
	double alphay = ay - cy;
	double betax = dx - k*bx;
	double betay = dy - k*by;
	
	return (betax*alphax + betay*alphay)/(1.0-k*k);

	//return (d*(a-c) - (d*b)*(b*(a-c)))/(1.0-pow(b*d, 2.0));
}

#endif // _DHoughFind_

