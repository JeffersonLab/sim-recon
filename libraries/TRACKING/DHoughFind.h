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
		void AddPoint(const DVector2 &point);
		void AddPoint(const double &x, const double &y);
		void AddPoints(const vector<DVector2> &points);
		unsigned int GetNPoints(void){return points.size();}
		void ClearPoints(void);
		void FindIndexes(const DVector2 &pos, int &ix, int &iy);
		double FindBeta(double xlo, double ylo, double widthx, double widthy, DVector2 &pos, DVector2 &step);
		double FindBeta(const DVector2 &a, const DVector2 &b, const DVector2 &c, const DVector2 &d);
		void PrintHist(void);
		
	protected:
		vector<DVector2> points;
		double xmin, xmax, ymin, ymax;
		unsigned int Nbinsx, Nbinsy;
		double bin_widthx, bin_widthy, bin_size;
		unsigned int imax_binx, imax_biny;
		double max_bin_content;
	
		double *hist;
	
	private:

};

#endif // _DHoughFind_

