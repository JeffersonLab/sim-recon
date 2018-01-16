#include <iostream>
#include <math.h>
#include <string.h>
#include <stdlib.h>

using namespace std;

#include "Math/SpecFunc.h"

#include "wave.h"
#include "3j.h"

namespace {
  #if ROOT_VERSION_CODE < ROOT_VERSION(5,28,0)
  long long
  fac(long long n)
  {
    long long result = 1;
    if (n > 1)
      for (int i = 1; i <= n; i++)
	result *= i;
    return result;
  }


  // Calculate the 3j-symbol from the explicit expression given in http://dlmf.nist.gov/34.2
  // First their delta (34.2.5):

  double
  delta(long j1, long j2, long j3)
  {
    // Valid argument range is assumed.  This is assured in the 3j function.
    double result = (sqrt(fac(j1 + j2 - j3) * fac(j1 - j2 + j3) * fac(-j1 + j2 + j3))
		     / sqrt(fac(j1 + j2 + j3 + 1)));
    //cout << "delta = " << result << endl;
    return result;
  }


  // The complete prefactor in (34.2.4).  Again, rangechecks are left to the main function.
  double
  prefactor(long j1, long j2, long j3, long m1, long m2, long m3)
  {
    int sign = ((j1 - j2 - m3) & 0x1) ? -1 : 1;  // odd / even
    double result = sign*delta(j1, j2, j3)*sqrt(fac(j1+m1)*fac(j1-m1)*fac(j2+m2)*fac(j2-m2)*fac(j3+m3)*fac(j3-m3));
    //cout << "prefactor = " << result << endl;
    return result;
  }
  #endif


  double
  theta(int m)
  {
    if (m == 0)
      return .5;
    else
      return sqrt(.5);
  }
}

double
getCoefficient(int eps, int L, int M, int l1, int m1, int l2, int m2)
{
  double threeJ1 = threeJ(L, l1, l2, 0, 0, 0);

  int sign = ((m1 + M) & 0x1) ? -1 : 1;
  int sign1 = ((m1 + 1) & 0x1) ? -1 : 1;
  int sign2 = ((m2 + 1) & 0x1) ? -1 : 1;
  int sign12 = ((m1 + m2) & 0x1) ? -1 : 1;
  double parentheses = (threeJ(L, l1, l2, -M, -m1, m2)
			+ threeJ(L, l1, l2, -M, m1, m2) * eps * sign1
			+ threeJ(L, l1, l2, -M, -m1, -m2) * eps * sign2
			+ threeJ(L, l1, l2, -M, m1, -m2) * sign12) * sign;

  return(theta(m1)*theta(m2)*sqrt((2*l1+1)*(2*l2+1))
	 *threeJ1*parentheses);
}


double
threeJ(long j1, long j2, long j3, long m1, long m2, long m3)
{
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,28,0)
  return ROOT::Math::wigner_3j(2*j1, 2*j2, 2*j3, 2*m1, 2*m2, 2*m3);
#else
  if (j1 < 0 || j2 < 0 || j3 < 0)
    return 0;
  if (j1 > j2 + j3)
    return 0;
  if (j2 > j3 + j1)
    return 0;
  if (j3 > j1 + j2)
    return 0;
  if (abs(j1 - j2) > j3)
    return 0;
  if (abs(j2 - j3) > j1)
    return 0;
  if (abs(j3 - j1) > j2)
    return 0;
  if (m1 + m2 + m3 != 0)
    return 0;
  if (abs(m1) > j1 || abs(m2) > j2 || abs(m3) > j3)
    return 0;

  // The 3j-symbol is != 0.  Find the valid range for the summation in loc.cit. (34.2.4)
  long minS = 0;
  minS = max(minS, -(j3 - j2 + m1));
  minS = max(minS, -(j3 - j1 - m2));
  long maxS = j1 + j2 - j3;  // needs refinement
  maxS = min(maxS, j1 - m1);
  maxS = min(maxS, j2 + m2);

  if (minS > maxS)
    return 0;

  double sum = 0;
  for (int s = minS; s <= maxS; s++)
    {
      double add = ((s & 0x1 ? -1. : 1.)
		    / fac(s) / fac(j1+j2-j3-s)
		    / fac(j1-m1-s) / fac(j2+m2-s)
		    / fac(j3-j2+m1+s) / fac(j3-j1-m2+s));
      //cout << add << endl;
      sum+=add;
    }

  return prefactor(j1,j2,j3,m1,m2,m3)*sum;
#endif
}


// Returns the list of non-zero moments for the given waveset.
std::vector<std::pair<size_t, size_t> >
listOfMoments(const waveset& ws)
{
  // Find maximum L, M.
  size_t maxL = 0, maxM = 0;
  for (size_t iWs = 0; iWs < ws.size(); iWs++)
    {
      const vector<wave>& w = ws[iWs].waves;
      for (size_t iW = 0; iW < w.size(); iW++)
	{
	  maxL = std::max(w[iW].l, maxL);
	  maxM = std::max(w[iW].m, maxM);
	}
    }

  std::vector<std::pair<size_t, size_t> > result;
  for (size_t L = 0; L <= 2*maxL; L++)
    for (size_t M = 0; M <= 2*maxM; M++)
      {
	for (size_t iWs = 0; iWs < ws.size(); iWs++)
	  {
	    int eps = ws[iWs].reflectivity;

	    const vector<wave>& w = ws[iWs].waves;
	    for (size_t iW1 = 0; iW1 < w.size(); iW1++)
	      {
		const wave& w1 = w[iW1];
		for (size_t iW2 = 0; iW2 < w.size(); iW2++)
		  {
		    const wave& w2 = w[iW2];
		    double coeff = getCoefficient(eps, L, M, w1.l, w1.m, w2.l, w2.m);
		    if (coeff != 0)
		      {
			result.push_back(std::pair<size_t,size_t>(L, M));
			goto nextM;
		      }
		  }
	      }
	  }
      nextM: ;
      }
  return result;
}


double
decomposeMoment(const std::pair<size_t, size_t>& LM, const waveset& ws, const vector<double>& x)
{
  return decomposeMoment(LM.first, LM.second, ws, x);
}

double
decomposeMoment(const std::pair<size_t, size_t>& LM, const waveset& ws, const double* x)
{
  return decomposeMoment(LM.first, LM.second, ws, x);
}

// Calculates moment H(LM) from the waveset ws with corresponding fit results x.
double
decomposeMoment(int L, int M, const waveset& ws, const vector<double>& x)
{
  return decomposeMoment(L, M, ws, &x[0]);
}

double
decomposeMoment(int L, int M, const waveset& ws, const double* x)
{
  double result = 0;
  for (size_t iWs = 0; iWs < ws.size(); iWs++)
    {
      int eps = ws[iWs].reflectivity;

      const vector<wave>& w = ws[iWs].waves;
      for (size_t iW1 = 0; iW1 < w.size(); iW1++)
	{
	  const wave& w1 = w[iW1];
	  for (size_t iW2 = 0; iW2 < w.size(); iW2++)
	    {
	      const wave& w2 = w[iW2];
	      double coeff = getCoefficient(eps, L, M, w1.l, w1.m, w2.l, w2.m);
	      result += coeff*(x[w1.getIndex()]*x[w2.getIndex()] + x[w1.getIndex()+1]*x[w2.getIndex()+1]);
	      //result[wavePair(w1.getIndex(), w2.getIndex())] = coeff;
	    }
	}
    }

  return result;
}

double
decomposeMomentError(const std::pair<size_t, size_t>& LM,
		     const waveset& ws, const vector<double>& x, const vector< vector< double > >& covMat)
{
  return decomposeMomentError(LM.first, LM.second, ws, x, covMat);
}

double
decomposeMomentError(int L, int M, const waveset& ws, const vector<double>& x, const vector< vector< double > >& covMat)
{
  double resultSquare = 0;
  for (size_t iWs = 0; iWs < ws.size(); iWs++)
    {
      int eps = ws[iWs].reflectivity;

      const vector<wave>& w = ws[iWs].waves;
      for (size_t iW1 = 0; iW1 < w.size(); iW1++)
  	{
  	  const wave& w1 = w[iW1];
  	  for (size_t iW2 = 0; iW2 < w.size(); iW2++)
  	    {
  	      const wave& w2 = w[iW2];
  	      double coeff = getCoefficient(eps, L, M, w1.l, w1.m, w2.l, w2.m);

  	      double re1 = x[w1.getIndex()];
  	      double im1 = x[w1.getIndex() + 1];
  	      double re2 = x[w2.getIndex()];
  	      double im2 = x[w2.getIndex() + 1];

  	      double errRe1_2 = covMat[w1.getIndex()][w1.getIndex()];
  	      double errIm1_2 = covMat[w1.getIndex() + 1][w1.getIndex() + 1];
  	      double errRe2_2 = covMat[w2.getIndex()][w2.getIndex()];
  	      double errIm2_2 = covMat[w2.getIndex() + 1][w2.getIndex() + 1];

  	      double covReRe = covMat[w1.getIndex()][w2.getIndex()];

  	      // The covariance between the imaginary parts is only non-zero if neither was fixed
  	      double covImIm = covMat[w1.getIndex() + 1][w2.getIndex() + 1];

  	      resultSquare += fabs(coeff)*(re2*re2*errRe1_2
  					   + re1*re1*errRe2_2
  					   + im2*im2*errIm1_2
  					   + im1*im1*errIm2_2
  					   + 2*re1*re2*covReRe
  					   + 2*im1*im2*covImIm);
  	    }
  	}
    }

  return sqrt(resultSquare);
}
