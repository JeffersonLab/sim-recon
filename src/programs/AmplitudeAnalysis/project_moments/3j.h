#ifndef THREEJ_H__
#define THREEJ_H__

#include <vector>
#include <map>

#include "wave.h"

class wavePair : public std::pair<size_t, size_t>
{
 public:
  wavePair(size_t a, size_t b)
    : std::pair<size_t,size_t>(a,b)
    {}


  bool operator<(const wavePair& o)
  {
    if (o.first < this->first)
      return true;
    return o.second < this->second;
  }
};

double threeJ(long j1, long j2, long j3, long m1, long m2, long m3);

std::vector<std::pair<size_t, size_t> > listOfMoments(const waveset& ws);

double getCoefficient(int eps, int L, int M, int l1, int m1, int l2, int m2);
double decomposeMoment(const std::pair<size_t, size_t>& LM, const waveset& ws, const vector<double>& x);
double decomposeMoment(const std::pair<size_t, size_t>& LM, const waveset& ws, const double* x);
double decomposeMoment(int L, int M, const waveset& ws, const vector<double>& x);
double decomposeMoment(int L, int M, const waveset& ws, const double* x);

double decomposeMomentError(const std::pair<size_t, size_t>& LM, const waveset& ws, const vector<double>& x, const vector< vector< double > >& covMat);
double decomposeMomentError(int L, int M, const waveset& ws, const vector<double>& x, const vector< vector< double > >& covMat);


#endif
