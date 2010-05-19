#include <math.h>
#include <emmintrin.h>
#include <iostream>
using namespace std;

class DMatrix2x1{
  public:
  DMatrix2x1(){
    mA.v=_mm_setzero_pd();
  };
  DMatrix2x1(const double x, const double y){
    mA.v=_mm_setr_pd(x,y);
  }
  DMatrix2x1(__m128d v){
    mA.v=v;
  }
  ~DMatrix2x1(){};

  __m128d GetV() const{
    return mA.v;
  }

  // Set the two components of the matrix
  void Set(const double c1,const double c2){
    mA.v=_mm_setr_pd(c1,c2);
  }

  // Access by row 
  double &operator() (int row){
    return mA.d[row];
  } 
  double operator() (int row) const{
    return mA.d[row];
  } 

  // Matrix subtraction
  DMatrix2x1 operator-(const DMatrix2x1 &m2) const{
    return DMatrix2x1(_mm_sub_pd(GetV(),m2.GetV()));
  }

  
  void Print(){
    cout << "DMatrix2x1:" <<endl;
    cout << "     |      0    |" <<endl;
    cout << "----------------------" <<endl;
    for (unsigned int i=0;i<2;i++){
      cout <<"   "<<i<<" |" <<  setw(11)<<setprecision(4) << mA.d[i] << endl;
    }      
  }

  

 private:
  union dvec{
    __m128d v;
    double d[2];
  }mA;    
};

// Scale 2x1 matrix by a floating point number
inline DMatrix2x1 operator*(const double c,const DMatrix2x1 &M){
  return DMatrix2x1(_mm_mul_pd(M.GetV(),_mm_set1_pd(c)));
}
