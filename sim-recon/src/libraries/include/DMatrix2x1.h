#include <align_16.h>

#ifndef USE_SSE2 

// Matrix class without SIMD instructions

class DMatrix2x1{
public:
  DMatrix2x1(){
    mA[0]=mA[1]=0.;
  }
  DMatrix2x1(double x, double y){
    mA[0]=x;
    mA[1]=y;
  }
  ~DMatrix2x1(){};

  // Set the two components of the matrix
  void Set(double c1, double c2){
    mA[0]=c1;
    mA[1]=c2;
  }
  // Access by row 
  double &operator() (int row){
    return mA[row];
  } 
  double operator() (int row) const{
    return mA[row];
  } 

  // Matrix subtraction
  DMatrix2x1 operator-(const DMatrix2x1 &m2) const{
    return DMatrix2x1(mA[0]-m2(0),mA[1]-m2(1));
  } 

  // Matrix addition
  DMatrix2x1 operator+(const DMatrix2x1 &m2) const{
    return DMatrix2x1(mA[0]+m2(0),mA[1]+m2(1));
  }

  void Print(){
     cout << "DMatrix2x1:" <<endl;
     cout << "     |      0    |" <<endl;
    cout << "----------------------" <<endl;
    for (unsigned int i=0;i<2;i++){
      cout <<"   "<<i<<" |" <<  setw(11)<<setprecision(4) << mA[i] << endl;
    }      
  }


private:
  double mA[2];

};


// Scale 2x1 matrix by a floating point number
inline DMatrix2x1 operator*(const double c,const DMatrix2x1 &M){ 
  return DMatrix2x1(c*M(0),c*M(1));
}




#else

// Matrix class with SIMD instructions

class DMatrix2x1{
  public:
  DMatrix2x1()
  : mA( ALIGNED_16_BLOCK_PTR(union dvec, 1, mA) )
  {
    mA->v=_mm_setzero_pd();
  };
  DMatrix2x1(double x, double y)
  : mA( ALIGNED_16_BLOCK_PTR(union dvec, 1, mA) )
  {
    mA->v=_mm_setr_pd(x,y);
  }
  DMatrix2x1(__m128d v)
  : mA( ALIGNED_16_BLOCK_PTR(union dvec, 1, mA) )
  {
    mA->v=v;
  }
  DMatrix2x1(const DMatrix2x1& dm)
  : mA( ALIGNED_16_BLOCK_PTR(union dvec, 1, mA) )
  {
    mA->v=dm.mA->v;
  }
  ~DMatrix2x1(){};

  __m128d GetV() const{
    return mA->v;
  }

  // Assignment
  DMatrix2x1& operator=(const DMatrix2x1& dm){
    mA->v=dm.mA->v;
    return *this;
  }
          
  // Set the two components of the matrix
  void Set(double c1, double c2){
    mA->v=_mm_setr_pd(c1,c2);
  }

  // Access by row 
  double &operator() (int row){
    return mA->d[row];
  } 
  double operator() (int row) const{
    return mA->d[row];
  } 

  // Matrix subtraction
  DMatrix2x1 operator-(const DMatrix2x1 &m2) const{
    return DMatrix2x1(_mm_sub_pd(GetV(),m2.GetV()));
  }

  // Matrix addition
  DMatrix2x1 operator+(const DMatrix2x1 &m2) const{
    return DMatrix2x1(_mm_add_pd(GetV(),m2.GetV()));
  }
  
  void Print(){
    cout << "DMatrix2x1:" <<endl;
    cout << "     |      0    |" <<endl;
    cout << "----------------------" <<endl;
    for (unsigned int i=0;i<2;i++){
      cout <<"   "<<i<<" |" <<  setw(11)<<setprecision(4) << mA->d[i] << endl;
    }      
  }

  

 private:
  union dvec{
    __m128d v;
    double d[2];
  };
  ALIGNED_16_BLOCK(union dvec, 1, mA)
};

// Scale 2x1 matrix by a floating point number
inline DMatrix2x1 operator*(double c, const DMatrix2x1 &M){
  return DMatrix2x1(_mm_mul_pd(M.GetV(),_mm_set1_pd(c)));
}
#endif
