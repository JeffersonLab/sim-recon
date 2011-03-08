#include <align_16.h>

#ifndef USE_SSE2
#error Woops, DMatrix4x2 is not yet implemented for non-SSE2 architectures

#else

class DMatrix4x2{
  public:
  DMatrix4x2()
  : mA( ALIGNED_16_BLOCK_PTR(union dvec, 2, mA) )
  {
    for (unsigned int j=0;j<2;j++){
      mA[0].v[j]=_mm_setzero_pd();
      mA[1].v[j]=_mm_setzero_pd();
    }
  }
  DMatrix4x2(double a11, double a12,
             double a21, double a22,
             double a31, double a32,
             double a41, double a42)
  : mA( ALIGNED_16_BLOCK_PTR(union dvec, 2, mA) )
  {
    mA[0].d[0]=a11;
    mA[0].d[1]=a21;
    mA[0].d[2]=a31;
    mA[0].d[3]=a41;
    mA[1].d[0]=a12;
    mA[1].d[1]=a22;
    mA[1].d[2]=a32;
    mA[1].d[3]=a42;
  }
  DMatrix4x2(const __m128d aa, const __m128d ab, const __m128d ba, const __m128d bb)
  : mA( ALIGNED_16_BLOCK_PTR(union dvec, 2, mA) )
  {
    mA[0].v[0]=aa;
    mA[0].v[1]=ba;
    mA[1].v[0]=ab;
    mA[1].v[1]=bb;
  }
  DMatrix4x2(const DMatrix4x2& dm)
  : mA( ALIGNED_16_BLOCK_PTR(union dvec, 2, mA) )
  {
    mA[0].v[0]=dm.mA[0].v[0];
    mA[0].v[1]=dm.mA[0].v[1];
    mA[1].v[0]=dm.mA[1].v[0];
    mA[1].v[1]=dm.mA[1].v[1];
  }
  ~DMatrix4x2(){};
  
  __m128d GetV(const int pair,const int col) const{
    return mA[col].v[pair];
  }
  void SetV(int pair,int col,__m128d v){
    mA[col].v[pair]=v;
  }
  
  // Assignment
  DMatrix4x2& operator=(const DMatrix4x2& dm){
    mA[0].v[0]=dm.mA[0].v[0];
    mA[0].v[1]=dm.mA[0].v[1];
    mA[1].v[0]=dm.mA[1].v[0];
    mA[1].v[1]=dm.mA[1].v[1];
    return *this;
  }

  double &operator() (int row,int col){
    return mA[col].d[row];
  } 
  double operator() (int row,int col) const{
      return mA[col].d[row];
  }

  DMatrix4x2 operator*(const DMatrix2x2 &m2){
    ALIGNED_16_BLOCK_WITH_PTR(__m128d, 4, p)
    __m128d &a11=p[0]; // row,col
    __m128d &a12=p[1];
    __m128d &a21=p[2];
    __m128d &a22=p[3];
    a11=_mm_set1_pd(m2(0,0));
    a12=_mm_set1_pd(m2(0,1));
    a21=_mm_set1_pd(m2(1,0));
    a22=_mm_set1_pd(m2(1,1));
    return DMatrix4x2(_mm_add_pd(_mm_mul_pd(GetV(0,0),a11),
				 _mm_mul_pd(GetV(0,1),a21)),
		      _mm_add_pd(_mm_mul_pd(GetV(0,0),a12),
				 _mm_mul_pd(GetV(0,1),a22)),
		      _mm_add_pd(_mm_mul_pd(GetV(1,0),a11),
				 _mm_mul_pd(GetV(1,1),a21)),
		      _mm_add_pd(_mm_mul_pd(GetV(1,0),a12),
				 _mm_mul_pd(GetV(1,1),a22)));
  }

  void Print(){
    cout << "DMatrix4x2:" <<endl;
    cout << "     |      0    |      1    |" <<endl;
    cout << "-----|-----------|-----------|" <<endl;
    for (unsigned int i=0;i<4;i++){
      cout <<"   "<<i<<" | " << mA[0].d[i] <<" "<< mA[1].d[i]<< endl;
    }      
  }
  
 private:
  union dvec{
    __m128d v[2];
    double d[4];
  };
  ALIGNED_16_BLOCK(union dvec, 2, mA)
};
#endif
