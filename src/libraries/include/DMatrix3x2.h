#include <align_16.h>

#ifndef USE_SSE2

// Matrix class without SIMD instructions
class DMatrix3x2{
 public:
  DMatrix3x2(){
    for (unsigned int i=0;i<2;i++){
      for (unsigned int j=0;j<3;j++){
	mA[i][j]=0.;
      }
    }
  }
  DMatrix3x2(double A1, double A2, double B1, double B2,
	     double C1, double C2){
    mA[0][0]=A1;
    mA[0][1]=A2;
    mA[1][0]=B1;
    mA[1][1]=B2;
    mA[2][0]=C1;
    mA[2][1]=C2;
  }
  ~DMatrix3x2(){};

  double &operator() (int row,int col){
    return mA[row][col];
  } 
  double operator() (int row,int col) const{
    return mA[row][col];
  }

  // Assignment operator
  DMatrix3x2 &operator=(const DMatrix3x2 &m2){
    for (unsigned int i=0;i<3;i++){
      mA[i][0]=m2(i,0);
      mA[i][1]=m2(i,1);
    }
    return *this;
  }

  // Matrix multiplication:  (3x2) x (2x1)
  DMatrix3x1 operator*(const DMatrix2x1 &m2){
    return DMatrix3x1(mA[0][0]*m2(0)+mA[0][1]*m2(1),
		      mA[1][0]*m2(0)+mA[1][1]*m2(1),
		      mA[2][0]*m2(0)+mA[2][1]*m2(1)
		      );

  }

  // Matrix multiplication:  (3x2) x (2x2)
  DMatrix3x2 operator*(const DMatrix2x2 &m2){
    return DMatrix3x2(mA[0][0]*m2(0,0)+mA[0][1]*m2(1,0),
		      mA[0][0]*m2(0,1)+mA[0][1]*m2(1,1),
		      mA[1][0]*m2(0,0)+mA[1][1]*m2(1,0),
		      mA[1][0]*m2(0,1)+mA[1][1]*m2(1,1),
		      mA[2][0]*m2(0,0)+mA[2][1]*m2(1,0),
		      mA[2][0]*m2(0,1)+mA[2][1]*m2(1,1)
		     );
  }
  void Print(){
    cout << "DMatrix3x2:" <<endl;
    cout << "     |      0    |      1    |" <<endl;
    cout << "-----|-----------|-----------|" <<endl;
    for (unsigned int i=0;i<3;i++){
      cout <<"   "<<i<<" |"<<setw(11)<<setprecision(4) << mA[i][0] <<" "
	   <<setw(11)<<setprecision(4)<< mA[i][1]<< endl;
    }      
  }
 private:
  double mA[3][2];

};

#else

// Matrix class with SIMD instructions

class DMatrix3x2{
 public:
  DMatrix3x2()
  : mA( ALIGNED_16_BLOCK_PTR(union dvec, 2, mA) )
  {
    mA[0].v[0]=_mm_setzero_pd();
    mA[1].v[0]=_mm_setzero_pd();
    mA[0].v[1]=_mm_setzero_pd();
    mA[1].v[1]=_mm_setzero_pd();
  }
  DMatrix3x2(__m128d aa, __m128d ab, __m128d ba, __m128d bb)
  : mA( ALIGNED_16_BLOCK_PTR(union dvec, 2, mA) )
  {
    // check_alignment();
    mA[0].v[0]=aa;
    mA[0].v[1]=ba;
    mA[1].v[0]=ab;
    mA[1].v[1]=bb;
  }
  DMatrix3x2(double a11, double a12,
             double a21, double a22,
             double a31, double a32)
  : mA( ALIGNED_16_BLOCK_PTR(union dvec, 2, mA) )
  {
    mA[0].d[0]=a11;
    mA[0].d[1]=a21;
    mA[0].d[2]=a31;
    mA[0].d[3]=0;
    mA[1].d[0]=a12;
    mA[1].d[1]=a22;
    mA[1].d[2]=a32;
    mA[1].d[3]=0;
  }
  DMatrix3x2(const DMatrix3x2& dm)
  : mA( ALIGNED_16_BLOCK_PTR(union dvec, 2, mA) )
  {
    mA[0].v[0]=dm.mA[0].v[0];
    mA[0].v[1]=dm.mA[0].v[1];
    mA[1].v[0]=dm.mA[1].v[0];
    mA[1].v[1]=dm.mA[1].v[1];
  }
  ~DMatrix3x2(){};
  
  __m128d GetV(int pair, int col) const{
    return mA[col].v[pair];
  }
  
  // Assignment
  DMatrix3x2& operator=(DMatrix3x2& dm){
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

  // Matrix multiplication:  (3x2) x (2x2)
  DMatrix3x2 operator*(const DMatrix2x2 &m2){
    ALIGNED_16_BLOCK_WITH_PTR(__m128d, 4, p)
    __m128d &a11=p[0];
    __m128d &a12=p[1];
    __m128d &a21=p[2];
    __m128d &a22=p[3];
    a11=_mm_set1_pd(m2(0,0)); // row,col
    a12=_mm_set1_pd(m2(0,1));
    a21=_mm_set1_pd(m2(1,0));
    a22=_mm_set1_pd(m2(1,1));
    return DMatrix3x2(_mm_add_pd(_mm_mul_pd(GetV(0,0),a11),
				 _mm_mul_pd(GetV(0,1),a21)),
		      _mm_add_pd(_mm_mul_pd(GetV(0,0),a12),
				 _mm_mul_pd(GetV(0,1),a22)),
		      _mm_add_pd(_mm_mul_pd(GetV(1,0),a11),
				 _mm_mul_pd(GetV(1,1),a21)),
		      _mm_add_pd(_mm_mul_pd(GetV(1,0),a12),
				 _mm_mul_pd(GetV(1,1),a22)));
  }

  void Print(){
    cout << "DMatrix3x2:" <<endl;
    cout << "     |      0    |      1    |" <<endl;
    cout << "-----|-----------|-----------|" <<endl;
    for (unsigned int i=0;i<3;i++){
      cout <<"   "<<i<<" |"<<setw(11)<<setprecision(4) << mA[0].d[i] <<" "
	   <<setw(11)<<setprecision(4)<< mA[1].d[i]<< endl;
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
