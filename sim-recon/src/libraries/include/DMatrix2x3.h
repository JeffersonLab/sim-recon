#include <align_16.h>

#ifndef USE_SSE2

// Matrix class without SIMD instructions

class DMatrix2x3{
public:
  DMatrix2x3(){
    for (unsigned int i=0;i<2;i++){
      for (unsigned int j=0;j<3;j++){
	mA[i][j]=0.;
      }
    }
  }
  DMatrix2x3(double c11, double c12, double c13,
	     double c21, double c22, double c23){
    mA[0][0]=c11;
    mA[0][1]=c12;
    mA[0][2]=c13;
    mA[1][0]=c21;
    mA[1][1]=c22;
    mA[1][2]=c23;
  }
  ~DMatrix2x3(){};

  double &operator() (int row, int col){
    return mA[row][col];
  } 
  double operator() (int row, int col) const{
    return mA[row][col];
  }
  
  // Matrix multiplication:  (2x3) x (3x2)
  DMatrix2x2 operator*(const DMatrix3x2 &m2){
    return DMatrix2x2(mA[0][0]*m2(0,0)+mA[0][1]*m2(1,0)+mA[0][2]*m2(2,0),
		       mA[0][0]*m2(0,1)+mA[0][1]*m2(1,1)+mA[0][2]*m2(2,1),
		       mA[1][0]*m2(0,0)+mA[1][1]*m2(1,0)+mA[1][2]*m2(2,0),
		       mA[1][0]*m2(0,1)+mA[1][1]*m2(1,1)+mA[1][2]*m2(2,1)
		       );
  }
  // Matrix multiplication:  (2x3) x (3x3)
  DMatrix2x3 operator*(const DMatrix3x3 &m2){
    return DMatrix2x3(mA[0][0]*m2(0,0)+mA[0][1]*m2(1,0)+mA[0][2]*m2(2,0),
		       mA[0][0]*m2(0,1)+mA[0][1]*m2(1,1)+mA[0][2]*m2(2,1),
		       mA[0][0]*m2(0,2)+mA[0][1]*m2(1,2)+mA[0][2]*m2(2,2),
		       mA[1][0]*m2(0,0)+mA[1][1]*m2(1,0)+mA[1][2]*m2(2,0),
		       mA[1][0]*m2(0,1)+mA[1][1]*m2(1,1)+mA[1][2]*m2(2,1),
		       mA[1][0]*m2(0,2)+mA[1][1]*m2(1,2)+mA[1][2]*m2(2,2)
		       );


  }
  void Print(){
    cout << "DMatrix2x3:" <<endl;
    cout << "     |      0    |      1    |      2    |" <<endl;
    cout << "-----------------------------------------------" <<endl;
    
    for (unsigned int i=0;i<2;i++){ 
      cout <<"   "<< i << " |";
      for (unsigned int j=0;j<3;j++){
	cout << setw(11)<<setprecision(4) << mA[i][j] <<" "; 
      } 
      cout << endl;
    }      
  }


private:
  double mA[2][3];
};

#else

// Matrix class with SIMD instructions

class DMatrix2x3{
 public:
  DMatrix2x3()
  : mA( ALIGNED_16_BLOCK_PTR(union dvec, 3, mA) )
  {
    mA[0].v=_mm_setzero_pd();
    mA[1].v=_mm_setzero_pd();
    mA[2].v=_mm_setzero_pd();
  }
  DMatrix2x3(__m128d c1, __m128d c2, __m128d c3)
  : mA( ALIGNED_16_BLOCK_PTR(union dvec, 3, mA) )
  {
    mA[0].v=c1;
    mA[1].v=c2;
    mA[2].v=c3;
  }
  DMatrix2x3(const DMatrix2x3& dm)
  : mA( ALIGNED_16_BLOCK_PTR(union dvec, 3, mA) )
  {
    mA[0].v=dm.mA[0].v;
    mA[1].v=dm.mA[1].v;
    mA[2].v=dm.mA[2].v;
  }
  ~DMatrix2x3(){};

  __m128d GetV(int col) const{
    return mA[col].v;
  }

  // Assignment
  DMatrix2x3& operator=(const DMatrix2x3& dm){
    mA[0].v=dm.mA[0].v;
    mA[1].v=dm.mA[1].v;
    mA[2].v=dm.mA[2].v;
    return *this;
  }

  double &operator() (int row, int col){
    return mA[col].d[row];
  } 
  double operator() (int row, int col) const{
    return mA[col].d[row];
  }

// Preprocessor macro for multiplying two __m128d elements together
#define MUL1(i,j) _mm_mul_pd(GetV((i)),_mm_set1_pd(m2((i),(j))))


  // Matrix multiplication:  (2x3) x (3x2)
  DMatrix2x2 operator*(const DMatrix3x2 &m2){
    return DMatrix2x2(_mm_add_pd(MUL1(0,0),_mm_add_pd(MUL1(1,0),MUL1(2,0))),
		      _mm_add_pd(MUL1(0,1),_mm_add_pd(MUL1(1,1),MUL1(2,1))));

  }

  // Matrix multiplication:  (2x3) x (3x3)
  DMatrix2x3 operator*(const DMatrix3x3 &m2){
    return DMatrix2x3(_mm_add_pd(MUL1(0,0),_mm_add_pd(MUL1(1,0),MUL1(2,0))),
		      _mm_add_pd(MUL1(0,1),_mm_add_pd(MUL1(1,1),MUL1(2,1))),
		      _mm_add_pd(MUL1(0,2),_mm_add_pd(MUL1(1,2),MUL1(2,2))));
  }



  
  void Print(){
    cout << "DMatrix2x3:" <<endl;
    cout << "     |      0    |      1    |      2    |" <<endl;
    cout << "-----------------------------------------------" <<endl;
    
    for (unsigned int i=0;i<2;i++){ 
      cout <<"   "<< i << " |";
      for (unsigned int j=0;j<3;j++){
	cout << setw(11)<<setprecision(4) << mA[j].d[i] <<" "; 
      } 
      cout << endl;
    }      
  }

 private:
  union dvec{
    __m128d v;
    double d[2];
  };
  ALIGNED_16_BLOCK(union dvec, 3, mA)
};
#endif

