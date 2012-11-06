#include <align_16.h>

#ifndef USE_SSE2

// Matrix class without SIMD instructions

class DMatrix2x4{
public:
  DMatrix2x4(){
    for (unsigned int i=0;i<2;i++){
      for (unsigned int j=0;j<4;j++){
	mA[i][j]=0.;
      }
    }
  }
  DMatrix2x4(double c11, double c12, double c13, double c14,
	     double c21, double c22, double c23, double c24){
    mA[0][0]=c11;
    mA[0][1]=c12;
    mA[0][2]=c13;
    mA[0][3]=c14;
    mA[1][0]=c21;
    mA[1][1]=c22;
    mA[1][2]=c23;
    mA[1][3]=c24;
  }
  ~DMatrix2x4(){};

  
  // Assignment
  DMatrix2x4& operator=(const DMatrix2x4 &m2){
    for (unsigned int i=0;i<4;i++){
      mA[0][i]=m2(0,i);
      mA[1][i]=m2(1,i);
    }
    return *this;
  }

  double &operator() (int row, int col){
    return mA[row][col];
  } 
  double operator() (int row, int col) const{
    return mA[row][col];
  }
  // matrix multiplication: (2x4)x(4x2)
  DMatrix2x2 operator*(const DMatrix4x2 &m2){
    return 
      DMatrix2x2(mA[0][0]*m2(0,0)+mA[0][1]*m2(1,0)+mA[0][2]*m2(2,0)+mA[0][3]*m2(3,0),
		 mA[0][0]*m2(0,1)+mA[0][1]*m2(1,1)+mA[0][2]*m2(2,1)+mA[0][3]*m2(3,1),
		 mA[1][0]*m2(0,0)+mA[1][1]*m2(1,0)+mA[1][2]*m2(2,0)+mA[1][3]*m2(3,0),
		 mA[1][0]*m2(0,1)+mA[1][1]*m2(1,1)+mA[1][2]*m2(2,1)+mA[1][3]*m2(3,1));
  }
  // Matrix multiplication:  (2x4) x (4x4)
  DMatrix2x4 operator*(const DMatrix4x4 &m2){
    return 
      DMatrix2x4(mA[0][0]*m2(0,0)+mA[0][1]*m2(1,0)+mA[0][2]*m2(2,0)+mA[0][3]*m2(3,0),
		 mA[0][0]*m2(0,1)+mA[0][1]*m2(1,1)+mA[0][2]*m2(2,1)+mA[0][3]*m2(3,1),
		 mA[0][0]*m2(0,2)+mA[0][1]*m2(1,2)+mA[0][2]*m2(2,2)+mA[0][3]*m2(3,2),
		 mA[0][0]*m2(0,3)+mA[0][1]*m2(1,3)+mA[0][2]*m2(2,3)+mA[0][3]*m2(3,3),
		 mA[1][0]*m2(0,0)+mA[1][1]*m2(1,0)+mA[1][2]*m2(2,0)+mA[1][3]*m2(3,0),
		 mA[1][0]*m2(0,1)+mA[1][1]*m2(1,1)+mA[1][2]*m2(2,1)+mA[1][3]*m2(3,1), 
		 mA[1][0]*m2(0,2)+mA[1][1]*m2(1,2)+mA[1][2]*m2(2,2)+mA[1][3]*m2(3,2),
		 mA[1][0]*m2(0,3)+mA[1][1]*m2(1,3)+mA[1][2]*m2(2,3)+mA[1][3]*m2(3,3));
  }

  void Print(){
    cout << "DMatrix2x4:" <<endl;
    cout << "     |      0    |      1    |      2    |      3    |" <<endl;
    cout << "------------------------------------------------------" <<endl;
    
    for (unsigned int i=0;i<2;i++){ 
      cout <<"   "<< i << " |";
      for (unsigned int j=0;j<4;j++){
	cout << setw(11)<<setprecision(4) << mA[i][j] <<" "; 
      } 
      cout << endl;
    }      
  }


private:
  double mA[2][4];
};

#else

class DMatrix2x4{
 public:
  DMatrix2x4()
  : mA( ALIGNED_16_BLOCK_PTR(union dvec, 4, mA) )
  {
    mA[0].v=_mm_setzero_pd();
    mA[1].v=_mm_setzero_pd();
    mA[2].v=_mm_setzero_pd();
    mA[3].v=_mm_setzero_pd();
  }
  DMatrix2x4(__m128d c1, __m128d c2, __m128d c3, __m128d c4)
  : mA( ALIGNED_16_BLOCK_PTR(union dvec, 4, mA) )
  {
    mA[0].v=c1;
    mA[1].v=c2;
    mA[2].v=c3;
    mA[3].v=c4;
  }
  DMatrix2x4(const DMatrix2x4& dm)
  : mA( ALIGNED_16_BLOCK_PTR(union dvec, 4, mA) )
  {
    mA[0].v=dm.mA[0].v;
    mA[1].v=dm.mA[1].v;
    mA[2].v=dm.mA[2].v;
    mA[3].v=dm.mA[3].v;
  }
  ~DMatrix2x4(){};

  __m128d GetV(int col) const{
    return mA[col].v;
  }

  // Assignment
  DMatrix2x4& operator=(const DMatrix2x4& dm){
    mA[0].v=dm.mA[0].v;
    mA[1].v=dm.mA[1].v;
    mA[2].v=dm.mA[2].v;
    mA[3].v=dm.mA[3].v;
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

  DMatrix2x2 operator*(const DMatrix4x2 &m2){
    return DMatrix2x2(_mm_add_pd(MUL1(0,0),
				 _mm_add_pd(MUL1(1,0),
					    _mm_add_pd(MUL1(2,0),
						       MUL1(3,0)))),
		      _mm_add_pd(MUL1(0,1),
				 _mm_add_pd(MUL1(1,1),
					    _mm_add_pd(MUL1(2,1),
						       MUL1(3,1)))));

  }

  DMatrix2x4 operator*(const DMatrix4x4 &m2){
    return DMatrix2x4(_mm_add_pd(MUL1(0,0),
				 _mm_add_pd(MUL1(1,0),
					    _mm_add_pd(MUL1(2,0),
						       MUL1(3,0)))),
		      _mm_add_pd(MUL1(0,1),
				 _mm_add_pd(MUL1(1,1),
					    _mm_add_pd(MUL1(2,1),
						       MUL1(3,1)))),
		      _mm_add_pd(MUL1(0,2),
				 _mm_add_pd(MUL1(1,2),
					    _mm_add_pd(MUL1(2,2),
						       MUL1(3,2)))),
		      _mm_add_pd(MUL1(0,3),
				 _mm_add_pd(MUL1(1,3),
					    _mm_add_pd(MUL1(2,3),
						       MUL1(3,3)))));
  }



  
  void Print(){
    cout << "DMatrix2x4:" <<endl;
    cout << "         |      0    |      1    |      2    |      3    |" <<endl;
    cout << "----------------------------------------------------------" <<endl;
    
    for (unsigned int i=0;i<2;i++){
      for (unsigned int j=0;j<4;j++){
	cout << mA[j].d[i] <<" "; 
      } 
      cout << endl;
    }      
  }

 private:
  union dvec{
    __m128d v;
    double d[2];
  };
  ALIGNED_16_BLOCK(union dvec, 4, mA)
};

#endif
