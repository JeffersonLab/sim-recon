#include <align_16.h>

#ifndef USE_SSE2

// Matrix class without SIMD instructions

class DMatrix2x5{
public:
  DMatrix2x5(){
    for (unsigned int i=0;i<2;i++){
      for (unsigned int j=0;j<5;j++){
	mA[i][j]=0.;
      }
    }
  }
  DMatrix2x5(double a1, double a2, double a3, double a4, double a5,
	     double b1, double b2, double b3, double b4, double b5){
    mA[0][0]=a1;
    mA[0][1]=a2;
    mA[0][2]=a3;
    mA[0][3]=a4;
    mA[0][4]=a5;  
    mA[1][0]=b1;
    mA[1][1]=b2;
    mA[1][2]=b3;
    mA[1][3]=b4;
    mA[1][4]=b5;
  }
	      
  ~DMatrix2x5(){};

  double &operator() (int row, int col){
    return mA[row][col];
  } 
  double operator() (int row, int col) const{
    return mA[row][col];
  }  
  // Matrix multiplication:  (2x5) x (5x2)
#define MUL5(i,j) mA[(i)][0]*m2(0,(j))+mA[(i)][1]*m2(1,(j))+mA[(i)][2]*m2(2,(j))+mA[(i)][3]*m2(3,(j))+mA[(i)][4]*m2(4,(j))

  DMatrix2x2 operator*(const DMatrix5x2 &m2) const{
    return DMatrix2x2(MUL5(0,0),MUL5(0,1),MUL5(1,0),MUL5(1,1));
  }

  // Matrix multiplication: (2x5) x (5x5)
  DMatrix2x5 operator*(const DMatrix5x5 &m2) const{
    return 
      DMatrix2x5(MUL5(0,0),MUL5(0,1),MUL5(0,2),MUL5(0,3),MUL5(0,4),
		  MUL5(1,0),MUL5(1,1),MUL5(1,2),MUL5(1,3),MUL5(1,4)
		  );
  }

  void Print(){
    cout << "DMatrix2x5:" <<endl;
    cout << "     |      0    |      1    |      2    |      3    |      4    |" <<endl;
    cout << "-------------------------------------------------------------------" <<endl;  
    for (unsigned int i=0;i<2;i++){   
      cout <<"   "<< i << " |";
      for (unsigned int j=0;j<5;j++){
	cout  <<setw(11)<<setprecision(4)<< mA[i][j] <<" "; 
      } 
      cout << endl;
    }      
  }

private:
  double mA[2][5];

};

#else

// Matrix class with SIMD instructions

class DMatrix2x5{
 public:
  DMatrix2x5()
  : mA( ALIGNED_16_BLOCK_PTR(union dvec, 5, mA) )
  {
    mA[0].v=_mm_setzero_pd();
    mA[1].v=_mm_setzero_pd();
    mA[2].v=_mm_setzero_pd();
    mA[3].v=_mm_setzero_pd();
    mA[4].v=_mm_setzero_pd();
  }
  DMatrix2x5(__m128d c1, __m128d c2, __m128d c3, __m128d c4, __m128d c5)
  : mA( ALIGNED_16_BLOCK_PTR(union dvec, 5, mA) )
  {
    mA[0].v=c1;
    mA[1].v=c2;
    mA[2].v=c3;
    mA[3].v=c4;
    mA[4].v=c5;
  }
  DMatrix2x5(const DMatrix2x5& dm)
  : mA( ALIGNED_16_BLOCK_PTR(union dvec, 5, mA) )
  {
    mA[0].v=dm.mA[0].v;
    mA[1].v=dm.mA[1].v;
    mA[2].v=dm.mA[2].v;
    mA[3].v=dm.mA[3].v;
    mA[4].v=dm.mA[4].v;
  }
  ~DMatrix2x5(){};

  __m128d GetV(const int col) const{
    return mA[col].v;
  }

  // Assignment
  DMatrix2x5& operator=(const DMatrix2x5& dm){
    mA[0].v=dm.mA[0].v;
    mA[1].v=dm.mA[1].v;
    mA[2].v=dm.mA[2].v;
    mA[3].v=dm.mA[3].v;
    mA[4].v=dm.mA[4].v;
    return *this;
  }

  double &operator() (int row, int col){
    return mA[col].d[row];
  } 
  double operator() (int row, int col) const{
    return mA[col].d[row];
  }

// Preprocessor macro for multiplying two __m128d elements together
//#define MUL1(i,j) _mm_mul_pd(GetV((i)),_mm_set1_pd(m2((i),(j))))

  DMatrix2x2 operator*(const DMatrix5x2 &m2) const{
    return DMatrix2x2(_mm_add_pd(MUL1(0,0),
				 _mm_add_pd(MUL1(1,0),
					    _mm_add_pd(MUL1(2,0),
						       _mm_add_pd(MUL1(3,0),MUL1(4,0))))),
		      _mm_add_pd(MUL1(0,1),
				 _mm_add_pd(MUL1(1,1),
					    _mm_add_pd(MUL1(2,1),
						       _mm_add_pd(MUL1(3,1),MUL1(4,1))))));

  }  

  DMatrix2x5 operator*(const DMatrix5x5 &m2) const{
    return 
      DMatrix2x5(_mm_add_pd(MUL1(0,0),
			    _mm_add_pd(MUL1(1,0),
				       _mm_add_pd(MUL1(2,0),
						  _mm_add_pd(MUL1(3,0),MUL1(4,0))))),
		 _mm_add_pd(MUL1(0,1),
			    _mm_add_pd(MUL1(1,1),
				       _mm_add_pd(MUL1(2,1),
						  _mm_add_pd(MUL1(3,1),MUL1(4,1))))),
		 _mm_add_pd(MUL1(0,2),
			    _mm_add_pd(MUL1(1,2),
				       _mm_add_pd(MUL1(2,2),
						  _mm_add_pd(MUL1(3,2),MUL1(4,2))))),
		 _mm_add_pd(MUL1(0,3),
			    _mm_add_pd(MUL1(1,3),
				       _mm_add_pd(MUL1(2,3),
						  _mm_add_pd(MUL1(3,3),MUL1(4,3))))),
		 _mm_add_pd(MUL1(0,4),
			    _mm_add_pd(MUL1(1,4),
				       _mm_add_pd(MUL1(2,4),
						  _mm_add_pd(MUL1(3,4),MUL1(4,4))))));
  }



  
  void Print(){
    cout << "DMatrix2x5:" <<endl;
    cout << "     |      0    |      1    |      2    |      3    |      4    |" <<endl;
    cout << "------------------------------------------------------------------" <<endl;
    
    for (unsigned int i=0;i<2;i++){  
      cout <<"   "<< i << " |";
      for (unsigned int j=0;j<5;j++){
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
  ALIGNED_16_BLOCK(union dvec, 5, mA)
};
#endif

