#include <align_16.h>

#ifndef USE_SSE2

// Matrix class without SIMD instructions

class DMatrix5x2{
public:
  DMatrix5x2(){
    for (unsigned int i=0;i<5;i++){
      for (unsigned int j=0;j<2;j++){
	mA[i][j]=0.;
      }
    }
  }
  DMatrix5x2(double A1,double A2,double B1,double B2,double C1,double C2,
	     double D1,double D2,double E1,double E2){
    mA[0][0]=A1;
    mA[0][1]=A2;
    mA[1][0]=B1;
    mA[1][1]=B2;  
    mA[2][0]=C1;
    mA[2][1]=C2;
    mA[3][0]=D1;
    mA[3][1]=D2; 
    mA[4][0]=E1;
    mA[4][1]=E2;
  }
  ~DMatrix5x2(){};

  // Access by indices
  double &operator() (int row,int col){
    return mA[row][col];
  } 
  double operator() (int row,int col) const{
    return mA[row][col];
  }

  // Assignment operator
  DMatrix5x2 &operator=(const DMatrix5x2 &m2){
    for (unsigned int i=0;i<5;i++){
      mA[i][0]=m2(i,0);
      mA[i][1]=m2(i,1);
    }
    return *this;
  }
  


  // multiplication:  (5x2) x (2x1)
  DMatrix5x1 operator*(const DMatrix2x1 &m2){
    return DMatrix5x1(mA[0][0]*m2(0)+mA[0][1]*m2(1),
		      mA[1][0]*m2(0)+mA[1][1]*m2(1),
		      mA[2][0]*m2(0)+mA[2][1]*m2(1),
		      mA[3][0]*m2(0)+mA[3][1]*m2(1), 
		      mA[4][0]*m2(0)+mA[4][1]*m2(1));
  }
  // multiplication:  (5x2) x (2x2)
  DMatrix5x2 operator*(const DMatrix2x2 &m2){
    return DMatrix5x2(mA[0][0]*m2(0,0)+mA[0][1]*m2(1,0),
		       mA[0][0]*m2(0,1)+mA[0][1]*m2(1,1),

		       mA[1][0]*m2(0,0)+mA[1][1]*m2(1,0), 
		       mA[1][0]*m2(0,1)+mA[1][1]*m2(1,1),

		       mA[2][0]*m2(0,0)+mA[2][1]*m2(1,0),
		       mA[2][0]*m2(0,1)+mA[2][1]*m2(1,1),

		       mA[3][0]*m2(0,0)+mA[3][1]*m2(1,0), 
		       mA[3][0]*m2(0,1)+mA[3][1]*m2(1,1), 

		       mA[4][0]*m2(0,0)+mA[4][1]*m2(1,0),
		       mA[4][0]*m2(0,1)+mA[4][1]*m2(1,1)
		       );
  }

  void Print(){
    cout << "DMatrix5x2:" <<endl;
    cout << "     |      0    |      1    |" <<endl;
    cout << "----------------------------------" <<endl;
    for (unsigned int i=0;i<5;i++){
      cout <<"   "<<i<<" |"<< setw(11)<<setprecision(4) << mA[i][0] 
	   << setw(11)<<setprecision(4) << mA[i][1]<< endl;
    }      
  }

private:
  double mA[5][2];
};


#else

// Matrix class with SIMD instructions

class DMatrix5x2{
 public:
    DMatrix5x2()
    : mA( ALIGNED_16_BLOCK_PTR(union dvec, 2, mA) )
    {
      for (unsigned int j=0;j<3;j++){
	mA[0].v[j]=_mm_setzero_pd();
	mA[1].v[j]=_mm_setzero_pd();
      }
    }
    DMatrix5x2(__m128d aa, __m128d ab,
               __m128d ba, __m128d bb,
               __m128d ca, __m128d cb)
    : mA( ALIGNED_16_BLOCK_PTR(union dvec, 2, mA) )
    {
      mA[0].v[0]=aa;
      mA[0].v[1]=ba;
      mA[0].v[2]=ca;
      mA[1].v[0]=ab;
      mA[1].v[1]=bb; 
      mA[1].v[2]=cb;
    }
    DMatrix5x2(const DMatrix5x2& dm)
    : mA( ALIGNED_16_BLOCK_PTR(union dvec, 2, mA) )
    {
      mA[0].v[0]=dm.mA[0].v[0];
      mA[0].v[1]=dm.mA[0].v[1];
      mA[0].v[2]=dm.mA[0].v[2];
      mA[1].v[0]=dm.mA[1].v[0];
      mA[1].v[1]=dm.mA[1].v[1];
      mA[1].v[2]=dm.mA[1].v[2];
    }
    ~DMatrix5x2(){};

    __m128d GetV(int pair, int col) const{
      return mA[col].v[pair];
    }
    void SetV(int pair, int col, __m128d v){
      mA[col].v[pair]=v;
    }
      
    // Assignment
    DMatrix5x2& operator=(const DMatrix5x2& dm){
      mA[0].v[0]=dm.mA[0].v[0];
      mA[0].v[1]=dm.mA[0].v[1];
      mA[0].v[2]=dm.mA[0].v[2];
      mA[1].v[0]=dm.mA[1].v[0];
      mA[1].v[1]=dm.mA[1].v[1];
      mA[1].v[2]=dm.mA[1].v[2];
      return *this;
    }

    double &operator() (int row,int col){
      return mA[col].d[row];
    } 
    double operator() (int row,int col) const{
      return mA[col].d[row];
    }
    
    DMatrix5x1 operator*(const DMatrix2x1 &m2){
      ALIGNED_16_BLOCK_WITH_PTR(__m128d, 2, p)
      __m128d &a=p[0];
      __m128d &b=p[1];
      a=_mm_set1_pd(m2(0));
      b=_mm_set1_pd(m2(1));
      return DMatrix5x1(_mm_add_pd(_mm_mul_pd(GetV(0,0),a),
				   _mm_mul_pd(GetV(0,1),b)),
			_mm_add_pd(_mm_mul_pd(GetV(1,0),a),
				   _mm_mul_pd(GetV(1,1),b)),
			_mm_add_pd(_mm_mul_pd(GetV(2,0),a),
				   _mm_mul_pd(GetV(2,1),b)));
    }

    DMatrix5x2 operator*(const DMatrix2x2 &m2){
      ALIGNED_16_BLOCK_WITH_PTR(__m128d, 4, p)
      __m128d &a11=p[0]; // row,col
      __m128d &a12=p[1];
      __m128d &a21=p[2];
      __m128d &a22=p[3];
      a11=_mm_set1_pd(m2(0,0));
      a12=_mm_set1_pd(m2(0,1));
      a21=_mm_set1_pd(m2(1,0));
      a22=_mm_set1_pd(m2(1,1));
      return DMatrix5x2(_mm_add_pd(_mm_mul_pd(GetV(0,0),a11),
				   _mm_mul_pd(GetV(0,1),a21)),
			_mm_add_pd(_mm_mul_pd(GetV(0,0),a12),
				   _mm_mul_pd(GetV(0,1),a22)),
			_mm_add_pd(_mm_mul_pd(GetV(1,0),a11),
				   _mm_mul_pd(GetV(1,1),a21)),
			_mm_add_pd(_mm_mul_pd(GetV(1,0),a12),
				   _mm_mul_pd(GetV(1,1),a22)),	
			_mm_add_pd(_mm_mul_pd(GetV(2,0),a11),
				   _mm_mul_pd(GetV(2,1),a21)),
			_mm_add_pd(_mm_mul_pd(GetV(2,0),a12),
				   _mm_mul_pd(GetV(2,1),a22)));

    }




    void Print(){
      cout << "DMatrix5x2:" <<endl;
      cout << "     |      0    |      1    |" <<endl;
      cout << "----------------------------------" <<endl;
      for (unsigned int i=0;i<5;i++){
	cout <<"   "<<i<<" |"<< setw(11)<<setprecision(4) << mA[0].d[i] 
	     << setw(11)<<setprecision(4) << mA[1].d[i]<< endl;
      }      
    }
    
  private:
    union dvec{
      __m128d v[3];
      double d[6];
    };
    ALIGNED_16_BLOCK(union dvec, 2, mA)
  };
#endif
