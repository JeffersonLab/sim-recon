#include <align_16.h>

#ifndef USE_SSE2

// Matrix class without SIMD instructions

class DMatrix2x2{
 public:
  DMatrix2x2(){
    mA[0][0]=mA[0][1]=mA[1][0]=mA[1][1]=0.;
  }
  ~DMatrix2x2(){};
  DMatrix2x2(double d11, double d12, double d21,
	     double d22){
    mA[0][0]=d11;
    mA[0][1]=d12;
    mA[1][0]=d21;
    mA[1][1]=d22;
  } 
  // Access by row and column
  double &operator() (int row, int col){
    return mA[row][col];
  }  
  double operator() (int row, int col) const{
    return mA[row][col];
  }
  // Assignment operator
  DMatrix2x2 &operator=(const DMatrix2x2 &m1){
    mA[0][0]=m1(0,0);
    mA[1][0]=m1(1,0);
    mA[0][1]=m1(0,1);
    mA[1][1]=m1(1,1);
    return *this;
  }
    
  // unary minus
  DMatrix2x2 operator-() const{
    return DMatrix2x2(-mA[0][0],-mA[0][1],-mA[1][0],-mA[1][1]);
  }  
  // Matrix multiplication:  (2x2) x (2x2)
  DMatrix2x2 operator*(const DMatrix2x2 &m2){
    return DMatrix2x2(mA[0][0]*m2(0,0)+mA[0][1]*m2(1,0),
		       mA[0][0]*m2(0,1)+mA[0][1]*m2(1,1),
		       mA[1][0]*m2(0,0)+mA[1][1]*m2(1,0),
		       mA[1][0]*m2(0,1)+mA[1][1]*m2(1,1));
  }
  // Matrix multiplication:  (2x2) x (2x1)
  DMatrix2x1 operator*(const DMatrix2x1 &m2){
    return DMatrix2x1(mA[0][0]*m2(0)+mA[0][1]*m2(1),
		       mA[1][0]*m2(0)+mA[1][1]*m2(1));

  }
  // Matrix addition
  DMatrix2x2 operator+(const DMatrix2x2 &m2){
    return DMatrix2x2(mA[0][0]+m2(0,0),mA[0][1]+m2(0,1),mA[1][0]+m2(1,0),
		       mA[1][1]+m2(1,1));
  }  
  // Matrix subtraction
  DMatrix2x2 operator-(const DMatrix2x2 &m2){
    return DMatrix2x2(mA[0][0]-m2(0,0),mA[0][1]-m2(0,1),mA[1][0]-m2(1,0),
		       mA[1][1]-m2(1,1));
  }
  // Matrix inversion
  DMatrix2x2 Invert(){
    double one_over_det=1./(mA[0][0]*mA[1][1]-mA[0][1]*mA[1][0]);
    return DMatrix2x2(one_over_det*mA[1][1],-one_over_det*mA[0][1],
		       -one_over_det*mA[1][0],one_over_det*mA[0][0]);

  }
  
  // Compute the determinant
  double Determinant(){
    return mA[0][0]*mA[1][1]-mA[0][1]*mA[1][0];
  }


  //Compute the chi2 contribution for a pair of hits with residual R and covariance "this"
  double Chi2(const DMatrix2x1 &R) const{
      return ( (R(0)*R(0)*mA[1][1]-R(0)*R(1)*(mA[0][1]+mA[1][0])
		+R(1)*R(1)*mA[0][0])/
	       (mA[0][0]*mA[1][1]-mA[0][1]*mA[1][0]));
    }

   void Print(){
     cout << "DMatrix2x2:" <<endl;
     cout << "     |      0    |      1    |" <<endl;
     cout << "----------------------------------" <<endl;
     for (unsigned int i=0;i<2;i++){
	cout <<"   "<<i<<" |"<< setw(11)<<setprecision(4) << mA[i][0] <<" "
	     <<  setw(11)<<setprecision(4)<< mA[i][1]<< endl;
     }      
    }


 private:
  double mA[2][2];

};

#else

// Matrix class with SIMD instructions

class DMatrix2x2{
  public:
    DMatrix2x2()
    : mA( ALIGNED_16_BLOCK_PTR(union dvec, 2, mA) )
    {
      mA[0].v=_mm_setzero_pd();
      mA[1].v=_mm_setzero_pd();
    }
    DMatrix2x2(double d11, double d12, double d21,
	       double d22)
    : mA( ALIGNED_16_BLOCK_PTR(union dvec, 2, mA) )
    {
      mA[0].v=_mm_setr_pd(d11,d21);
      mA[1].v=_mm_setr_pd(d12,d22);
    }
    DMatrix2x2(__m128d v1, __m128d v2)
    : mA( ALIGNED_16_BLOCK_PTR(union dvec, 2, mA) )
    {
      mA[0].v=v1;
      mA[1].v=v2;
    }
    DMatrix2x2(const DMatrix2x2& dm)
    : mA( ALIGNED_16_BLOCK_PTR(union dvec, 2, mA) )
    {
      mA[0].v=dm.mA[0].v;
      mA[1].v=dm.mA[1].v;
    }
    ~DMatrix2x2(){};

    __m128d GetV(int col) const{
      return mA[col].v;
    }

    // Access by row and column
    double &operator() (int row, int col){
      return mA[col].d[row];
    }  
    double operator() (int row, int col) const{
      return mA[col].d[row];
    }

    // Assignment operator
    DMatrix2x2 &operator=(const DMatrix2x2 &m1){
      //mA[0].v=m1.GetV(0);
      //mA[1].v=m1.GetV(1);   
      mA[0].v=m1.mA[0].v;
      mA[1].v=m1.mA[1].v;
      return *this;
    }
    
    // unary minus
    DMatrix2x2 operator-() const{
      ALIGNED_16_BLOCK_WITH_PTR(__m128d, 1, p);
      __m128d &zero=p[0];
      zero=_mm_setzero_pd();
      return DMatrix2x2(_mm_sub_pd(zero,mA[0].v),
			_mm_sub_pd(zero,mA[1].v));
    }
    // Matrix multiplication:  (2x2) x (2x2)
    DMatrix2x2 operator*(const DMatrix2x2 &m2){
      return 
	DMatrix2x2(_mm_add_pd(_mm_mul_pd(mA[0].v,
					 _mm_set1_pd(m2(0,0))),
			      _mm_mul_pd(mA[1].v,
					 _mm_set1_pd(m2(1,0)))), 
		   _mm_add_pd(_mm_mul_pd(mA[0].v,
					 _mm_set1_pd(m2(0,1))),
			      _mm_mul_pd(mA[1].v,
					 _mm_set1_pd(m2(1,1)))));
    }
    // Matrix multiplication:  (2x2) x (2x1)
    DMatrix2x1 operator*(const DMatrix2x1 &m2){
      return DMatrix2x1(_mm_add_pd(_mm_mul_pd(GetV(0),_mm_set1_pd(m2(0))),
				   _mm_mul_pd(GetV(1),_mm_set1_pd(m2(1)))));
    }

    // Matrix addition
    DMatrix2x2 operator+(const DMatrix2x2 &m2){
      return DMatrix2x2(_mm_add_pd(mA[0].v,m2.mA[0].v),
			_mm_add_pd(mA[1].v,m2.mA[1].v));
    }		     
    // matrix subtraction
    DMatrix2x2 operator-(const DMatrix2x2 &m2){
      return DMatrix2x2(_mm_sub_pd(mA[0].v,m2.mA[0].v),
			_mm_sub_pd(mA[1].v,m2.mA[1].v));
    }
    
    // Matrix inversion
    DMatrix2x2 Invert(){
      ALIGNED_16_BLOCK_WITH_PTR(__m128d, 1, p);
      __m128d &scale=p[0];
      scale=_mm_set1_pd(1./(mA[0].d[0]*mA[1].d[1]-mA[1].d[0]*mA[0].d[1]));
      return DMatrix2x2( _mm_mul_pd(scale,_mm_setr_pd(mA[1].d[1],-mA[0].d[1])),
		      _mm_mul_pd(scale,_mm_setr_pd(-mA[1].d[0],mA[0].d[0])));

    } 
    // Compute the determinant
    double Determinant(){
    return mA[0].d[0]*mA[1].d[1]-mA[0].d[1]*mA[1].d[0];
  }

   
    //Compute the chi2 contribution for a pair of hits with residual R and covariance "this"
    double Chi2(const DMatrix2x1 &R) const{
      return ( (R(0)*R(0)*mA[1].d[1]-R(0)*R(1)*(mA[0].d[1]+mA[1].d[0])
		+R(1)*R(1)*mA[0].d[0])/
	       (mA[0].d[0]*mA[1].d[1]-mA[0].d[1]*mA[1].d[0]));
    }
    
    void Print(){
      cout << "DMatrix2x2:" <<endl;
      cout << "     |      0    |      1    |" <<endl;
      cout << "----------------------------------" <<endl;
      for (unsigned int i=0;i<2;i++){
	cout <<"   "<<i<<" |"<< setw(11)<<setprecision(4) << mA[0].d[i] <<" "
	     <<  setw(11)<<setprecision(4)<< mA[1].d[i]<< endl;
      }      
    }
    
    // private:
    union dvec{
      __m128d v;
      double d[2];
    };
    ALIGNED_16_BLOCK(union dvec, 2, mA)
};
#endif
