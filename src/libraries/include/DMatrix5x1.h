#include <align_16.h>

#ifndef USE_SSE2

// Matrix class without SIMD instructions
class DMatrix5x1{
public:
  DMatrix5x1(){
    for (unsigned int i=0;i<5;i++) mA[i]=0.;
  }
  DMatrix5x1(double a1, double a2, double a3, double a4, double a5){
    mA[0]=a1;
    mA[1]=a2;
    mA[2]=a3;
    mA[3]=a4;
    mA[4]=a5;
  }
  ~DMatrix5x1(){};

  // Access by row number
  double &operator() (int row){
    return mA[row];
  } 
  double operator() (int row) const{
    return mA[row];
  }
  // Copy constructor
  DMatrix5x1(const DMatrix5x1 &m2){
    for (unsigned int i=0;i<5;i++){
      mA[i]=m2(i);
    }
  }  
  // Assignment operator
  DMatrix5x1 &operator=(const DMatrix5x1 &m2){
    for (unsigned int i=0;i<5;i++){
      mA[i]=m2(i);
    }
    return *this;
  }
  // Matrix addition
  DMatrix5x1 operator+(const DMatrix5x1 &m2) const{
    return DMatrix5x1(mA[0]+m2(0),mA[1]+m2(1),mA[2]+m2(2),mA[3]+m2(3),
		      mA[4]+m2(4));

  }
  DMatrix5x1 &operator+=(const DMatrix5x1 &m2){
    for (unsigned int i=0;i<5;i++){
      mA[i]+=m2(i);
    }
    return *this;
  }

  // Matrix subtraction
  DMatrix5x1 operator-(const DMatrix5x1 &m2) const{
    return DMatrix5x1(mA[0]-m2(0),mA[1]-m2(1),mA[2]-m2(2),mA[3]-m2(3),
		      mA[4]-m2(4));

  }

  bool IsFinite(){
     return (isfinite(mA[0]) && isfinite(mA[1]) && isfinite(mA[2]) && isfinite(mA[3]) && isfinite(mA[4]));
  }

  
  void Print(){
      cout << "DMatrix5x1:" <<endl;
      cout << "     |      0    |" <<endl;
      cout << "----------------------" <<endl;
      for (unsigned int i=0;i<5;i++){
	cout <<"   "<<i<<" |" <<  setw(11)<<setprecision(6) << mA[i] << endl;
      }      
    }

private:
  double mA[5];
};

// Scale 5x1 matrix by a floating point number
inline DMatrix5x1 operator*(const double c,const DMatrix5x1 &M){ 
  return DMatrix5x1(c*M(0),c*M(1),c*M(2),c*M(3),c*M(4));
}

#else

// Matrix class with SIMD instructions

 class DMatrix5x1{
  public:
    DMatrix5x1()
    : mA( ALIGNED_16_BLOCK_PTR(union dvec, 1, mA) )
    {
      for (unsigned int j=0;j<3;j++){
	mA->v[j]=_mm_setzero_pd();
      }
    }
    DMatrix5x1(__m128d v1, __m128d v2, __m128d v3)
    : mA( ALIGNED_16_BLOCK_PTR(union dvec, 1, mA) )
    {
      mA->v[0]=v1;
      mA->v[1]=v2;
      mA->v[2]=v3;
    }
    DMatrix5x1(double a1, double a2, double a3, double a4, double a5)
    : mA( ALIGNED_16_BLOCK_PTR(union dvec, 1, mA) )
    {
      mA->d[0]=a1;
      mA->d[1]=a2;
      mA->d[2]=a3;
      mA->d[3]=a4;
      mA->d[4]=a5;
      mA->d[5]=0.;
    }
    ~DMatrix5x1(){};

    __m128d GetV(int pair) const{
      return mA->v[pair];
    }

    // Copy constructor
    DMatrix5x1(const DMatrix5x1 &m2)
    : mA( ALIGNED_16_BLOCK_PTR(union dvec, 1, mA) )
    {
      mA->v[0]=m2.GetV(0);
      mA->v[1]=m2.GetV(1);
      mA->v[2]=m2.GetV(2);
    }
    // Access by row number
    double &operator() (int row){
      return mA->d[row];
    } 
    double operator() (int row) const{
      return mA->d[row];
    }
    // Assignment operator
    DMatrix5x1 &operator=(const DMatrix5x1 &m2){
      mA->v[0]=m2.GetV(0);
      mA->v[1]=m2.GetV(1);
      mA->v[2]=m2.GetV(2);
      return *this;
    }
    
    // Matrix addition
    DMatrix5x1 operator+(const DMatrix5x1 &m2) const{
      return DMatrix5x1(_mm_add_pd(GetV(0),m2.GetV(0)),
			_mm_add_pd(GetV(1),m2.GetV(1)),
			_mm_add_pd(GetV(2),m2.GetV(2)));
    }
    DMatrix5x1 &operator+=(const DMatrix5x1 &m2){
      mA->v[0]=_mm_add_pd(GetV(0),m2.GetV(0));
      mA->v[1]=_mm_add_pd(GetV(1),m2.GetV(1)); 
      mA->v[2]=_mm_add_pd(GetV(2),m2.GetV(2));
      return *this;
    }

    // Matrix subtraction
    DMatrix5x1 operator-(const DMatrix5x1 &m2) const{
      return DMatrix5x1(_mm_sub_pd(GetV(0),m2.GetV(0)),
			_mm_sub_pd(GetV(1),m2.GetV(1)),
			_mm_sub_pd(GetV(2),m2.GetV(2)));
    }
      
    void Print(){
      cout << "DMatrix5x1:" <<endl;
      cout << "     |      0    |" <<endl;
      cout << "----------------------" <<endl;
      for (unsigned int i=0;i<5;i++){
	cout <<"   "<<i<<" |" <<  setw(11)<<setprecision(4) << mA->d[i] << endl;
      }      
    }
    
  private:
    union dvec{
      __m128d v[3];
      double d[6];
    };
    ALIGNED_16_BLOCK(union dvec, 1, mA)
    
  };

// Scale 5x1 matrix by a floating point number
inline DMatrix5x1 operator*(const double c,const DMatrix5x1 &M){ 
  ALIGNED_16_BLOCK_WITH_PTR(__m128d, 1, p)
  __m128d &scale=p[0];
  scale=_mm_set1_pd(c);
  return DMatrix5x1(_mm_mul_pd(scale,M.GetV(0)),
                    _mm_mul_pd(scale,M.GetV(1)),
		    _mm_mul_pd(scale,M.GetV(2)));
}

#endif
