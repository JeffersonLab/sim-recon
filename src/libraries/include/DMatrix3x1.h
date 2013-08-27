#include <align_16.h>

#ifndef USE_SSE2

// Matrix class without SIMD instructions
class DMatrix3x1{
public:
  DMatrix3x1(){
    for (unsigned int i=0;i<3;i++) mA[i]=0.;
  }
  DMatrix3x1(double a1, double a2, double a3){
    mA[0]=a1;
    mA[1]=a2;
    mA[2]=a3;
  }
  ~DMatrix3x1(){};

  // Access by row number
  double &operator() (int row){
    return mA[row];
  } 
  double operator() (int row) const{
    return mA[row];
  }
  // Copy constructor
  DMatrix3x1(const DMatrix3x1 &m2){
    for (unsigned int i=0;i<3;i++){
      mA[i]=m2(i);
    }
  }  
  // Assignment operator
  DMatrix3x1 &operator=(const DMatrix3x1 &m2){
    for (unsigned int i=0;i<3;i++){
      mA[i]=m2(i);
    }
    return *this;
  }
  // Matrix addition
  DMatrix3x1 operator+(const DMatrix3x1 &m2) const{
    return DMatrix3x1(mA[0]+m2(0),mA[1]+m2(1),mA[2]+m2(2));

  }
  DMatrix3x1 &operator+=(const DMatrix3x1 &m2){
    for (unsigned int i=0;i<3;i++){
      mA[i]+=m2(i);
    }
    return *this;
  }

  // Matrix subtraction
  DMatrix3x1 operator-(const DMatrix3x1 &m2) const{
    return DMatrix3x1(mA[0]-m2(0),mA[1]-m2(1),mA[2]-m2(2));

  }

  // Square of the magnitude
  double Mag2() const{
    return mA[0]*mA[0]+mA[1]*mA[1]+mA[2]*mA[2];
  }
 

  
  void Print(){
      cout << "DMatrix3x1:" <<endl;
      cout << "     |      0    |" <<endl;
      cout << "----------------------" <<endl;
      for (unsigned int i=0;i<3;i++){
	cout <<"   "<<i<<" |" <<  setw(11)<<setprecision(4) << mA[i] << endl;
      }      
    }

private:
  double mA[3];
};

// Scale 3x1 matrix by a floating point number
inline DMatrix3x1 operator*(const double c,const DMatrix3x1 &M){ 
  return DMatrix3x1(c*M(0),c*M(1),c*M(2));
}

#else

// Matrix class with SIMD instructions

 class DMatrix3x1{
  public:
    DMatrix3x1()
    : mA( ALIGNED_16_BLOCK_PTR(union dvec, 1, mA) )
    {
      for (unsigned int j=0;j<2;j++){
	mA->v[j]=_mm_setzero_pd();
      }
    }
    DMatrix3x1(__m128d v1, __m128d v2)
    : mA( ALIGNED_16_BLOCK_PTR(union dvec, 1, mA) )
    {
      mA->v[0]=v1;
      mA->v[1]=v2;
    }
    DMatrix3x1(double a1, double a2, double a3)
    : mA( ALIGNED_16_BLOCK_PTR(union dvec, 1, mA) )
    {
      mA->d[0]=a1;
      mA->d[1]=a2;
      mA->d[2]=a3;
    }
    ~DMatrix3x1(){};

    __m128d GetV(int pair) const{
      return mA->v[pair];
    }

    // Copy constructor
    DMatrix3x1(const DMatrix3x1 &m2)
    : mA( ALIGNED_16_BLOCK_PTR(union dvec, 1, mA) )
    {
      mA->v[0]=m2.GetV(0);
      mA->v[1]=m2.GetV(1);
    }
    // Access by row number
    double &operator() (int row){
      return mA->d[row];
    } 
    double operator() (int row) const{
      return mA->d[row];
    }
    // Assignment operator
    DMatrix3x1 &operator=(const DMatrix3x1 &m2){
      mA->v[0]=m2.GetV(0);
      mA->v[1]=m2.GetV(1);
      return *this;
    }
    
    // Matrix addition
    DMatrix3x1 operator+(const DMatrix3x1 &m2) const{
      return DMatrix3x1(_mm_add_pd(GetV(0),m2.GetV(0)),
			_mm_add_pd(GetV(1),m2.GetV(1)));
    }
    DMatrix3x1 &operator+=(const DMatrix3x1 &m2){
      mA->v[0]=_mm_add_pd(GetV(0),m2.GetV(0));
      mA->v[1]=_mm_add_pd(GetV(1),m2.GetV(1)); 
      return *this;
    }

    // Matrix subtraction
    DMatrix3x1 operator-(const DMatrix3x1 &m2) const{
      return DMatrix3x1(_mm_sub_pd(GetV(0),m2.GetV(0)),
			_mm_sub_pd(GetV(1),m2.GetV(1)));
    }

    double Mag2() const{
      return mA->d[0]*mA->d[0]+mA->d[1]*mA->d[1]+mA->d[2]*mA->d[2];
    }
    
    void Print(){
      cout << "DMatrix3x1:" <<endl;
      cout << "     |      0    |" <<endl;
      cout << "----------------------" <<endl;
      for (unsigned int i=0;i<3;i++){
	cout <<"   "<<i<<" |" <<  setw(11)<<setprecision(4) << mA->d[i] << endl;
      }      
    }
    
  private:
    union dvec{
      __m128d v[2];
      double d[4];
    };
    ALIGNED_16_BLOCK(union dvec, 1, mA)
    
  };

// Scale 3x1 matrix by a floating point number
inline DMatrix3x1 operator*(const double c,const DMatrix3x1 &M){ 
  ALIGNED_16_BLOCK_WITH_PTR(__m128d, 1, p)
  __m128d &scale=p[0];
  scale=_mm_set1_pd(c);
  return DMatrix3x1(_mm_mul_pd(scale,M.GetV(0)),
                    _mm_mul_pd(scale,M.GetV(1)));
}

#endif
