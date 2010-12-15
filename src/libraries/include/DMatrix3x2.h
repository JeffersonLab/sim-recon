#ifndef USE_SIMD 

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
  DMatrix3x2(const double A1,const double A2,const double B1,const double B2,
	      const double C1,const double C2){
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
  DMatrix3x2(){
    mA[0].v[0]=_mm_setzero_pd();
    mA[1].v[0]=_mm_setzero_pd();   
    mA[0].v[1]=_mm_setzero_pd();
    mA[1].v[1]=_mm_setzero_pd();
  }
  DMatrix3x2(__m128d aa,__m128d ab,__m128d ba,__m128d bb){
    mA[0].v[0]=aa;
    mA[0].v[1]=ba;
    mA[1].v[0]=ab;
    mA[1].v[1]=bb;
  }
  ~DMatrix3x2(){};
  
  __m128d GetV(int pair,int col) const{
    return mA[col].v[pair];
  }
  
  double &operator() (int row,int col){
    return mA[col].d[row];
  } 
  double operator() (int row,int col) const{
      return mA[col].d[row];
  }

  // Matrix multiplication:  (3x2) x (2x2)
  DMatrix3x2 operator*(const DMatrix2x2 &m2){
    __m128d a11=_mm_set1_pd(m2(0,0)); // row,col
    __m128d a12=_mm_set1_pd(m2(0,1));
    __m128d a21=_mm_set1_pd(m2(1,0));
    __m128d a22=_mm_set1_pd(m2(1,1));

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
  union dvec mA[2];
  
};
#endif
