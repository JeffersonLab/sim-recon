class DMatrix4x2{
  public:
  DMatrix4x2(){
    for (unsigned int j=0;j<2;j++){
      mA[0].v[j]=_mm_setzero_pd();
      mA[1].v[j]=_mm_setzero_pd();
    }
  }
  DMatrix4x2(__m128d aa,__m128d ab,__m128d ba,__m128d bb){
    mA[0].v[0]=aa;
    mA[0].v[1]=ba;
    mA[1].v[0]=ab;
    mA[1].v[1]=bb;
  }
  ~DMatrix4x2(){};
  
  __m128d GetV(int pair,int col) const{
    return mA[col].v[pair];
  }
  void SetV(int pair,int col,__m128d v){
    mA[col].v[pair]=v;
  }
  
  double &operator() (int row,int col){
    return mA[col].d[row];
  } 
  double operator() (int row,int col) const{
      return mA[col].d[row];
  }

  DMatrix4x2 operator*(const DMatrix2x2 &m2){
    __m128d a11=_mm_set1_pd(m2(0,0)); // row,col
    __m128d a12=_mm_set1_pd(m2(0,1));
    __m128d a21=_mm_set1_pd(m2(1,0));
    __m128d a22=_mm_set1_pd(m2(1,1));

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
  union dvec mA[2];
  
};
