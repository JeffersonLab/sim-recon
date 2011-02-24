#ifndef USE_SSE2

#else

class DMatrix2x4{
 public:
  DMatrix2x4(){
    for (unsigned int i=0;i<4;i++){
      mA[0].v=_mm_setzero_pd();
    }
  }
  DMatrix2x4(__m128d c1,__m128d c2,__m128d c3,__m128d c4){
    mA[0].v=c1;
    mA[1].v=c2;
    mA[2].v=c3;
    mA[3].v=c4;
  }
  ~DMatrix2x4(){};

  __m128d GetV(int col) const{
    return mA[col].v;
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
  union dvec mA[4];

};

#endif
