class DMatrix2x3{
 public:
  DMatrix2x3(){
    for (unsigned int i=0;i<3;i++){
      mA[0].v=_mm_setzero_pd();
    }
  }
  DMatrix2x3(__m128d c1,__m128d c2,__m128d c3){
    mA[0].v=c1;
    mA[1].v=c2;
    mA[2].v=c3;
  }
  ~DMatrix2x3(){};

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
  union dvec mA[3];

};


