class DMatrix2x5{
 public:
  DMatrix2x5(){
    for (unsigned int i=0;i<5;i++){
      mA[i].v=_mm_setzero_pd();
    }
  }
  DMatrix2x5(__m128d c1,__m128d c2,__m128d c3,__m128d c4,__m128d c5){
    mA[0].v=c1;
    mA[1].v=c2;
    mA[2].v=c3;
    mA[3].v=c4;
    mA[4].v=c5;
  }
  ~DMatrix2x5(){};

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
    cout << "         |      0    |      1    |      2    |      3    |      4    |" <<endl;
    cout << "----------------------------------------------------------------------" <<endl;
    
    for (unsigned int i=0;i<2;i++){
      for (unsigned int j=0;j<5;j++){
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
  union dvec mA[5];

};


