class DMatrix3x3{
 public:
  DMatrix3x3(){
    for (unsigned int i=0;i<3;i++){
      mA[i].v[0]=_mm_setzero_pd();
      mA[i].v[1]=_mm_setzero_pd();
    }
  }
  DMatrix3x3(__m128d m11,__m128d m12,__m128d m13,__m128d m21,__m128d m22,
	     __m128d m23){

    mA[0].v[0]=m11;
    mA[1].v[0]=m12;
    mA[2].v[0]=m13; 
    mA[0].v[1]=m21;
    mA[1].v[1]=m22;
    mA[2].v[1]=m23;
  }

  ~DMatrix3x3(){};

  __m128d GetV(int pair, int col) const{
      return mA[col].v[pair];
  }

  double &operator() (int row, int col){
    return mA[col].d[row];
  } 
  double operator() (int row, int col) const{
    return mA[col].d[row];
  }

  // unary minus
  DMatrix3x3 operator-(){
    __m128d zero=_mm_setzero_pd();
    return DMatrix3x3(_mm_sub_pd(zero,GetV(0,0)),
		      _mm_sub_pd(zero,GetV(0,1)),
		      _mm_sub_pd(zero,GetV(0,2)),
		      _mm_sub_pd(zero,GetV(1,0)),
		      _mm_sub_pd(zero,GetV(1,1)),
		      _mm_sub_pd(zero,GetV(1,2)));
  }

  // Matrix subtraction
  DMatrix3x3 operator-(const DMatrix3x3 &m2){
    return DMatrix3x3(_mm_sub_pd(GetV(0,0),m2.GetV(0,0)),
		      _mm_sub_pd(GetV(0,1),m2.GetV(0,1)),
		      _mm_sub_pd(GetV(0,2),m2.GetV(0,2)),
		      _mm_sub_pd(GetV(1,0),m2.GetV(1,0)),
		      _mm_sub_pd(GetV(1,1),m2.GetV(1,1)),
		      _mm_sub_pd(GetV(1,2),m2.GetV(1,2)));
  }



  
  // Matrix multiplication:  (3x3) x (3x2)
  DMatrix3x2 operator*(const DMatrix3x2 &m2){
    __m128d m11=_mm_set1_pd(m2(0,0));
    __m128d m12=_mm_set1_pd(m2(0,1)); 
    __m128d m21=_mm_set1_pd(m2(1,0));
    __m128d m22=_mm_set1_pd(m2(1,1));  
    __m128d m31=_mm_set1_pd(m2(2,0));
    __m128d m32=_mm_set1_pd(m2(2,1)); 
    return 
      DMatrix3x2(_mm_add_pd(_mm_mul_pd(GetV(0,0),m11),
			    _mm_add_pd(_mm_mul_pd(GetV(0,1),m21),
				       _mm_mul_pd(GetV(0,2),m31))),
		 _mm_add_pd(_mm_mul_pd(GetV(0,0),m12),
			    _mm_add_pd(_mm_mul_pd(GetV(0,1),m22),
				       _mm_mul_pd(GetV(0,2),m32))),
		 _mm_add_pd(_mm_mul_pd(GetV(1,0),m11),
			    _mm_add_pd(_mm_mul_pd(GetV(1,1),m21),
				       _mm_mul_pd(GetV(1,2),m31))),
		 _mm_add_pd(_mm_mul_pd(GetV(1,0),m12),
			    _mm_add_pd(_mm_mul_pd(GetV(1,1),m22),
				       _mm_mul_pd(GetV(1,2),m32))));
  }


  // Matrix inversion
  DMatrix3x3 Invert(){
    double b11=mA[1].d[1]*mA[2].d[2]-mA[1].d[2]*mA[2].d[1];
    double b12=mA[2].d[0]*mA[1].d[2]-mA[1].d[0]*mA[2].d[2];
    double b13=mA[1].d[0]*mA[2].d[1]-mA[1].d[1]*mA[2].d[0];
    double b21=mA[2].d[1]*mA[0].d[2]-mA[0].d[1]*mA[2].d[2];
    double b22=mA[0].d[0]*mA[2].d[2]-mA[0].d[2]*mA[2].d[0];
    double b23=mA[2].d[0]*mA[0].d[1]-mA[0].d[0]*mA[2].d[1];
    double b31=mA[0].d[1]*mA[1].d[2]-mA[1].d[1]*mA[0].d[2];
    double b32=mA[1].d[0]*mA[0].d[2]-mA[0].d[0]*mA[1].d[2];
    double b33=mA[0].d[0]*mA[1].d[1]-mA[0].d[1]*mA[1].d[0];
    __m128d one_over_detA=_mm_set1_pd(1./(mA[0].d[0]*b11+mA[1].d[0]*b21+mA[2].d[0]*b31));
    return DMatrix3x3(_mm_mul_pd(one_over_detA,_mm_setr_pd(b11,b21)),
		      _mm_mul_pd(one_over_detA,_mm_setr_pd(b12,b22)),
		      _mm_mul_pd(one_over_detA,_mm_setr_pd(b13,b23)),
		      _mm_mul_pd(one_over_detA,_mm_setr_pd(b31,0.)),
		      _mm_mul_pd(one_over_detA,_mm_setr_pd(b32,0.)),
		      _mm_mul_pd(one_over_detA,_mm_setr_pd(b33,0.)));
				
  }

  void Print(){
    cout << "DMatrix3x3:" <<endl;
    cout << "     |      0    |      1    |      2    |" <<endl;
    cout << "----------------------------------------------" <<endl;
    
    for (unsigned int i=0;i<3;i++){
      cout << "   " << i <<" |";
      for (unsigned int j=0;j<3;j++){
	cout <<setw(11)<<setprecision(4)<< mA[j].d[i] <<" "; 
      } 
      cout << endl;
    }      
  }



 private:
  union dvec{
    __m128d v[2];
    double d[4];
  };
  union dvec mA[3];
};
