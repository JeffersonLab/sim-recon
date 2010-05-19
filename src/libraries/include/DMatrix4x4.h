class DMatrix4x4{
 public:
  DMatrix4x4(){
    for (unsigned int i=0;i<4;i++){
      mA[i].v[0]=_mm_setzero_pd();
      mA[i].v[1]=_mm_setzero_pd();
    }
  }
  DMatrix4x4(__m128d m11,__m128d m12,__m128d m13,__m128d m14,
	     __m128d m21,__m128d m22,__m128d m23,__m128d m24){

    mA[0].v[0]=m11;
    mA[1].v[0]=m12;
    mA[2].v[0]=m13;
    mA[3].v[0]=m14; 
    mA[0].v[1]=m21;
    mA[1].v[1]=m22;
    mA[2].v[1]=m23;
    mA[3].v[1]=m24;
  }

  DMatrix4x4(const DMatrix2x2 &m1,const DMatrix2x2 &m2,
	     const DMatrix2x2 &m3,const DMatrix2x2 &m4){
    mA[0].v[0]=m1.GetV(0);
    mA[1].v[0]=m1.GetV(1);
    mA[0].v[1]=m3.GetV(0);
    mA[1].v[1]=m3.GetV(1);
    mA[2].v[0]=m2.GetV(0);
    mA[3].v[0]=m2.GetV(1);
    mA[2].v[1]=m4.GetV(0);
    mA[3].v[1]=m4.GetV(1);
  }
  ~DMatrix4x4(){}; 
  __m128d GetV(int pair, int col) const{
      return mA[col].v[pair];
  }
  void SetV(int pair,int col,__m128d v){
    mA[col].v[pair]=v;
  }
  double &operator() (int row, int col){
    return mA[col].d[row];
  } 
  double operator() (int row, int col) const{
    return mA[col].d[row];
  }

  DMatrix4x4 operator-(){
    __m128d zero=_mm_setzero_pd();
    return DMatrix4x4(_mm_sub_pd(zero,GetV(0,0)),
		      _mm_sub_pd(zero,GetV(0,1)),
		      _mm_sub_pd(zero,GetV(0,2)),
		      _mm_sub_pd(zero,GetV(0,3)),
		      _mm_sub_pd(zero,GetV(1,0)),
		      _mm_sub_pd(zero,GetV(1,1)),
		      _mm_sub_pd(zero,GetV(1,2)),
		      _mm_sub_pd(zero,GetV(1,3)));
  }

  DMatrix4x4 operator-(const DMatrix4x4 &m2){
    return DMatrix4x4(_mm_sub_pd(GetV(0,0),m2.GetV(0,0)),
		      _mm_sub_pd(GetV(0,1),m2.GetV(0,1)),
		      _mm_sub_pd(GetV(0,2),m2.GetV(0,2)),
		      _mm_sub_pd(GetV(0,3),m2.GetV(0,3)),
		      _mm_sub_pd(GetV(1,0),m2.GetV(1,0)),
		      _mm_sub_pd(GetV(1,1),m2.GetV(1,1)),
		      _mm_sub_pd(GetV(1,2),m2.GetV(1,2)),
		      _mm_sub_pd(GetV(1,3),m2.GetV(1,3)));
  }


  DMatrix4x2 operator*(const DMatrix4x2 &m2){
    __m128d m11=_mm_set1_pd(m2(0,0));
    __m128d m12=_mm_set1_pd(m2(0,1)); 
    __m128d m21=_mm_set1_pd(m2(1,0));
    __m128d m22=_mm_set1_pd(m2(1,1));  
    __m128d m31=_mm_set1_pd(m2(2,0));
    __m128d m32=_mm_set1_pd(m2(2,1)); 
    __m128d m41=_mm_set1_pd(m2(3,0));
    __m128d m42=_mm_set1_pd(m2(3,1));
    return 
      DMatrix4x2(_mm_add_pd(_mm_mul_pd(GetV(0,0),m11),
			    _mm_add_pd(_mm_mul_pd(GetV(0,1),m21),
				       _mm_add_pd(_mm_mul_pd(GetV(0,2),m31),
						  _mm_mul_pd(GetV(0,3),m41)))),
		 _mm_add_pd(_mm_mul_pd(GetV(0,0),m12),
			    _mm_add_pd(_mm_mul_pd(GetV(0,1),m22),
				       _mm_add_pd(_mm_mul_pd(GetV(0,2),m32),
						  _mm_mul_pd(GetV(0,3),m42)))),
		 _mm_add_pd(_mm_mul_pd(GetV(1,0),m11),
			    _mm_add_pd(_mm_mul_pd(GetV(1,1),m21),
				       _mm_add_pd(_mm_mul_pd(GetV(1,2),m31),
						  _mm_mul_pd(GetV(1,3),m41)))),
		 _mm_add_pd(_mm_mul_pd(GetV(1,0),m12),
			    _mm_add_pd(_mm_mul_pd(GetV(1,1),m22),
				       _mm_add_pd(_mm_mul_pd(GetV(1,2),m32),
						  _mm_mul_pd(GetV(1,3),m42)))));

		 }


  DMatrix4x4 &operator=(const DMatrix4x4 &m1){
    for (unsigned int i=0;i<4;i++){
      mA[i].v[0]=m1.GetV(0,i);
      mA[i].v[1]=m1.GetV(1,i);
    }    
    return *this;
  }
  
  DMatrix4x4 Invert(){
    DMatrix2x2 F(GetV(0,0),GetV(0,1));
    DMatrix2x2 Finv=F.Invert();
    DMatrix2x2 G(GetV(0,2),GetV(0,3)); 
    DMatrix2x2 H(GetV(1,0),GetV(1,1));
    DMatrix2x2 J(GetV(1,2),GetV(1,3));
    DMatrix2x2 Jinv=J.Invert();
    DMatrix2x2 FF=(F-G*Jinv*H).Invert();
    DMatrix2x2 JJ=(J-H*Finv*G).Invert();
    return DMatrix4x4(FF,-FF*G*Jinv,-JJ*H*Finv,JJ);     
  }
  

  
  void Print(){
    cout << "DMatrix4x4:" <<endl;
    cout << "         |      0    |      1    |      2    |      3    |" <<endl;
    cout << "----------------------------------------------------------" <<endl;
    
    for (unsigned int i=0;i<4;i++){
      for (unsigned int j=0;j<4;j++){
	cout << mA[j].d[i] <<" "; 
      } 
      cout << endl;
    }      
  }
  
 private:
  union dvec{
    __m128d v[2];
    double d[4];
  };
  union dvec mA[4];
};
