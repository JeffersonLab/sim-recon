class DMatrix2x2{
  public:
    DMatrix2x2(){
      mA[0].v=_mm_setzero_pd();
      mA[1].v=_mm_setzero_pd();
    }
    DMatrix2x2(const double d11,const double d12,const double d21,
	       const double d22){
      mA[0].v=_mm_setr_pd(d11,d21);
      mA[1].v=_mm_setr_pd(d12,d22);
    }
    DMatrix2x2(__m128d v1,__m128d v2){
      mA[0].v=v1;
      mA[1].v=v2;
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
      __m128d zero=_mm_setzero_pd();
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
      __m128d scale=_mm_set1_pd(1./(mA[0].d[0]*mA[1].d[1]-mA[1].d[0]*mA[0].d[1]));
      return DMatrix2x2( _mm_mul_pd(scale,_mm_setr_pd(mA[1].d[1],-mA[0].d[1])),
		      _mm_mul_pd(scale,_mm_setr_pd(-mA[1].d[0],mA[0].d[0])));

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
    union dvec mA[2];
};
