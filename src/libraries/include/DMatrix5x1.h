 class DMatrix5x1{
  public:
    DMatrix5x1(){
      for (unsigned int j=0;j<3;j++){
	mA.v[j]=_mm_setzero_pd();
      }
    }
    DMatrix5x1(__m128d v1,__m128d v2,__m128d v3){
      mA.v[0]=v1;
      mA.v[1]=v2;
      mA.v[2]=v3;
    }
    DMatrix5x1(double a1, double a2, double a3, double a4, double a5){
      mA.d[0]=a1;
      mA.d[1]=a2;
      mA.d[2]=a3;
      mA.d[3]=a4;
      mA.d[4]=a5;
      mA.d[5]=0.;
    }
    ~DMatrix5x1(){};

    __m128d GetV(int pair) const{
      return mA.v[pair];
    }

    // Copy constructor
    DMatrix5x1(const DMatrix5x1 &m2){
      mA.v[0]=m2.GetV(0);
      mA.v[1]=m2.GetV(1);
      mA.v[2]=m2.GetV(2);
    }
      
    // Access by row number
    double &operator() (int row){
      return mA.d[row];
    } 
    double operator() (int row) const{
      return mA.d[row];
    }
    // Assignment operator
    DMatrix5x1 &operator=(const DMatrix5x1 &m2){
      mA.v[0]=m2.GetV(0);
      mA.v[1]=m2.GetV(1);
      mA.v[2]=m2.GetV(2);
      return *this;
    }
    
    // Matrix addition
    DMatrix5x1 operator+(const DMatrix5x1 &m2) const{
      return DMatrix5x1(_mm_add_pd(GetV(0),m2.GetV(0)),
			_mm_add_pd(GetV(1),m2.GetV(1)),
			_mm_add_pd(GetV(2),m2.GetV(2)));
    }
    DMatrix5x1 &operator+=(const DMatrix5x1 &m2){
      mA.v[0]=_mm_add_pd(GetV(0),m2.GetV(0));
      mA.v[1]=_mm_add_pd(GetV(1),m2.GetV(1)); 
      mA.v[2]=_mm_add_pd(GetV(2),m2.GetV(2));
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
	cout <<"   "<<i<<" |" <<  setw(11)<<setprecision(4) << mA.d[i] << endl;
      }      
    }
    
  private:
    union dvec{
      __m128d v[3];
      double d[6];
    };
    union dvec mA;
    
  };

// Scale 5x1 matrix by a floating point number
inline DMatrix5x1 operator*(const double c,const DMatrix5x1 &M){ 
  __m128d scale=_mm_set1_pd(c);
  return DMatrix5x1(_mm_mul_pd(scale,M.GetV(0)),_mm_mul_pd(scale,M.GetV(1)),
		    _mm_mul_pd(scale,M.GetV(2)));
}
