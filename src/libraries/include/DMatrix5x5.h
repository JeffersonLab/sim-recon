class DMatrix5x5{
 public:
  DMatrix5x5(){
    for (unsigned int i=0;i<5;i++){
      for (unsigned int j=0;j<3;j++){
	mA[i].v[j]=_mm_setzero_pd();
      }
    }
  }
  DMatrix5x5(__m128d m11,__m128d m12,__m128d m13,__m128d m14,__m128d m15,
	     __m128d m21,__m128d m22,__m128d m23,__m128d m24,__m128d m25,
	     __m128d m31,__m128d m32,__m128d m33,__m128d m34,__m128d m35
	     ){
    mA[0].v[0]=m11;
    mA[1].v[0]=m12;  
    mA[2].v[0]=m13;
    mA[3].v[0]=m14;  
    mA[4].v[0]=m15;
    mA[0].v[1]=m21;
    mA[1].v[1]=m22;  
    mA[2].v[1]=m23;
    mA[3].v[1]=m24;  
    mA[4].v[1]=m25;
    mA[0].v[2]=m31;
    mA[1].v[2]=m32;  
    mA[2].v[2]=m33;
    mA[3].v[2]=m34;  
    mA[4].v[2]=m35;
  }
  // Constructor using block matrices from matrix inversion 
  DMatrix5x5(const DMatrix2x2 &A,const DMatrix2x3 &B,
	     const DMatrix3x2 &C,const DMatrix3x3 &D){
    mA[0].v[0]=A.GetV(0);
    mA[1].v[0]=A.GetV(1);
    for (unsigned int i=0;i<3;i++){
      mA[2+i].v[0]=B.GetV(i);
      for (unsigned int j=0;j<2;j++){
	mA[2+i].v[j+1]=D.GetV(j,i);
      }
    }
    mA[0].v[1]=C.GetV(0,0);
    mA[0].v[2]=C.GetV(1,0);
    mA[1].v[1]=C.GetV(0,1);
    mA[1].v[2]=C.GetV(1,1);
  }
  ~DMatrix5x5(){};
  
  __m128d GetV(int pair, int col) const{
    return mA[col].v[pair];
  }
  void SetV(int pair,int col,__m128d v){
    mA[col].v[pair]=v;
  }
  
  // Copy constructor
  DMatrix5x5(const DMatrix5x5 &m2){
    for (unsigned int i=0;i<5;i++){
      for (unsigned int j=0;j<3;j++){
	mA[i].v[j]=m2.GetV(j,i);
      }
    }
  }
  
  // Access by indices
  double &operator() (int row, int col){
    return mA[col].d[row];
  } 
  double operator() (int row, int col) const{
    return mA[col].d[row];
  }
  
  // Assignment operator
  DMatrix5x5 &operator=(const DMatrix5x5 &m1){
    for (unsigned int i=0;i<5;i++){
      mA[i].v[0]=m1.GetV(0,i);
      mA[i].v[1]=m1.GetV(1,i);
      mA[i].v[2]=m1.GetV(2,i);
    }    
    return *this;
  }

  // Matrix addition
  DMatrix5x5 operator+(const DMatrix5x5 &m2) const{
#define ADD(i,j) _mm_add_pd(GetV((i),(j)),m2.GetV((i),(j)))

    return DMatrix5x5(ADD(0,0),ADD(0,1),ADD(0,2),ADD(0,3),ADD(0,4),ADD(1,0),ADD(1,1),ADD(1,2),
		      ADD(1,3),ADD(1,4),ADD(2,0),ADD(2,1),ADD(2,2),ADD(2,3),ADD(2,4));
  }
  DMatrix5x5 &operator+=(const DMatrix5x5 &m2){
    for (unsigned int i=0;i<5;i++){
      for (unsigned int j=0;j<3;j++){
	mA[i].v[j]=ADD(j,i);
      }
    }
    return *this;
  }


  // Matrix subtraction
  DMatrix5x5 operator-(const DMatrix5x5 &m2) const{
#define SUB(i,j) _mm_sub_pd(GetV((i),(j)),m2.GetV((i),(j)))

    return DMatrix5x5(SUB(0,0),SUB(0,1),SUB(0,2),SUB(0,3),SUB(0,4),SUB(1,0),SUB(1,1),SUB(1,2),
		      SUB(1,3),SUB(1,4),SUB(2,0),SUB(2,1),SUB(2,2),SUB(2,3),SUB(2,4));
  }


  // Matrix multiplication:  (5x5) x (5x1)
  DMatrix5x1 operator*(const DMatrix5x1 &m2){
    __m128d a1=_mm_set1_pd(m2(0));
    __m128d a2=_mm_set1_pd(m2(1));
    __m128d a3=_mm_set1_pd(m2(2));
    __m128d a4=_mm_set1_pd(m2(3));
    __m128d a5=_mm_set1_pd(m2(4));

    return 
      DMatrix5x1(_mm_add_pd(_mm_mul_pd(GetV(0,0),a1),
			    _mm_add_pd(_mm_mul_pd(GetV(0,1),a2),
				       _mm_add_pd(_mm_mul_pd(GetV(0,2),a3),
						  _mm_add_pd(_mm_mul_pd(GetV(0,3),a4),
							     _mm_mul_pd(GetV(0,4),a5))))),
		 _mm_add_pd(_mm_mul_pd(GetV(1,0),a1),
			    _mm_add_pd(_mm_mul_pd(GetV(1,1),a2),
				       _mm_add_pd(_mm_mul_pd(GetV(1,2),a3),
						  _mm_add_pd(_mm_mul_pd(GetV(1,3),a4),
							     _mm_mul_pd(GetV(1,4),a5))))),
		 _mm_add_pd(_mm_mul_pd(GetV(2,0),a1),
			    _mm_add_pd(_mm_mul_pd(GetV(2,1),a2),
				       _mm_add_pd(_mm_mul_pd(GetV(2,2),a3),
						  _mm_add_pd(_mm_mul_pd(GetV(2,3),a4),
							     _mm_mul_pd(GetV(2,4),a5))))));


  }

  // Matrix multiplication:  (5x5) x (5x2)
  DMatrix5x2 operator*(const DMatrix5x2 &m2){
    __m128d m11=_mm_set1_pd(m2(0,0));
    __m128d m12=_mm_set1_pd(m2(0,1)); 
    __m128d m21=_mm_set1_pd(m2(1,0));
    __m128d m22=_mm_set1_pd(m2(1,1));  
    __m128d m31=_mm_set1_pd(m2(2,0));
    __m128d m32=_mm_set1_pd(m2(2,1)); 
    __m128d m41=_mm_set1_pd(m2(3,0));
    __m128d m42=_mm_set1_pd(m2(3,1)); 
    __m128d m51=_mm_set1_pd(m2(4,0));
    __m128d m52=_mm_set1_pd(m2(4,1));
    return 
      DMatrix5x2(_mm_add_pd(_mm_mul_pd(GetV(0,0),m11),
			    _mm_add_pd(_mm_mul_pd(GetV(0,1),m21),
				       _mm_add_pd(_mm_mul_pd(GetV(0,2),m31),
						  _mm_add_pd(_mm_mul_pd(GetV(0,3),m41),
							     _mm_mul_pd(GetV(0,4),m51))))),
		 _mm_add_pd(_mm_mul_pd(GetV(0,0),m12),
			    _mm_add_pd(_mm_mul_pd(GetV(0,1),m22),
				       _mm_add_pd(_mm_mul_pd(GetV(0,2),m32),
						  _mm_add_pd(_mm_mul_pd(GetV(0,3),m42),
							     _mm_mul_pd(GetV(0,4),m52))))),
		 _mm_add_pd(_mm_mul_pd(GetV(1,0),m11),
			    _mm_add_pd(_mm_mul_pd(GetV(1,1),m21),
				       _mm_add_pd(_mm_mul_pd(GetV(1,2),m31),
						  _mm_add_pd(_mm_mul_pd(GetV(1,3),m41),
							     _mm_mul_pd(GetV(1,4),m51))))),
		 _mm_add_pd(_mm_mul_pd(GetV(1,0),m12),
			    _mm_add_pd(_mm_mul_pd(GetV(1,1),m22),
				       _mm_add_pd(_mm_mul_pd(GetV(1,2),m32),
						  _mm_add_pd(_mm_mul_pd(GetV(1,3),m42),
							     _mm_mul_pd(GetV(1,4),m52))))), 
		 _mm_add_pd(_mm_mul_pd(GetV(2,0),m11),
			    _mm_add_pd(_mm_mul_pd(GetV(2,1),m21),
				       _mm_add_pd(_mm_mul_pd(GetV(2,2),m31),
						  _mm_add_pd(_mm_mul_pd(GetV(2,3),m41),
							     _mm_mul_pd(GetV(2,4),m51))))),
		 _mm_add_pd(_mm_mul_pd(GetV(2,0),m12),
			    _mm_add_pd(_mm_mul_pd(GetV(2,1),m22),
				       _mm_add_pd(_mm_mul_pd(GetV(2,2),m32),
						  _mm_add_pd(_mm_mul_pd(GetV(2,3),m42),
							     _mm_mul_pd(GetV(2,4),m52))))));
  }

  // Matrix multiplication: (5x5) x (5x5)
  DMatrix5x5 operator*(const DMatrix5x5 &m2){
    __m128d temp[5][5];
    for (unsigned int i=0;i<5;i++){
      for (unsigned int j=0;j<5;j++){
	temp[i][j]=_mm_set1_pd(m2(i,j));
      }
    }
    // Preprocessor macro for multiplying two __m128d elements together
#define MUL(i,j,k) _mm_mul_pd(GetV((i),(j)),temp[(j)][(k)])
    
    return DMatrix5x5(_mm_add_pd(MUL(0,0,0),
				 _mm_add_pd(MUL(0,1,0),
					    _mm_add_pd(MUL(0,2,0),
						       _mm_add_pd(MUL(0,3,0),
								  MUL(0,4,0))))),
		      _mm_add_pd(MUL(0,0,1),
				 _mm_add_pd(MUL(0,1,1),
					    _mm_add_pd(MUL(0,2,1),
						       _mm_add_pd(MUL(0,3,1),
								  MUL(0,4,1))))),
		      _mm_add_pd(MUL(0,0,2),
				 _mm_add_pd(MUL(0,1,2),
					    _mm_add_pd(MUL(0,2,2),
						       _mm_add_pd(MUL(0,3,2),
								  MUL(0,4,2))))),
		      _mm_add_pd(MUL(0,0,3),
				 _mm_add_pd(MUL(0,1,3),
					    _mm_add_pd(MUL(0,2,3),
						       _mm_add_pd(MUL(0,3,3),
								  MUL(0,4,3))))),  
		      _mm_add_pd(MUL(0,0,4),
				 _mm_add_pd(MUL(0,1,4),
					    _mm_add_pd(MUL(0,2,4),
						       _mm_add_pd(MUL(0,3,4),
								  MUL(0,4,4))))),
		      _mm_add_pd(MUL(1,0,0),
				 _mm_add_pd(MUL(1,1,0),
					    _mm_add_pd(MUL(1,2,0),
						       _mm_add_pd(MUL(1,3,0),
								  MUL(1,4,0))))),
		      _mm_add_pd(MUL(1,0,1),
				 _mm_add_pd(MUL(1,1,1),
					    _mm_add_pd(MUL(1,2,1),
						       _mm_add_pd(MUL(1,3,1),
								  MUL(1,4,1))))),
		      _mm_add_pd(MUL(1,0,2),
				 _mm_add_pd(MUL(1,1,2),
					    _mm_add_pd(MUL(1,2,2),
						       _mm_add_pd(MUL(1,3,2),
								  MUL(1,4,2))))),
		      _mm_add_pd(MUL(1,0,3),
				 _mm_add_pd(MUL(1,1,3),
					    _mm_add_pd(MUL(1,2,3),
						       _mm_add_pd(MUL(1,3,3),
								  MUL(1,4,3))))),  
		      _mm_add_pd(MUL(1,0,4),
				 _mm_add_pd(MUL(1,1,4),
					    _mm_add_pd(MUL(1,2,4),
						       _mm_add_pd(MUL(1,3,4),
								  MUL(1,4,4))))),
		      _mm_add_pd(MUL(2,0,0),
				 _mm_add_pd(MUL(2,1,0),
					    _mm_add_pd(MUL(2,2,0),
						       _mm_add_pd(MUL(2,3,0),
								  MUL(2,4,0))))),
		      _mm_add_pd(MUL(2,0,1),
				 _mm_add_pd(MUL(2,1,1),
					    _mm_add_pd(MUL(2,2,1),
						       _mm_add_pd(MUL(2,3,1),
								  MUL(2,4,1))))),
		      _mm_add_pd(MUL(2,0,2),
				 _mm_add_pd(MUL(2,1,2),
					    _mm_add_pd(MUL(2,2,2),
						       _mm_add_pd(MUL(2,3,2),
								  MUL(2,4,2))))),
		      _mm_add_pd(MUL(2,0,3),
				 _mm_add_pd(MUL(2,1,3),
					    _mm_add_pd(MUL(2,2,3),
						       _mm_add_pd(MUL(2,3,3),
								  MUL(2,4,3))))),  
		      _mm_add_pd(MUL(2,0,4),
				 _mm_add_pd(MUL(2,1,4),
					    _mm_add_pd(MUL(2,2,4),
						       _mm_add_pd(MUL(2,3,4),
								  MUL(2,4,4))))));    
  }

  // Find the transpose of this matrix
  DMatrix5x5 Transpose(){
#define SWAP(i,j,k,m) _mm_setr_pd(mA[(j)].d[(i)],mA[(m)].d[(k)])

    return DMatrix5x5(SWAP(0,0,0,1),SWAP(1,0,1,1),SWAP(2,0,2,1),SWAP(3,0,3,1),SWAP(4,0,4,1),
		      SWAP(0,2,0,3),SWAP(1,2,1,3),SWAP(2,2,2,3),SWAP(3,2,3,3),SWAP(4,2,4,3),
		      SWAP(0,4,0,5),SWAP(1,4,1,5),SWAP(2,4,2,5),SWAP(3,4,3,5),SWAP(4,4,4,5));
  }

  // Matrix inversion by blocks   
  DMatrix5x5 Invert(){
    DMatrix2x2 A(GetV(0,0),GetV(0,1));
    DMatrix3x2 C(GetV(1,0),GetV(1,1),GetV(2,0),GetV(2,1));
    DMatrix3x2 CAinv=C*A.Invert();
    DMatrix3x3 D(GetV(1,2),GetV(1,3),GetV(1,4),GetV(2,2),GetV(2,3),GetV(2,4));
    DMatrix2x3 B(GetV(0,2),GetV(0,3),GetV(0,4));
    DMatrix2x3 BDinv=B*D.Invert();
    DMatrix2x2 AA=(A-BDinv*C).Invert();
    DMatrix3x3 DD=(D-CAinv*B).Invert();
    return DMatrix5x5(AA,-AA*BDinv,-DD*CAinv,DD);
  }

  // Zero out the matrix
  DMatrix5x5 Zero(){
    for (unsigned int i=0;i<5;i++){
      for (unsigned int j=0;j<3;j++){
	mA[i].v[j]=_mm_setzero_pd();
      }
    }
    return *this;
  }


  void Print(){
    cout << "DMatrix5x5:" <<endl;
    cout << "     |      0    |      1    |      2    |      3    |      4    |" <<endl;
    cout << "----------------------------------------------------------------------" <<endl;
    
    for (unsigned int i=0;i<5;i++){
      cout <<"   "<< i << " |";
      for (unsigned int j=0;j<5;j++){
	cout << setw(11)<<setprecision(4)<<mA[j].d[i] <<" "; 
      } 
      cout << endl;
    }      
  }
    
  private:
  union dvec{
    __m128d v[3];
      double d[6];
  };
  union dvec mA[5];
    
};
  
// Scale 5x5 matrix by a floating point number
inline DMatrix5x5 operator*(const double c,const DMatrix5x5 &M){
  __m128d scale=_mm_set1_pd(c);

#define SCALE(i,j) _mm_mul_pd(scale,M.GetV((i),(j)))
  return DMatrix5x5(SCALE(0,0),SCALE(0,1),SCALE(0,2),SCALE(0,3),SCALE(0,4),
		    SCALE(1,0),SCALE(0,1),SCALE(1,2),SCALE(1,3),SCALE(1,4),
		    SCALE(2,0),SCALE(2,1),SCALE(2,2),SCALE(2,3),SCALE(2,4));

}



