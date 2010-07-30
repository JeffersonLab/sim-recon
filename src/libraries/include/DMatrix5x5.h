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
  // Constructor for a symmetric matrix
  DMatrix5x5(__m128d m11,__m128d m12,__m128d m13,__m128d m14,__m128d m15,
	     __m128d m23,__m128d m24,__m128d m25,__m128d m35
	     ){
    mA[0].v[0]=m11;
    mA[1].v[0]=m12;  
    mA[2].v[0]=m13;
    mA[3].v[0]=m14;  
    mA[4].v[0]=m15;
    mA[0].v[1]=_mm_setr_pd(mA[2].d[0],mA[3].d[0]);
    mA[1].v[1]=_mm_setr_pd(mA[2].d[1],mA[3].d[1]);  
    mA[2].v[1]=m23;
    mA[3].v[1]=m24;  
    mA[4].v[1]=m25;
    mA[0].v[2]=_mm_setr_pd(mA[4].d[0],0.);
    mA[1].v[2]=_mm_setr_pd(mA[4].d[1],0.);
    mA[2].v[2]=_mm_setr_pd(mA[4].d[2],0.);
    mA[3].v[2]=_mm_setr_pd(mA[4].d[3],0.);
    mA[4].v[2]=m35;
  }
  // Constructor for symmetric matrix by elements
  DMatrix5x5(double C11,double C12,double C13,double C14,double C15,double C22,double C23,double C24,double C25,
	     double C33,double C34,double C35,double C44,double C45,double C55){
    mA[0].v[0]=_mm_setr_pd(C11,C12);
    mA[1].v[0]=_mm_setr_pd(C12,C22);
    mA[2].v[0]=_mm_setr_pd(C13,C23);
    mA[3].v[0]=_mm_setr_pd(C14,C24);
    mA[4].v[0]=_mm_setr_pd(C15,C25);
    mA[0].v[1]=_mm_setr_pd(C13,C14);
    mA[1].v[1]=_mm_setr_pd(C23,C24);
    mA[2].v[1]=_mm_setr_pd(C33,C34);
    mA[3].v[1]=_mm_setr_pd(C34,C44);
    mA[4].v[1]=_mm_setr_pd(C35,C45);
    mA[0].v[2]=_mm_setr_pd(C15,0);
    mA[1].v[2]=_mm_setr_pd(C25,0);
    mA[2].v[2]=_mm_setr_pd(C35,0);
    mA[3].v[2]=_mm_setr_pd(C45,0);
    mA[4].v[2]=_mm_setr_pd(C55,0);
    
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
  // Method for adding symmetric matrices
  DMatrix5x5 AddSym(const DMatrix5x5 &m2) const{
    return DMatrix5x5(ADD(0,0),ADD(0,1),ADD(0,2),ADD(0,3),ADD(0,4),ADD(1,2),
		      ADD(1,3),ADD(1,4),ADD(2,4));
  }

  // Matrix subtraction
  DMatrix5x5 operator-(const DMatrix5x5 &m2) const{
#define SUB(i,j) _mm_sub_pd(GetV((i),(j)),m2.GetV((i),(j)))

    return DMatrix5x5(SUB(0,0),SUB(0,1),SUB(0,2),SUB(0,3),SUB(0,4),SUB(1,0),SUB(1,1),SUB(1,2),
		      SUB(1,3),SUB(1,4),SUB(2,0),SUB(2,1),SUB(2,2),SUB(2,3),SUB(2,4));
  }
  // method for subtracting a symmetric matrix from another symmetric matrix
  DMatrix5x5 SubSym(const DMatrix5x5 &m2) const{
    return DMatrix5x5(SUB(0,0),SUB(0,1),SUB(0,2),SUB(0,3),SUB(0,4),SUB(1,2),
		      SUB(1,3),SUB(1,4),SUB(2,4));
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

#ifdef USE_SSE3
  
  // The following code performs the matrix operation ABA^T, where B is a symmetric matrix
  DMatrix5x5 SandwichMultiply(const DMatrix5x5 &A){
    __m128d A1112=_mm_setr_pd(A(0,0),A(0,1));
    __m128d A1314=_mm_setr_pd(A(0,2),A(0,3));
    __m128d A15=_mm_set1_pd(A(0,4));
    union dvec temp;

    // BA^T column 1
    temp.v[0]=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(GetV(0,0),A1112),
						_mm_mul_pd(GetV(0,1),A1112)),
				    _mm_hadd_pd(_mm_mul_pd(GetV(1,0),A1314),
						_mm_mul_pd(GetV(1,1),A1314))),
			 _mm_mul_pd(GetV(0,4),A15));
    temp.v[1]=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(GetV(0,2),A1112),
						_mm_mul_pd(GetV(0,3),A1112)),
				    _mm_hadd_pd(_mm_mul_pd(GetV(1,2),A1314),
						_mm_mul_pd(GetV(1,3),A1314))),
			 _mm_mul_pd(GetV(1,4),A15));
    temp.v[2] =_mm_set1_pd(mA[4].d[0]*A(0,0)+mA[4].d[1]*A(0,1)+mA[4].d[2]*A(0,2)+mA[4].d[3]*A(0,3)+mA[4].d[4]*A(0,4));
    
    union dvec2{
      __m128d v;
      double d[2];
    }temp2;
    
    // C=ABA^T elements C11,C12
#define A2122 _mm_setr_pd(A(1,0),A(1,1))
#define A2324 _mm_setr_pd(A(1,2),A(1,3))
    temp2.v=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A1112,temp.v[0]),_mm_mul_pd(A2122,temp.v[0])),
				  _mm_hadd_pd(_mm_mul_pd(A1314,temp.v[1]),_mm_mul_pd(A2324,temp.v[1]))),
		       _mm_mul_pd(A.GetV(0,4),temp.v[2]));
    double C11=temp2.d[0];
    double C12=temp2.d[1];

    // C=ABA^T elements C13,C14
#define A3132 _mm_setr_pd(A(2,0),A(2,1))
#define A3334 _mm_setr_pd(A(2,2),A(2,3))
#define A4142 _mm_setr_pd(A(3,0),A(3,1))
#define A4344 _mm_setr_pd(A(3,2),A(3,3))
    temp2.v=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A3132,temp.v[0]),_mm_mul_pd(A4142,temp.v[0])),
				    _mm_hadd_pd(_mm_mul_pd(A3334,temp.v[1]),_mm_mul_pd(A4344,temp.v[1]))),
			 _mm_mul_pd(A.GetV(1,4),temp.v[2]));
    double C13=temp2.d[0];
    double C14=temp2.d[1];

    // BA^T column 2
#define A25 _mm_set1_pd(A(1,4))
    temp.v[0]=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(GetV(0,0),A2122),
						_mm_mul_pd(GetV(0,1),A2122)),
				    _mm_hadd_pd(_mm_mul_pd(GetV(1,0),A2324),
						_mm_mul_pd(GetV(1,1),A2324))),
			 _mm_mul_pd(GetV(0,4),A25));
    temp.v[1]=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(GetV(0,2),A2122),
						_mm_mul_pd(GetV(0,3),A2122)),
				    _mm_hadd_pd(_mm_mul_pd(GetV(1,2),A2324),
						_mm_mul_pd(GetV(1,3),A2324))),
			 _mm_mul_pd(GetV(1,4),A25));
    temp.v[2]=_mm_set1_pd(mA[4].d[0]*A(1,0)+mA[4].d[1]*A(1,1)+mA[4].d[2]*A(1,2)+mA[4].d[3]*A(1,3)+mA[4].d[4]*A(1,4));
    
    // C=ABA^T elements C22,C23
#define A2535 _mm_setr_pd(A(1,4),A(2,4))
    temp2.v=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A2122,temp.v[0]),_mm_mul_pd(A3132,temp.v[0])),
				    _mm_hadd_pd(_mm_mul_pd(A2324,temp.v[1]),_mm_mul_pd(A3334,temp.v[1]))),
			 _mm_mul_pd(A2535,temp.v[2]));
    double C22=temp2.d[0];
    double C23=temp2.d[1];

    // C=ABA^T elements C24,C25
#define A5152 _mm_setr_pd(A(4,0),A(4,1))
#define A5354 _mm_setr_pd(A(4,2),A(4,3)) 
#define A4555 _mm_setr_pd(A(3,4),A(4,4))
    temp2.v=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A4142,temp.v[0]),_mm_mul_pd(A5152,temp.v[0])),
				    _mm_hadd_pd(_mm_mul_pd(A4344,temp.v[1]),_mm_mul_pd(A5354,temp.v[1]))),
			 _mm_mul_pd(A4555,temp.v[2]));  
    double C24=temp2.d[0];
    double C25=temp2.d[1];

    // BA^T column 3
#define A35 _mm_set1_pd(A(2,4))
    temp.v[0]=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(GetV(0,0),A3132),
						_mm_mul_pd(GetV(0,1),A3132)),
				    _mm_hadd_pd(_mm_mul_pd(GetV(1,0),A3334),
						_mm_mul_pd(GetV(1,1),A3334))),
			 _mm_mul_pd(GetV(0,4),A35));
    temp.v[1]=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(GetV(0,2),A3132),
						_mm_mul_pd(GetV(0,3),A3132)),
				    _mm_hadd_pd(_mm_mul_pd(GetV(1,2),A3334),
						_mm_mul_pd(GetV(1,3),A3334))),
			 _mm_mul_pd(GetV(1,4),A35));
    temp.v[2]=_mm_set1_pd(mA[4].d[0]*A(2,0)+mA[4].d[1]*A(2,1)+mA[4].d[2]*A(2,2)+mA[4].d[3]*A(2,3)+mA[4].d[4]*A(2,4));
    

    // C=ABA^T elements C33,C34,C35
#define A3545 _mm_setr_pd(A(2,4),A(3,4))
    temp2.v=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A3132,temp.v[0]),_mm_mul_pd(A4142,temp.v[0])),
				    _mm_hadd_pd(_mm_mul_pd(A3334,temp.v[1]),_mm_mul_pd(A4344,temp.v[1]))),
			 _mm_mul_pd(A3545,temp.v[2]));
    double C33=temp2.d[0];
    double C34=temp2.d[1];
    double C35=A(4,0)*temp.d[0]+A(4,1)*temp.d[1]+A(4,2)*temp.d[2]+A(4,3)*temp.d[3]+A(4,4)*temp.d[4];

// BA^T column 4
#define A45 _mm_set1_pd(A(3,4))
    temp.v[0]=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(GetV(0,0),A4142),
						_mm_mul_pd(GetV(0,1),A4142)),
				    _mm_hadd_pd(_mm_mul_pd(GetV(1,0),A4344),
						_mm_mul_pd(GetV(1,1),A4344))),
			 _mm_mul_pd(GetV(0,4),A45));
    temp.v[1]=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(GetV(0,2),A4142),
						_mm_mul_pd(GetV(0,3),A4142)),
				    _mm_hadd_pd(_mm_mul_pd(GetV(1,2),A4344),
						_mm_mul_pd(GetV(1,3),A4344))),
			 _mm_mul_pd(GetV(1,4),A45));
    temp.v[2]=_mm_set1_pd(mA[4].d[0]*A(3,0)+mA[4].d[1]*A(3,1)+mA[4].d[2]*A(3,2)+mA[4].d[3]*A(3,3)+mA[4].d[4]*A(3,4));

    // C=ABA^T elements C44,C45
#define A4545 _mm_setr_pd(A(3,4),A(4,4))
    temp2.v=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A4142,temp.v[0]),_mm_mul_pd(A5152,temp.v[0])),
				    _mm_hadd_pd(_mm_mul_pd(A4344,temp.v[1]),_mm_mul_pd(A5354,temp.v[1]))),
			 _mm_mul_pd(A4555,temp.v[2]));
    double C44=temp2.d[0];
    double C45=temp2.d[1];
 
       // BA^T column 5
#define A55 _mm_set1_pd(A(4,4))
    temp.v[0]=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(GetV(0,0),A5152),
						_mm_mul_pd(GetV(0,1),A5152)),
				    _mm_hadd_pd(_mm_mul_pd(GetV(1,0),A5354),
						_mm_mul_pd(GetV(1,1),A5354))),
			 _mm_mul_pd(GetV(0,4),A55));
    temp.v[1]=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(GetV(0,2),A5152),
						_mm_mul_pd(GetV(0,3),A5152)),
				    _mm_hadd_pd(_mm_mul_pd(GetV(1,2),A5354),
						_mm_mul_pd(GetV(1,3),A5354))),
			 _mm_mul_pd(GetV(1,4),A55));
    temp.v[2]=_mm_set1_pd(mA[4].d[0]*A(4,0)+mA[4].d[1]*A(4,1)+mA[4].d[2]*A(4,2)+mA[4].d[3]*A(4,3)+mA[4].d[4]*A(4,4));


  // C=ABA^T elements C15,C55
#define A1555 _mm_setr_pd(A(0,4),A(4,4)) 
    temp2.v=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A1112,temp.v[0]),_mm_mul_pd(A5152,temp.v[0])),
				    _mm_hadd_pd(_mm_mul_pd(A1314,temp.v[1]),_mm_mul_pd(A5354,temp.v[1]))),
			 _mm_mul_pd(A1555,temp.v[2]));
    double C15=temp2.d[0];
    double C55=temp2.d[1];

    return DMatrix5x5(C11,C12,C13,C14,C15,C22,C23,C24,C25,C33,C34,C35,C44,C45,C55);
  }

  // Matrix multiplication. Requires the SSE3 instruction HADD (horizontal add)
  DMatrix5x5 operator*(const DMatrix5x5 &m2){
    __m128d A11A12=_mm_setr_pd(mA[0].d[0],mA[1].d[0]);
    __m128d A13A14=_mm_setr_pd(mA[2].d[0],mA[3].d[0]);
    __m128d A21A22=_mm_setr_pd(mA[0].d[1],mA[1].d[1]);
    __m128d A23A24=_mm_setr_pd(mA[2].d[1],mA[3].d[1]);
    __m128d A31A32=_mm_setr_pd(mA[0].d[2],mA[1].d[2]);
    __m128d A33A34=_mm_setr_pd(mA[2].d[2],mA[3].d[2]);
    __m128d A41A42=_mm_setr_pd(mA[0].d[3],mA[1].d[3]);
    __m128d A43A44=_mm_setr_pd(mA[2].d[3],mA[3].d[3]);

   return
     DMatrix5x5(_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A11A12,m2.GetV(0,0)),_mm_mul_pd(A21A22,m2.GetV(0,0))),
				      _mm_hadd_pd(_mm_mul_pd(A13A14,m2.GetV(1,0)),_mm_mul_pd(A23A24,m2.GetV(1,0)))),
			   _mm_mul_pd(mA[4].v[0],_mm_set1_pd(m2(4,0)))),
		_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A11A12,m2.GetV(0,1)),_mm_mul_pd(A21A22,m2.GetV(0,1))),
				      _mm_hadd_pd(_mm_mul_pd(A13A14,m2.GetV(1,1)),_mm_mul_pd(A23A24,m2.GetV(1,1)))),
			   _mm_mul_pd(mA[4].v[0],_mm_set1_pd(m2(4,1)))),
		_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A11A12,m2.GetV(0,2)),_mm_mul_pd(A21A22,m2.GetV(0,2))),
				      _mm_hadd_pd(_mm_mul_pd(A13A14,m2.GetV(1,2)),_mm_mul_pd(A23A24,m2.GetV(1,2)))),
			   _mm_mul_pd(mA[4].v[0],_mm_set1_pd(m2(4,2)))),
		_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A11A12,m2.GetV(0,3)),_mm_mul_pd(A21A22,m2.GetV(0,3))),
				      _mm_hadd_pd(_mm_mul_pd(A13A14,m2.GetV(1,3)),_mm_mul_pd(A23A24,m2.GetV(1,3)))),
			   _mm_mul_pd(mA[4].v[0],_mm_set1_pd(m2(4,3)))),
		_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A11A12,m2.GetV(0,4)),_mm_mul_pd(A21A22,m2.GetV(0,4))),
				      _mm_hadd_pd(_mm_mul_pd(A13A14,m2.GetV(1,4)),_mm_mul_pd(A23A24,m2.GetV(1,4)))),
			   _mm_mul_pd(mA[4].v[0],_mm_set1_pd(m2(4,4)))),
		
		_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A31A32,m2.GetV(0,0)),_mm_mul_pd(A41A42,m2.GetV(0,0))),
				      _mm_hadd_pd(_mm_mul_pd(A33A34,m2.GetV(1,0)),_mm_mul_pd(A43A44,m2.GetV(1,0)))),
			   _mm_mul_pd(mA[4].v[1],_mm_set1_pd(m2(4,0)))),
		_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A31A32,m2.GetV(0,1)),_mm_mul_pd(A41A42,m2.GetV(0,1))),
				      _mm_hadd_pd(_mm_mul_pd(A33A34,m2.GetV(1,1)),_mm_mul_pd(A43A44,m2.GetV(1,1)))),
			   _mm_mul_pd(mA[4].v[1],_mm_set1_pd(m2(4,1)))),	
		_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A31A32,m2.GetV(0,2)),_mm_mul_pd(A41A42,m2.GetV(0,2))),
				      _mm_hadd_pd(_mm_mul_pd(A33A34,m2.GetV(1,2)),_mm_mul_pd(A43A44,m2.GetV(1,2)))),
			   _mm_mul_pd(mA[4].v[1],_mm_set1_pd(m2(4,2)))),
		_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A31A32,m2.GetV(0,3)),_mm_mul_pd(A41A42,m2.GetV(0,3))),
				      _mm_hadd_pd(_mm_mul_pd(A33A34,m2.GetV(1,3)),_mm_mul_pd(A43A44,m2.GetV(1,3)))),
			   _mm_mul_pd(mA[4].v[1],_mm_set1_pd(m2(4,3)))),	
		_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A31A32,m2.GetV(0,4)),_mm_mul_pd(A41A42,m2.GetV(0,4))),
				      _mm_hadd_pd(_mm_mul_pd(A33A34,m2.GetV(1,4)),_mm_mul_pd(A43A44,m2.GetV(1,4)))),
			   _mm_mul_pd(mA[4].v[1],_mm_set1_pd(m2(4,4)))),	
		
		_mm_set_sd(mA[0].d[4]*m2(0,0)+mA[1].d[4]*m2(1,0)+mA[2].d[4]*m2(2,0)+mA[3].d[4]*m2(3,0)+mA[4].d[4]*m2(4,0)),
		_mm_set_sd(mA[0].d[4]*m2(0,1)+mA[1].d[4]*m2(1,1)+mA[2].d[4]*m2(2,1)+mA[3].d[4]*m2(3,1)+mA[4].d[4]*m2(4,1)),
		_mm_set_sd(mA[0].d[4]*m2(0,2)+mA[1].d[4]*m2(1,2)+mA[2].d[4]*m2(2,2)+mA[3].d[4]*m2(3,2)+mA[4].d[4]*m2(4,2)),
		_mm_set_sd(mA[0].d[4]*m2(0,3)+mA[1].d[4]*m2(1,3)+mA[2].d[4]*m2(2,3)+mA[3].d[4]*m2(3,3)+mA[4].d[4]*m2(4,3)),
		_mm_set_sd(mA[0].d[4]*m2(0,4)+mA[1].d[4]*m2(1,4)+mA[2].d[4]*m2(2,4)+mA[3].d[4]*m2(3,4)+mA[4].d[4]*m2(4,4))
		);
 }


#else 
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
#endif


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
  // Matrix inversion by blocks for a symmetric matrix
  DMatrix5x5 InvertSym(){
    DMatrix2x2 A(GetV(0,0),GetV(0,1));
    DMatrix3x2 C(GetV(1,0),GetV(1,1),GetV(2,0),GetV(2,1));
    DMatrix3x2 CAinv=C*A.Invert();
    DMatrix3x3 D(GetV(1,2),GetV(1,3),GetV(1,4),GetV(2,2),GetV(2,3),GetV(2,4));
    DMatrix2x3 B(GetV(0,2),GetV(0,3),GetV(0,4));
    DMatrix2x3 BDinv=B*D.InvertSym();
    DMatrix2x2 AA=(A-BDinv*C).Invert();
    DMatrix3x3 DD=(D-CAinv*B).InvertSym();
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



