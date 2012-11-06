class DMatrix1x3{
 public:
  DMatrix1x3(){
    for (unsigned int i=0;i<3;i++){
      mA[i]=0.;
    }
  }
  DMatrix1x3(double c1,double c2,double c3){
    mA[0]=c1;
    mA[1]=c2;
    mA[2]=c3;
  }
  ~DMatrix1x3(){};


  double &operator() (int col){
    return mA[col];
  } 
  double operator() (int col) const{
    return mA[col];
  }

  // Matrix multiplication:  (1x3) x (3x1)
  double operator*(const DMatrix3x1 &m2){
    return (mA[0]*m2(0)+mA[1]*m2(1)+mA[2]*m2(2));
  }

  // Matrix multiplication:  (1x3) x (3x3)
  DMatrix1x3 operator*(const DMatrix3x3 &m2){
    return 
      DMatrix1x3(mA[0]*m2(0,0)+mA[1]*m2(1,0)+mA[2]*m2(2,0),
		 mA[0]*m2(0,1)+mA[1]*m2(1,1)+mA[2]*m2(2,1),
		 mA[0]*m2(0,2)+mA[1]*m2(1,2)+mA[2]*m2(2,2));
  }



  
  void Print(){
    cout << "DMatrix1x3:" <<endl;
    cout << "     |      0    |      1    |      2     |" <<endl;
    cout << "------------------------------------------------------" <<endl;
    cout << "   0 |" ;
    for (unsigned int j=0;j<3;j++){
      cout <<setw(11)<<setprecision(4)<< mA[j]<<" "; 
    } 
    cout << endl;
  }      


 private:
  double mA[3];

};


