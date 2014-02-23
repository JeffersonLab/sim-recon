class DMatrix1x4{
 public:
  DMatrix1x4(){
    for (unsigned int i=0;i<4;i++){
      mA[i]=0.;
    }
  }
  DMatrix1x4(double c1,double c2,double c3,double c4){
    mA[0]=c1;
    mA[1]=c2;
    mA[2]=c3;
    mA[3]=c4;
  }
  ~DMatrix1x4(){};


  double &operator() (int col){
    return mA[col];
  } 
  double operator() (int col) const{
    return mA[col];
  }

  // Matrix multiplication:  (1x4) x (4x1)
  double operator*(const DMatrix4x1 &m2){
    return (mA[0]*m2(0)+mA[1]*m2(1)+mA[2]*m2(2)+mA[3]*m2(3));
  }

  // Matrix multiplication:  (1x4) x (4x4)
  DMatrix1x4 operator*(const DMatrix4x4 &m2){
    return 
      DMatrix1x4(mA[0]*m2(0,0)+mA[1]*m2(1,0)+mA[2]*m2(2,0)+mA[3]*m2(3,0),
		 mA[0]*m2(0,1)+mA[1]*m2(1,1)+mA[2]*m2(2,1)+mA[3]*m2(3,1),
		 mA[0]*m2(0,2)+mA[1]*m2(1,2)+mA[2]*m2(2,2)+mA[3]*m2(3,2),
		 mA[0]*m2(0,3)+mA[1]*m2(1,3)+mA[2]*m2(2,3)+mA[3]*m2(3,3));
  }



  
  void Print(){
    cout << "DMatrix1x4:" <<endl;
    cout << "     |      0    |      1    |      2    |      3    |" <<endl;
    cout << "------------------------------------------------------" <<endl;
    cout << "   0 |" ;
    for (unsigned int j=0;j<4;j++){
      cout <<setw(11)<<setprecision(4)<< mA[j]<<" "; 
    } 
    cout << endl;
  }      


 private:
  double mA[4];

};


