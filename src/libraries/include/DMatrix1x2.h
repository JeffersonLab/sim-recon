class DMatrix1x2{
 public:
  DMatrix1x2(){
    mA[0]=0.;
    mA[1]=0.;
  }
  DMatrix1x2(double c1,double c2){
    mA[0]=c1;
    mA[1]=c2;
  }
  ~DMatrix1x2(){};


  double &operator() (int col){
    return mA[col];
  } 
  double operator() (int col) const{
    return mA[col];
  }

  // Matrix multiplication:  (1x2) x (2x1)
  double operator*(const DMatrix2x1 &m2){
    return (mA[0]*m2(0)+mA[1]*m2(1));
  }

  // Matrix multiplication:  (1x2) x (2x2)
  DMatrix1x2 operator*(const DMatrix2x2 &m2){
    return 
      DMatrix1x2(mA[0]*m2(0,0)+mA[1]*m2(1,0),
		 mA[0]*m2(0,1)+mA[1]*m2(1,1));
  }



  
  void Print(){
    cout << "DMatrix1x2:" <<endl;
    cout << "     |      0    |      1    |" <<endl;
    cout << "------------------------------" <<endl;
    cout << "   0 |" ;
    for (unsigned int j=0;j<2;j++){
      cout <<setw(11)<<setprecision(4)<< mA[j]<<" "; 
    } 
    cout << endl;
  }      


 private:
  double mA[2];

};


