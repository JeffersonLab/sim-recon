class pseudoResidFunc:public residFunc {
 public:
  pseudoResidFunc(vector<const DFDCPseudo*> *pseudopoints,
				 MyTrajectory *trajectory);
  void resid(const HepVector *x, void *data, HepVector *f);
  void deriv(const HepVector *x, void *data, HepMatrix *J);
  void residAndDeriv(const HepVector *x, void *data, HepVector *f,
		     HepMatrix *J);
  inline unsigned int getn() {return n;};
  inline unsigned int getp() {return trajPtr->getNumberOfParams();};
 private:
  unsigned int n;
  vector<const DFDCPseudo*> *ppPtr;
  MyTrajectory *trajPtr;
  vector<double> delta;
};
