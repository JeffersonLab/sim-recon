class parabolaResidFunc:public residFunc {
 public:
  parabolaResidFunc(unsigned int nIn): n(nIn) {};
  void resid(const HepVector *x, void *data, HepVector *f);
  void deriv(const HepVector *x, void *data, HepMatrix *J);
  void residAndDeriv(const HepVector *x, void *data, HepVector *f,
		     HepMatrix *J);
  unsigned int getn() {return n;};
 private:
  unsigned int n;
};
