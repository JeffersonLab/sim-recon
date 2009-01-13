// Class for dealing with a map of the material in the detector

#ifndef _DMaterialMap_
#define _DMaterialMap_

#define NUM_Z_POINTS 1400
#define NUM_X_POINTS 260

class DMaterialMap{
 public:
  DMaterialMap(){};
  virtual ~DMaterialMap(){};

  double GetRadLen(double x, double y, double z) const;

 protected:
  double material_z[NUM_Z_POINTS];
  double material_x[NUM_X_POINTS];
  double radlen[NUM_Z_POINTS][NUM_X_POINTS];

};

#endif // _DMaterialMap_
