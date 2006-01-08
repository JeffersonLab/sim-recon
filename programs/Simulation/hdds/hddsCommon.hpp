/*  HDDS Common Classes
 *
 *  Author: richard.t.jones@uconn.edu
 *
 *  Original version - Richard Jones, January 6, 2006.
 *
 */

#ifndef SAW_HDDSCOMMON_DEF
#define SAW_HDDSCOMMON_DEF true

class Refsys
{
 /* The Refsys class is used to propagate coordinate system information
  * such as origin, orientation, and magnetic fields from parent volume
  * to daughter volume within the geomtry hierarchy.
  */
 public:
   DOMElement* fMother;		// current mother volume element
   int fMagField;		// flag indicating magnetic field
   double fOrigin[3];		// x,y,z coordinate of volume origin (cm)
   double fPhiOffset;		// azimuthal angle of volume origin (deg)
   double fRmatrix[3][3];	// rotation matrix (daughter -> mother)
   int fRotation;		// unique Rmatrix flag

   static int fVolumes;		// total number of volumes to far
   static XString fIdentifierList; // list of identifier strings (space-sep)

   struct VolIdent
   {
      XString fieldS;	
      int value;
      int step;
   };
   std::list<struct VolIdent> fIdentifier;	// identifier tag list 

   struct Partition
   {
      DOMElement* divEl;
      int ncopy;
      int iaxis;
      double start;
      double step;
   };
   struct Partition fPartition;		// prescription for axis partitioning

   Refsys();				// empty constructor
   Refsys(const Refsys& src);		// copy constructor
   Refsys& operator=(Refsys& src);	// copy operator
   Refsys& reset();			// reset origin, Rmatrix
   Refsys& reset(const Refsys& ref);	// reset origin, Rmatrix to ref
   Refsys& shift(const double vector[3]); // translate origin
   Refsys& shift(const Refsys& ref);	  // copy origin from ref
   Refsys& shift(const Refsys& ref,
                 const double vector[3]); // translate origin in ref frame
   Refsys& rotate(const double omega[3]); // rotate by vector omega (radians)
   Refsys& rotate(const Refsys& ref);	  // copy Rmatrix from ref
   Refsys& rotate(const Refsys& ref,
                  const double omega[3]); // rotate by omega in ref frame

   static int fRotations;	// non-trivial rotations defined so far
};

class Units
{
 /* The Units class provides conversion constants for creating readable
  * units-aware code.  For example, to ensure that the quantity "velocity"
  * is provided in units of km/hr regardless of the internal representation
  * simply write:
  *	distance_in_km_per_hr = distance * unit.km/unit.hr;
  * where unit is a Units class instance that has been previously set.
  * The user of the Units class should generally treat its data members
  * as constants and use Units::getConversions() to manage the values. 
  */
 public:
   double s,ns,ms;
   double min,hr,days,weeks;
   double m,km,cm,mm,um,nm;
   double in,ft,miles,mils;
   double rad,mrad,urad;
   double deg,arcmin,arcsec;
   double eV,KeV,MeV,GeV,TeV;
   double g,kg,mg;
   double m2,cm2,mm2;
   double b,mb,ub,nb,pb;
   double l,ml,cm3;
   double G,kG,Tesla;
   double percent;

   Units();				// empty constructor
   Units(Units& u);			// copy constructor
   void getConversions(DOMElement* el);	// get conversion constants from tag

 private:
   void set_1s(double tu);
   void set_1cm(double lu);
   void set_1rad(double au);
   void set_1deg(double au);
   void set_1KeV(double eu);
   void set_1MeV(double eu);
   void set_1GeV(double eu);
   void set_1g(double mu);
   void set_1cm2(double l2u);
   void set_1cm3(double l3u);
   void set_1G(double bfu);
};

class Substance
{
 /* The Substance class is used to collect and manage materials
  * property information for the simulation components.
  */
 public:
   Substance();
   Substance(DOMElement* elem);

   XString getName();		// return name of material
   XString getSymbol();		// return chem. symbol (if any)
   double getAtomicWeight();	// return A for a material
   double getAtomicNumber();	// return Z for a material
   double getDensity();		// return density [g/cm^3]
   double getRadLength();	// return radiation len. [cm]
   double getAbsLength();	// return nucl.abs.len. [cm]
   double getColLength();	// return nucl.col.len. [cm]
   double getMIdEdx();		// return min. dE/dx [MeV/g/cm^3]
   DOMElement* getDOMElement();	// return DOM element ptr

   int fUniqueID;		// user-assignable index
   class Brew
   {
    public:
      int natoms;		// number of atoms in chemical formula
      double wfact;		// fraction by weight in mixture
      Substance* sub;		// ptr to material description
   };
   std::list<Brew> fBrewList;

 protected:
   DOMElement* fMaterialEl;
   double fAtomicWeight;
   double fAtomicNumber;
   double fDensity;
   double fRadLen;
   double fAbsLen;
   double fColLen;
   double fMIdEdx;
};

class CodeWriter
{
 /* The CodeWriter class provides basic functionality for instantiating
  * the materials, solids and volume hierarchy that are derived from
  * the xml geometry specification.  An application for writing code
  * for a particular language/simulation system should extend this class
  * and add the code that does the actual output in methods that override
  * the createXXX() methods below.  The user methods probably should
  * invoke the corresponding methods of the base class in order to
  * obtain the correct geometry construction, or else they need to 
  * incorporate equivalent functionality themselves.
  */
 public:
   CodeWriter() {};
   void translate(DOMElement* el);		// invokes the code writer
   virtual void createHeader();
   virtual void createTrailer();
   virtual int createMaterial(DOMElement* el);	// generate code for materials
   virtual int createSolid(DOMElement* el,
                           const Refsys& ref);	// generate code for solids
   virtual int createRotation(Refsys& ref);	// generate code for rotations
   virtual int createVolume(DOMElement* el,
                            const Refsys& ref);	// generate code for placement
   virtual int createDivision(XString& divStr,
                              Refsys& ref);	// generate code for divisions
   virtual void createGetFunctions(DOMElement* el,
                              XString& ident);	// generate code for identifiers

 protected:
   bool fPending;       // indicates a volume positioning request is pending
   Substance fSubst;    // work area for latest material definition
   Refsys fRef;		// work area for latest reference system
};

#endif
