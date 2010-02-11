/*  HDDS Common Classes
 *
 *  Author: richard.t.jones@uconn.edu
 *
 *  Original version - Richard Jones, January 6, 2006.
 *
 */

#ifndef SAW_HDDSCOMMON_DEF
#define SAW_HDDSCOMMON_DEF true

#include <vector>
#include <list>
#include <map>
#include "XString.hpp"

class Refsys
{
 /* The Refsys class is used to propagate coordinate system information
  * such as origin, orientation, and magnetic fields from parent volume
  * to daughter volume within the geomtry hierarchy.
  *
  * The following conventions are used by hdds in describing transformations
  * between mother and daughter coordinate systems.  Let x be a vector
  * describing a point in reference system M, and x' be the coordinates of
  * the same point in daughter reference system D.
  *
  *        x = R x' + s ,  where R = Rz Ry Rx, and s is a shift of origin.
  *
  *              / 1         0             0        \
  *        Rx = |  0    cos(omega_x)   sin(omega_x)  |
  *              \ 0   -sin(omega_x)   cos(omega_x) /
  *
  *              /-sin(omega_y)   0    cos(omega_y) \
  *        Ry = |      0          1        0         |
  *              \ cos(omega_y)   0    sin(omega_y) /
  *
  *              / cos(omega_z)   sin(omega_z)   0  \
  *        Rz = | -sin(omega_z)   cos(omega_z)   0   |
  *              \     0               0         1  /
  *
  * This corresponds to transformation from daughter to mother coordinates
  * after placement of the daughter in the mother using an active rotation
  * by omega_x around Xhat, by omega_y around Yhat' and by omega_z around
  * Zhat'' performed in that order, followed by a shift of the origin by
  * vector s.  The three angles are specified in the hdds positioning
  * elements using the attribute rot="omega_x omega_y omega_z" and the
  * offset s is specified in various ways, depending on the positioning
  * tag.  For example, in the posXYZ tag the value of s is specified in
  * terms of its three Cartesian coordinates as X_Y_Z="s_x s_y s_z".
  *
  * The Refsys object carries the notion of a current mother volume and
  * daughter volume as a part of its state.  The transformation is stored
  * internally in the form of the 3x3 real orthogonal matrix R (see above)
  * denoted fRmatrix and the 3-vector s denoted fOrigin.  Successive steps
  * in the placement heirarchy are combined into one effective step using
  * the rules of linear algebra as follows.  To place a daughter inside
  * the existing daughter with a rotation R' followed by a shift s',
  *
  *         R --> R (Rz' Ry' Rx')   and  s --> s + R s'
  *
  * These actions are applied using the methods shift() and rotate().
  * Simulation geometry interfaces generally need to know only enough
  * information to place the daughter inside its immediate mother, which
  * means that before each placement step the transformation must be reset
  * back to R=1, s=0.  This functionality is provided by the reset() method.
  * There are also cases where the code generator would like to know the
  * accumulated transformation from the master reference system to that of
  * the current daughter, for example in the case of looking up a point in
  * a magnetic field map.  For that purpose, the Refsys object state also
  * stores a second set of transformation arrays, fMRmatrix and fMOrigin.
  * These are initialized by the copy constructor and are updated by the
  * rotate() and shift() methods, but are not reset by the reset() method.
  * The master reference system (MRS) coincides with the local coordinate
  * system of the root volume in the geometry tree.
  */
 public:
   DOMElement* fMother;        	// current mother volume element
   DOMElement* fRegion;        	// associated region, 0 if default
   double fPhiOffset;        	// azimuthal angle of volume origin (deg)
   double fOrigin[3];        	// x,y,z coordinate of volume origin (cm)
   double fRmatrix[3][3];       // rotation matrix (daughter -> mother)
   int fRotation;        	// unique Rmatrix identifier
   int fRegionID;        	// unique region identifier

   double fMOrigin[3];        	// same as fOrigin, but with MRS as mother
   double fMRmatrix[3][3];      // same as fRmatrix, but with MRS as mother

   struct VolIdent
   {
      int value;
      int step;
   };
   std::map<std::string,VolIdent> fIdentifier;          // identifier list 
   static std::map<std::string,VolIdent> fIdentifiers;  // master id list 

   static std::vector<std::map<std::string,std::vector<int> > >
          fIdentifierTable;			    // identifier lookup maps

   struct Partition
   {
      DOMElement* divEl;
      int ncopy;
      int iaxis;
      double start;
      double step;
   };
   Partition fPartition;		// prescription for axis partitioning

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

   void addIdentifier(XString ident,
                      int value,
                      int step);// adds a new identifier to the current list
   void incrementIdentifiers(); // increments all identifiers by one step
   void clearIdentifiers();	// resets the current identifier list

   int nextRotationID();	// generate unique rotation index sequence
   int nextRegionID();		// generate unique region index sequence
   int nextVolumeID();		// generate unique volume index sequence

   std::map<std::string,double> fPar; // key-value table for user needs

   static int fRotations;	// non-trivial rotations defined so far
   static int fRegions;        	// total number of regions so far
   static int fVolumes;        	// total number of volumes so far
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
   Substance(Substance& src);
   Substance(DOMElement* elem);
   ~Substance();

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

   Substance& operator=(const Substance& src);

   int fUniqueID;		// user-assignable index
   class Brew
   {
    public:
      int natoms;		// number of atoms in chemical formula
      double wfact;		// fraction by weight in mixture
      Substance* sub;		// ptr to material description of component
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
                           Refsys& ref);	// generate code for solids
   virtual int createRotation(Refsys& ref);	// generate code for rotations
   virtual int createRegion(DOMElement* el,
                            Refsys& ref);	// generate code for regions
   virtual int createVolume(DOMElement* el,
                            Refsys& ref);	// generate code for placement
   virtual int createDivision(XString& divStr,
                              Refsys& ref);	// generate code for divisions
   virtual void createSetFunctions(DOMElement* el,
                         const XString& ident);	// generate property setters
   virtual void createGetFunctions(DOMElement* el,
                         const XString& ident);	// generate identifier lookups
   virtual void createMapFunctions(DOMElement* el,
                         const XString& ident);	// generate field map functions
   virtual void createUtilityFunctions(DOMElement* el,
			 const XString& ident);	// generate utility functions

 protected:
   bool fPending;       // indicates a volume positioning request is pending
   Substance fSubst;    // work area for latest material definition
   Refsys fRef;		// work area for latest reference system

 private:
   void dump(DOMElement* el, int level);  // useful for debugging, keep me!
};

# ifdef LINUX_CPUTIME_PROFILING
#   include <sys/resource.h>
#   include <sys/time.h>
#   include <unistd.h>

class CPUtimer
{
 /* The CPUtimer class is useful for programs that want to monitor their
  * usage of cpu time.  Under unix there are three ways of measuring the
  * cpu time on a machine:
  *
  *	user time   : the time a cpu was running the user's code.
  *     system time : the time spent in the kernel (eg. doing i/o) on
  *                   behalf of the user's code (other than sleeping).
  *     real time   : the time passed on the wall clock since the user's
  *                   code began execution.
  *
  * All three of these measures are available through calls to the system
  * function getrusage() on a linux system.  This class may need a new
  * implementation to port it to other non-linux systems.
  */

 public:
   CPUtimer();

   double getUserTime();	// returns total user time (s)
   double getSystemTime();	// returns total system time (s)
   double getRealTime();	// returns total real time in (s)
   void resetClocks();		// resets the time references to zero
   double getUserDelta();	// returns total user time since last reset
   double getSystemDelta();	// returns total system time since last reset
   double getRealDelta();	// returns total real time since last reset
  
 private:
   struct rusage fRef;
   struct rusage fLast;
   struct timezone fTZ;
   struct timeval fClock;
   struct timeval fClock0;
   struct timeval fClockRef;

   int getRusage();
};
# endif

#endif
