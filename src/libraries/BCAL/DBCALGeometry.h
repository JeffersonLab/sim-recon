// $Id$
//
//    File: DBCALGeometry.h
// Created: Thu Nov 17 15:10:51 CST 2005
// Creator: gluexuser (on Linux hydra.phys.uregina.ca 2.4.20-8smp i686)
//

#ifndef _DBCALGeometry_
#define _DBCALGeometry_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

// create a single number channel id which is useful in algorithms
// if M L S are module layer sector the bit map looks like:
// MMMM MMMM LLLL SSSS

#define MODULE_SHIFT 8
#define MODULE_MASK 0xFF00
#define LAYER_SHIFT 4
#define LAYER_MASK 0x00F0
#define SECTOR_SHIFT 4
#define SECTOR_MASK 0X000F

class DBCALGeometry : public JObject{
    
public:
    
    JOBJECT_PUBLIC(DBCALGeometry);
    
    enum End { kUpstream, kDownstream };
    
    int NBCALMODS;         ///> number of modules
    int NBCALLAYS1;        ///> number of layers in first 10 ccm 
    int NBCALLAYS2;        ///> number of layers in last  15 cm 
    int NBCALSECS1;        ///> number of sectors in first 10cm of Mod 
    int NBCALSECS2;        ///> number of sectors in last 15cm of Mod 
    float BCALINNERRAD;    ///> innner radius of BCAL in cm
    float BCALMIDRAD;      ///> mid radius of BCAL in cm
    float BCALOUTERRAD;    ///> outer radius of BCAL in cm
    float BCALFIBERLENGTH; ///> BCAL Scintilator fiber lenth in cm
    float GLOBAL_CENTER;  ///> center of BCAL in gloobal coordinate system
    
    float ATTEN_LENGTH;   ///> attenuation length
    float C_EFFECTIVE;    ///> speed of light in fibers

    //Public Data Members for smear.cc
    static int NBcalMods() { return 48; }  
    static int NBcalLays1() { return 2; }
    static int NBcalLays2() { return 2; }
    static int NBcalSecs1() { return 4; }
    static int NBcalSecs2() { return 2; }
    static float BcalInnerRad() { return 64.3; }
    static float BcalMidRad() { return 76.3; }
    static float BcalOuterRad() { return 86.17; }
    static float BcalFiberLength() { return 390.0; }
    static float Global_Center() { return 212; }
    static float Atten_Length() { return 300.; }
    static float C_Effective() { return 16.75;}
    
    static inline int module( int cellId ) { return ( cellId & MODULE_MASK ) >> MODULE_SHIFT; };
    static inline int layer( int cellId ) { return ( cellId & LAYER_MASK ) >> LAYER_SHIFT; };
    static inline int sector( int cellId ) { return ( cellId & SECTOR_MASK ) >> SECTOR_SHIFT; };                
    
    static inline int cellId( int module, int layer, int sector ) {
        return ( ( module << MODULE_SHIFT ) | ( layer << LAYER_SHIFT ) | 
                 ( sector << SECTOR_SHIFT ) ); }

	void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "mod", "%d", NBCALMODS);
			AddString(items, "layn1", "%d", NBCALLAYS1);
			AddString(items, "layn2", "%d", NBCALLAYS2);
			AddString(items, "secn1", "%d", NBCALSECS1);
			AddString(items, "secn2", "%d", NBCALSECS2);
			AddString(items, "inr", "%6.3f", BCALINNERRAD);
			AddString(items, "midr", "%6.3f", BCALMIDRAD);
			AddString(items, "otr", "%6.3f", BCALOUTERRAD);
			AddString(items, "length", "%6.3f", BCALFIBERLENGTH);
			AddString(items, "cntr", "%3.2g", GLOBAL_CENTER);
			AddString(items, "atten len", "%3.2g", ATTEN_LENGTH);
			AddString(items, "c eff.", "%3.2g", C_EFFECTIVE);
	}
};

#endif // _DBCALGeometry_
