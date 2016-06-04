// $Id$
// $HeadURL$
//
//    File: Df125Config.h
// Created: Sun Sep  7 15:52:50 EDT 2014
// Creator: davidl (on Darwin harriet.local 13.3.0 i386)
//

#ifndef _Df125Config_
#define _Df125Config_

#include <DAQ/DDAQConfig.h>

class Df125Config:public DDAQConfig{
	public:
		JOBJECT_PUBLIC(Df125Config);

		Df125Config(uint32_t rocid, uint32_t slot_mask):DDAQConfig(rocid,slot_mask)
			,NSA(0xFFFF),NSB(0xFFFF),NSA_NSB(0xFFFF),NPED(0xFFFF),WINWIDTH(0xFFFF)
			,PL(0xFFFF),NW(0xFFFF),NPK(0xFFFF),P1(0xFFFF),P2(0xFFFF),PG(0xFFFF)
			,IE(0xFFFF),H(0xFFFF),TH(0xFFFF),TL(0xFFFF),IBIT(0xFFFF),ABIT(0xFFFF),PBIT(0xFFFF){}
		Df125Config(const Df125Config *c):DDAQConfig(c->rocid,c->slot_mask)
			,NSA(c->NSA),NSB(c->NSB),NSA_NSB(c->NSA_NSB),NPED(c->NPED),WINWIDTH(c->WINWIDTH)
			,PL(c->PL),NW(c->NW),NPK(c->NPK),P1(c->P1),P2(c->P2),PG(c->PG)
			,IE(c->IE),H(c->H),TH(c->TH),TL(c->TL),IBIT(c->IBIT),ABIT(c->ABIT),PBIT(c->PBIT){}
		
		uint16_t NSA;      // Num. samples before threshold crossing sample
		uint16_t NSB;      // Num. samples after  threshold crossing sample
		uint16_t NSA_NSB;  // NSA+NSB = total number of samples in integration window
		uint16_t NPED;     // Number of samples used to determine pedestal
		uint16_t WINWIDTH; // maximum integration window size (in samples)
		
		// See GlueX-doc-2274
		uint16_t PL;
		uint16_t NW;
		uint16_t NPK;
		uint16_t P1;
		uint16_t P2;
		uint16_t PG;
		uint16_t IE;
		uint16_t H;
		uint16_t TH;
		uint16_t TL;
		uint16_t IBIT;
		uint16_t ABIT;
		uint16_t PBIT;
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			DDAQConfig::toStrings(items);
			AddString(items, "NSA"      , "%d", NSA);
			AddString(items, "NSB"      , "%d", NSB);
			AddString(items, "NSA_NSB"  , "%d", NSA_NSB);
			AddString(items, "NPED"     , "%d", NPED);
			AddString(items, "WINWIDTH" , "%d", WINWIDTH);
			AddString(items, "PL"       , "%d", PL);
			AddString(items, "NW"       , "%d", NW);
			AddString(items, "NPK"      , "%d", NPK);
			AddString(items, "P1"       , "%d", P1);
			AddString(items, "P2"       , "%d", P2);
			AddString(items, "PG"       , "%d", PG);
			AddString(items, "IE"       , "%d", IE);
			AddString(items, "H"        , "%d", H);
			AddString(items, "TH"       , "%d", TH);
			AddString(items, "TL"       , "%d", TL);
			AddString(items, "IBIT"     , "%d", IBIT);
			AddString(items, "ABIT"     , "%d", ABIT);
			AddString(items, "PBIT"     , "%d", PBIT);
		}
		
};

#endif // _Df125Config_

