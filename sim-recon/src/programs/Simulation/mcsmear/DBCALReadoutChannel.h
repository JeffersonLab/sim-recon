// $Id$
//
//    File: DBCALReadoutChannel.h
// Created: Tue Apr 12 14:18:19 EDT 2011
// Creator: davidl (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DBCALReadoutChannel_
#define _DBCALReadoutChannel_

/// One object of this class is created for each BCAL electronic
/// channel. The objects are created once at program start up and
/// certain values filled in that exist for the life of the object.
///
/// The uphits and downhits containers are cleared at the start of 
/// each event and filled as the s_BcalSiPMUpHit_t and s_BcalSiPMDownHit_t
/// structs are created in the HDDM tree.
///
/// The s_BcalfADCUpHit_t and s_BcalfADCDownHit_t structs are made
/// later using the uphits and downhits. The value of NSiPM
/// is used to hold how many SiPMs are being added together for
/// this readout channel. The number of channels to generate dark
/// hits for the readout channel is calculated from:
///
/// Ndark_channels_up   = NSiPM - uphits.size();
/// Ndark_channels_down = NSiPM - downhits.size();
///
/// This is susceptable to a problem when more than one s_BcalSiPMUpHit_t
/// is generated for a given cell and they are not properly separated 
/// into multiple s_BcalfADCUpHit_t structs.
///
/// The value for <i>threshold</i> is calculated once at program startup
/// based on whether this is an "inner" or an "outer" cell. It is applied
/// to all summed energy values before creating a s_BcalfADCUpHit_t struct
/// (or a s_BcalfADCDownHit_t).

#include <vector>

#include <HDDM/hddm_s.h>

class DBCALReadoutChannel{
	public:
		DBCALReadoutChannel():NSiPM(0),threshold(0){}
		DBCALReadoutChannel(unsigned int NSiPM, double threshold, int module, int layer, int sector):NSiPM(NSiPM),threshold(threshold),module(module),layer(layer),sector(sector){}
		virtual ~DBCALReadoutChannel(){}
		
		// Assigned at object creation and never changed
		unsigned int NSiPM;
		double threshold;
		int module;
		int layer;
		int sector;
		
		// Assigned during event
		double Eup, Edown;
		double tup, tdown;
		
		
		// Cleared at the start of each event.
		std::vector<s_BcalSiPMUpHit_t*> uphits;
		std::vector<s_BcalSiPMDownHit_t*> downhits;
		
		void Clear(void){
			uphits.clear();
			downhits.clear();
			Eup = Edown = 0.0;
			tup = tdown = 0.0;
		}

	private:
};

#endif // _DBCALReadoutChannel_

