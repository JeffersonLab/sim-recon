// $Id$
//
//    File: DPSHit_factory.cc
// Created: Wed Oct 15 16:45:01 EDT 2014
// Creator: staylor (on Linux gluon05.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//


#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
using namespace std;

#include "DPSHit_factory.h"
#include "DPSDigiHit.h"
#include "DPSHit.h"
#include <DAQ/Df250PulseIntegral.h>
#include <DAQ/Df250Config.h>
using namespace jana;


//------------------
// init
//------------------
jerror_t DPSHit_factory::init(void)
{
	/// set the base conversion scales
	a_scale    = 0.0001; 
	t_scale    = 0.0625;   // 62.5 ps/count
	t_base     = 0.;    // ns

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DPSHit_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	/// Read in calibration constants
	jout << "In DPSHit_factory, loading constants..." << endl;

	// extract the PS Geometry
	vector<const DPSGeometry*> psGeomVect;
	eventLoop->Get( psGeomVect );
	if (psGeomVect.size() < 1)
		return OBJECT_NOT_AVAILABLE;
	const DPSGeometry& psGeom = *(psGeomVect[0]);

	// load scale factors
	map<string,double> scale_factors;
	if (eventLoop->GetCalib("/PHOTON_BEAM/pair_spectrometer/fine/digi_scales", scale_factors))
		jout << "Error loading /PHOTON_BEAM/pair_spectrometer/fine/digi_scales !" << endl;
	if (scale_factors.find("PS_ADC_ASCALE") != scale_factors.end())
		a_scale = scale_factors["PS_ADC_ASCALE"];
	else
		jerr << "Unable to get PS_ADC_ASCALE from /PHOTON_BEAM/pair_spectrometer/fine/digi_scales !" 
		     << endl;
	if (scale_factors.find("PS_ADC_TSCALE") != scale_factors.end())
		t_scale = scale_factors["PS_ADC_TSCALE"];
	else
		jerr << "Unable to get PS_ADC_TSCALE from /PHOTON_BEAM/pair_spectrometer/fine/digi_scales !" 
		     << endl;

	// load base time offset
	map<string,double> base_time_offset;
	if (eventLoop->GetCalib("/PHOTON_BEAM/pair_spectrometer/fine/base_time_offset",base_time_offset))
		jout << "Error loading /PHOTON_BEAM/pair_spectrometer/fine/base_time_offset !" << endl;
	if (base_time_offset.find("PS_BASE_TIME_OFFSET") != base_time_offset.end())
		t_base = base_time_offset["PS_BASE_TIME_OFFSET"];
	else
		jerr << "Unable to get PS_BASE_TIME_OFFSET from /PHOTON_BEAM/pair_spectrometer/fine/base_time_offset !" << endl;


        FillCalibTable(adc_pedestals, "/PHOTON_BEAM/pair_spectrometer/fine/adc_pedestals", psGeom);
        FillCalibTable(adc_gains, "/PHOTON_BEAM/pair_spectrometer/fine/adc_gain_factors", psGeom);
        FillCalibTable(adc_time_offsets, "/PHOTON_BEAM/pair_spectrometer/fine/adc_timing_offsets", psGeom);

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DPSHit_factory::evnt(JEventLoop *loop, int eventnumber)
{
	/// Generate DPSHit object for each DPSDigiHit object.
	/// This is where the first set of calibration constants
	/// is applied to convert from digitzed units into natural
	/// units.
	///
	/// Note that this code does NOT get called for simulated
	/// data in HDDM format. The HDDM event source will copy
	/// the precalibrated values directly into the _data vector.

	// extract the PS Geometry
	vector<const DPSGeometry*> psGeomVect;
	eventLoop->Get( psGeomVect );
	if (psGeomVect.size() < 1)
		return OBJECT_NOT_AVAILABLE;
	const DPSGeometry& psGeom = *(psGeomVect[0]);

	// First, make hits out of all fADC250 hits
	vector<const DPSDigiHit*> digihits;
	loop->Get(digihits);
	char str[256];

	for (unsigned int i=0; i < digihits.size(); i++){
		const DPSDigiHit *digihit = digihits[i];

		DPSHit *hit = new DPSHit;

		// Make sure channel id is in valid range
		if( (digihit->arm < 0) && (digihit->arm >= psGeom.NUM_ARMS) ) {
			sprintf(str, "DPSDigiHit arm out of range! arm=%d (should be 0-%d)", 
				digihit->arm, psGeom.NUM_ARMS);
			throw JException(str);
		}
		if( (digihit->arm <= 0) && (digihit->column > psGeom.NUM_FINE_COLUMNS) ) {
			sprintf(str, "DPSDigiHit column out of range! column=%d (should be 0-%d)", 
				digihit->column, psGeom.NUM_FINE_COLUMNS);
			throw JException(str);
		}
		
		// The translation table has PS channels labaled as paddles 1-16
		// The PSHit class labels hits as
		//   arm:     North/South (0/1)
		//   column:  1-185
		hit->arm     = digihit->arm;
		hit->column  = digihit->column;

                // Get pedestal.  Prefer associated event pedestal if it exists.
                // Otherwise, use the average pedestal from CCDB
                double pedestal = GetConstant(adc_pedestals,digihit,psGeom);
                const Df250PulseIntegral* PIobj = NULL;
                const Df250Config *configObj = NULL;
                digihit->GetSingle(PIobj);
                PIobj->GetSingle(configObj);
                if ((PIobj != NULL) && (configObj != NULL)) {
                        // the measured pedestal must be scaled by the ratio of the number
                        // of samples used to calculate the pedestal and the actual pulse
                        pedestal = static_cast<double>(configObj->NSA_NSB) * PIobj->pedestal;
                }

                // Apply calibration constants here
                double A = (double)digihit->pulse_integral;
                double T = (double)digihit->pulse_time;

		hit->dE = a_scale * (A - pedestal);
                hit->t = t_scale * (T - GetConstant(adc_time_offsets, digihit, psGeom)) + t_base;
                hit->sigma_t = 4.0;    // ns (what is the fADC time resolution?)
                hit->has_fADC = true;
                hit->has_TDC  = false; // will get set to true below if appropriate

		hit->AddAssociatedObject(digihit);
                
                _data.push_back(hit);
        }


	return NOERROR;
}



//------------------
// erun
//------------------
jerror_t DPSHit_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DPSHit_factory::fini(void)
{
	return NOERROR;
}

//------------------
// FillCalibTable
//------------------
void DPSHit_factory::FillCalibTable(ps_digi_constants_t &table, string table_name, 
                                     const DPSGeometry &psGeom)
{
	char str[256];

        // load constant table
        if(eventLoop->GetCalib(table_name, table))
		jout << "Error loading " + table_name + " !" << endl;

	
	// check that the size of the tables loaded are correct
	if( (int)table.size() != psGeom.NUM_FINE_COLUMNS ) {
		sprintf(str, "PS table loaded with wrong size! number of columns=%d (should be %d)", 
			(int)table.size(), psGeom.NUM_FINE_COLUMNS );
		cerr << str << endl;
		throw JException(str);
	}
	
	for( int column=0; column < psGeom.NUM_FINE_COLUMNS; column++) {
		if( (int)table[column].size() != psGeom.NUM_ARMS ) {
			sprintf(str, "PS table loaded with wrong size! column=%d number of arms=%d (should be %d)", 
				column, (int)table[column].size(), psGeom.NUM_ARMS );
			cerr << str << endl;
			throw JException(str);
		}
	}
}

//------------------------------------
// GetConstant
//   Allow a few different interfaces
//
//   PS Geometry as defined in the Translation Table:
//       arm:     North/South (0/1)
//       column:  1-145
//   Note the different counting schemes used
//------------------------------------
const double DPSHit_factory::GetConstant( const ps_digi_constants_t &the_table, 
					   const DPSGeometry::Arm in_arm, const int in_column,
					   const DPSGeometry &psGeom ) const
{
        char str[256];
        
        if( (in_arm != DPSGeometry::kNorth) && (in_arm != DPSGeometry::kSouth)) {
                sprintf(str, "Bad arm requested in DPSHit_factory::GetConstant()! requested=%d , should be 0-%d", 
			static_cast<int>(in_arm), static_cast<int>(DPSGeometry::kSouth));
                cerr << str << endl;
                throw JException(str);
        }
        if( (in_column <= 0) || (in_column > psGeom.NUM_COARSE_COLUMNS)) {
                sprintf(str, "Bad column # requested in DPSHit_factory::GetConstant()! requested=%d , should be 1-%d", in_column, psGeom.NUM_COARSE_COLUMNS);
                cerr << str << endl;
                throw JException(str);
        }

	// the tables are indexed by column, with the different values for the two arms
	// stored in the two fields of the pair
        if(in_arm == DPSGeometry::kNorth) {
		return the_table[in_column-1][in_arm];
        } else {
		return the_table[in_column-1][in_arm];
        }
}

const double DPSHit_factory::GetConstant( const ps_digi_constants_t &the_table, 
					   const DPSHit *in_hit, const DPSGeometry &psGeom ) const
{
        char str[256];
        
        if( (in_hit->arm != DPSGeometry::kNorth) && (in_hit->arm != DPSGeometry::kSouth)) {
                sprintf(str, "Bad arm requested in DPSHit_factory::GetConstant()! requested=%d , should be 0-%d", 
			static_cast<int>(in_hit->arm), static_cast<int>(DPSGeometry::kSouth));
                cerr << str << endl;
                throw JException(str);
        }
        if( (in_hit->column <= 0) || (in_hit->column > psGeom.NUM_COARSE_COLUMNS)) {
                sprintf(str, "Bad column # requested in DPSHit_factory::GetConstant()! requested=%d , should be 1-%d", in_hit->column, psGeom.NUM_COARSE_COLUMNS);
                cerr << str << endl;
                throw JException(str);
        }

	// the tables are indexed by column, with the different values for the two arms
	// stored in the two fields of the pair
        if(in_hit->arm == DPSGeometry::kNorth) {
		return the_table[in_hit->column-1][in_hit->arm];
        } else {
		return the_table[in_hit->column-1][in_hit->arm];
        }
}

const double DPSHit_factory::GetConstant( const ps_digi_constants_t &the_table, 
					   const DPSDigiHit *in_digihit, const DPSGeometry &psGeom) const
{
        char str[256];

        if( (in_digihit->arm != DPSGeometry::kNorth) && (in_digihit->arm != DPSGeometry::kSouth)) {
                sprintf(str, "Bad arm requested in DPSHit_factory::GetConstant()! requested=%d , should be 0-%d", 
			static_cast<int>(in_digihit->arm), static_cast<int>(DPSGeometry::kSouth));
                cerr << str << endl;
                throw JException(str);
        }
        if( (in_digihit->column <= 0) || (in_digihit->column > psGeom.NUM_COARSE_COLUMNS)) {
                sprintf(str, "Bad column # requested in DPSHit_factory::GetConstant()! requested=%d , should be 1-%d", in_digihit->column, psGeom.NUM_COARSE_COLUMNS);
                cerr << str << endl;
                throw JException(str);
        }

	// the tables are indexed by column, with the different values for the two arms
	// stored in the two fields of the pair
        if(in_digihit->arm == DPSGeometry::kNorth) {
		return the_table[in_digihit->column-1][in_digihit->arm];
        } else {
		return the_table[in_digihit->column-1][in_digihit->arm];
        }
}

