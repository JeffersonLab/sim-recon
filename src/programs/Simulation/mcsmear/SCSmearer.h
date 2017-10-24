// Smearing class for start counter (SC)

#ifndef _SCSMEARER_H_
#define _SCSMEARER_H_

#include "Smearer.h"


class sc_config_t 
{
  public:
	sc_config_t(JEventLoop *loop);

    double GetPaddleTimeResolution(int sector, double sc_local_z)  { 
        double time_resolution = 0.;

        if(sc_local_z < SC_BOUNDARY1[sector]) {
            time_resolution = SC_SECTION1_P0[sector] + SC_SECTION1_P1[sector]*sc_local_z;
        } else if(sc_local_z < SC_BOUNDARY2[sector]) {
            time_resolution = SC_SECTION2_P0[sector] + SC_SECTION2_P1[sector]*sc_local_z;
        } else {
            time_resolution = SC_SECTION3_P0[sector] + SC_SECTION3_P1[sector]*sc_local_z;
        }
        
        // max sure that we aren't getting some ridiculously large resolution
        if(time_resolution > SC_MAX_RESOLUTION[sector])
            time_resolution = SC_MAX_RESOLUTION[sector];
    
        // If these resolutions come from data, apply correction factors to remove any other contributions
        time_resolution = (time_resolution - SC_MC_CORRECTION_P0) / SC_MC_CORRECTION_P1;

        // convert ps to ns
        time_resolution /= 1000.;
        cout << " time resolution = " << time_resolution << endl;
        return time_resolution;
	}

	double GetEfficiencyCorrectionFactor(int sector) {
		return paddle_efficiencies.at(sector-1);
	}
	
	double START_SIGMA;
	double START_PHOTONS_PERMEV;
	double START_PADDLE_THRESHOLD;

	vector<double> paddle_efficiencies;

    // Start counter geometry parameters
    vector<double> SC_START_Z;

    // Start counter resolution parameters
    vector<double> SC_MAX_RESOLUTION;
    vector<double> SC_BOUNDARY1, SC_BOUNDARY2;
    vector<double> SC_SECTION1_P0, SC_SECTION1_P1;
    vector<double> SC_SECTION2_P0, SC_SECTION2_P1;
    vector<double> SC_SECTION3_P0, SC_SECTION3_P1;

    double SC_MC_CORRECTION_P0, SC_MC_CORRECTION_P1;

};


class SCSmearer : public Smearer
{
  public:
	SCSmearer(JEventLoop *loop, mcsmear_config_t *in_config) : Smearer(loop, in_config) {
		sc_config = new sc_config_t(loop);
	}
	~SCSmearer() {
		delete sc_config;
	}
	
	void SmearEvent(hddm_s::HDDM *record);

  protected:
    hddm_s::StcTruthPointList::iterator FindMatchingTruthPoint(hddm_s::StcTruthHitList::iterator hiter, hddm_s::StcTruthPointList &truthPoints);
	
  private:
  	sc_config_t  *sc_config;
};



#endif // _SCSMEARER_H_
