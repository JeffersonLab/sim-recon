/*
 * hitBCal - registers hits for barrel calorimeter
 *
 *	This is a part of the hits package for the
 *	HDGeant simulation program for Hall D.
 *
 *	version 1.0 	-Richard Jones July 16, 2001
 *
 *  changes:
 *          Wed Feb 22 08:06:31 EST 2012 D. Lawrence
 *          - Convert to .cc suffix to compile with C++ compiler
 *          - Modify to store full enrgy-weighted time spectra
 *            in HDDM file
 *
 *          Wed Jun 20 13:19:56 EDT 2007 B. Zihlmann 
 *          add ipart to the function hitBarrelEMcal
 *
 *          Fri Jul 13 08:54:48 EDT 2007 M. Shepherd
 *          remove attenuation, condense up and dowstream hits
 *          pass up true z position and hit time
 */

//-----------------
// A brief overview on how this works:
//
// Most of the hitXXX.c(c) files work in a similar way: They
// have routines like hitBarrelEMcal that get called from
// savehits.F which itself gets called from gustep.F. Thus,
// hitBarrelEMcal gets called for every step of every particle
// that is in a volume of the BCAL. The routine accumulates 
// the data for the event, creating multiple hits etc. as needed
// and stores this in a binTree_t structure. The binTree_t
// is essentially a C-structure that implements a doubly-linked
// list with each node containing a pointer to a structure that
// actually holds the data. These structures tend to be HDDM
// defined data structures. The bintree mechanism uses a "mark"
// to identify a node and is detector system specific. Bits in
// the "mark" are used to identify the detector element
// (e.g. module, layer, sector) as well as the type of data
// (truth or hit). The bintree stores the pointers as void*
// so that it can be used for any type of data structure. This
// means though, that one must know the type and cast it appropriately
// when extracting a pointer from the bintree.
//
// After all of the particles in the event have been tracked,
// GEANT calls guout which in turn calls hddmOutput. This calls 
// the individual "pick" routines like pickBarrelEMcal. The
// pickBarrelEMcal routine goes through everything in the bintree
// (which happen to be stored in seperate s_BarrelEMcal_t structures)
// and combines them into a single s_BarrelEMcal_t structure
// that is in the hddm_s structure that will be written to
// the output file.
//
// For the case of the BCAL timing spectra, a similar method is
// used except global variable STL containers are used instead of
// the bintree mechanism.
//-----------------


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <map>
#include <set>
#include <sstream>
using namespace std;

extern "C" {
#include <hddm_s.h>
#include <geant3.h>
#include <bintree.h>

#include "calibDB.h"
}

#include <DHistogram.h>

// Use the following to switch on and off the single hit
// and full spectra methods of writing out the BCAL info.
#define WRITE_OUT_BCAL_CELL_HITS 1
#define WRITE_OUT_BCAL_TIME_SPECTRA 1


extern float BGGATE1;
extern float BGGATE2;


static float THRESH_MEV	     = 1.;
static float TWO_HIT_RESOL   = 50.;
static int   MAX_HITS 	     = 100;

binTree_t* barrelEMcalTree = 0;
static int cellCount = 0;
static int showerCount = 0;
static int initialized = 0;

// The following will be changed in initializeBarrelEMcal
static double ATTEN_LENGTH = 300.0;
static double C_EFFECTIVE = 16.75;
static double SiPM_tbin_width = 0.100;
static double atten_full_length = exp(-390.0/ATTEN_LENGTH); // attenuation factor for length of full BCAL module
static double THRESH_ATTENUATED_GEV = (THRESH_MEV/1000.0)*atten_full_length;
static int NHISTO_BINS = 0; 

extern s_HDDM_t* thisInputEvent;

// Prevent name mangling so these routines keep their
// C-style names in the object
extern "C"{
	void recordbcalentry_(int *mech, int *itra, int*istak, int *ipart, float *vect, float *getot);
	void hitBarrelEMcal (float xin[4], float xout[4],
                     	float pin[5], float pout[5], float dEsum,
                     	int track, int stack, int history, int ipart);
	void hitbarrelemcal_(float* xin, float* xout,
                     	float* pin, float* pout, float* dEsum,
                     	int* track, int* stack, int* history, int* ipart);
	s_BarrelEMcal_t* pickBarrelEMcal();
}

//..........................
// bcal_index is a utility class that encapsulates the
// module, layer, sector, and end in a single object that
// can be used as a key to index an STL map. 
//..........................
class bcal_index{
	public:
		enum EndType{
			kUp,
			kDown
		};
		
		bcal_index(unsigned int module, unsigned int layer, unsigned int sector, unsigned int incident_id, EndType end):
			module(module),layer(layer),sector(sector),incident_id(incident_id),end(end){}
	
		unsigned int module;
		unsigned int layer;
		unsigned int sector;
		unsigned int incident_id;
		EndType end;
		
		bool operator<(const bcal_index &idx) const{
			if(module<idx.module)return true;
			if(module>idx.module)return false;
			if(layer<idx.layer)return true;
			if(layer>idx.layer)return false;
			if(sector<idx.sector)return true;
			if(sector>idx.sector)return false;
			if(incident_id<idx.incident_id)return true;
			if(incident_id>idx.incident_id)return false;
			if((end==kUp) && (idx.end==kDown))return true;
			return false;
		}

		// For debugging
		string ToString(void) const {
			stringstream ss;
			ss << "module:"<<module<<" layer:"<<layer<<" sector:"<<sector<<" incident_id:"<<incident_id;
			return ss.str();
		}
};

// This holds the energy-weighted time spectra for
// each SiPM in the event. At the begining of the 
// event it is empty and entries are added as cells
// are hit.
map<bcal_index, DHistogram*> SiPMspectra;

// Keep a pool of DHistogram objects that recycled
// event to minimize memory usage.
vector<DHistogram*> histo_pool;

// This is used to hold the integrated energy in GeV for
// each cell. This is only for the SiPMspectra hits. This 
// map will be filled with entries for each cell at initialization
// time, but only those elements corresponding to entries in
// histos_with_data will be reset at the end of the event.
// n.b. an entry will be made for each end, even though they
// will be identical for each cell. It turns out to make it 
// simpler in the pickBarrelEMcal (and elsewhere) to just
// keep the duplicate information.
map<bcal_index, double> E_CellTruth;

// The IncidentParticle_t class is used to hold information on
// a single particle entering the BCAL. This is used when timing
// spectra are written out so that sampling fluctuations (and
// any other smearing that depends on incident particle parameters)
// may be done. A list of up to MAX_INCIDENT_PARTICLES is kept
// in the BCAL_INCIDENT_PARTICLES vector. All particles in the
// vector are written to the output file with an id that is also
// written with each SiPM spectrum.
class IncidentParticle_t{
	public:
		IncidentParticle_t(const float *v, const float getot, int ptype, int track):x(v[0]),y(v[1]),z(v[2]),px(v[3]*v[6]),py(v[4]*v[6]),pz(v[5]*v[6]),E(getot),ptype(ptype),track(track){}
		float x,y,z;
		float px, py, pz;
		float E;
		int ptype, track;
		float dPhi(const IncidentParticle_t &pos){float a=(pos.x*x + pos.y*y)/sqrt((x*x + y*y)*(pos.x*pos.x+pos.y*pos.y)); return a<1.0 ? fabs(acos(a)):0.0;}
		float dZ(const IncidentParticle_t &pos){return fabs(pos.z-z);}
};
#define MAX_INCIDENT_PARTICLES 100
vector<IncidentParticle_t> BCAL_INCIDENT_PARTICLES;
int BCAL_INCIDENT_PARTICLE_COUNT = 0;
vector<int> INCIDENT_ID; // holds map of tracks to incident particle id
bool SHOWED_INCIDENT_PARTICLE_LONG_WARNING = false;
bool SHOWED_INCIDENT_PARTICLE_SHORT_WARNING = false;

//----------------------
// GetSiPMHistogram
//----------------------
DHistogram* GetSiPMHistogram(const bcal_index &idx)
{
	/// Look to see if an entry in the SiPMspectra map
	/// already exists for the specified index. If so,
	/// return it. If not, a histogram is assigned to
	/// it from the pool. If the pool is empty, then
	/// a new histogram is created.
	map<bcal_index, DHistogram*>::iterator iter = SiPMspectra.find(idx);

	// If we found it, return immediately
	if(iter != SiPMspectra.end()) return iter->second;

	// If we didn't find it, try getting one from the pool
	DHistogram *h;
	if(histo_pool.size()>0){
		// Pool has one, take it out to use it
		h = histo_pool.back();
		histo_pool.pop_back();
	}else{
		// Pool does not have one. Create a new one. This
		// will be pushed onto the pool at the end of the
		// event.
		h = new DHistogram(NHISTO_BINS, BGGATE1, BGGATE2);
	}
	
	// Assign this histo to the SiPMspectra map
	SiPMspectra[idx] = h;

	// Return a pointer to the histogram
	return h;
}

//----------------------
// initializeBarrelEMcal
//----------------------
void initializeBarrelEMcal(void)
{
   mystr_t strings[50];
   float values[50];
   int nvalues = 50;

	// Get parameters from BCAL/bcal_parms
   int status = GetConstants("BCAL/bcal_parms", &nvalues, values, strings);    
   if (!status) {
     int ncounter = 0;
     for (int i=0;i<(int)nvalues;i++){
       if (!strcmp(strings[i].str,"BCAL_THRESH_MEV")) {
         THRESH_MEV  = values[i];
         ncounter++;
       }
       if (!strcmp(strings[i].str,"BCAL_TWO_HIT_RESOL")) {
         TWO_HIT_RESOL  = values[i];
         ncounter++;
       }
       if (!strcmp(strings[i].str,"BCAL_MAX_HITS")) {
         MAX_HITS  = (int)values[i];
         ncounter++;
       }
     }
     if (ncounter==3){
       printf("BCAL/bcal_parms: ALL parameters loaded from Data Base\n");
     } else if (ncounter<3){
       printf("BCAL/bcal_parms: NOT ALL necessary parameters found in Data Base %d out of 3\n",ncounter);
     } else {
       printf("BCAL/bcal_parms: SOME parameters found more than once in Data Base\n");
     }
   }

	// Get parameters from BCAL/mc_parms
	nvalues = 50;
   status = GetConstants("BCAL/mc_parms", &nvalues, values, strings);
   if (!status) {
     int ncounter = 0;
     for (int i=0;i<(int)nvalues;i++){
       if (!strcmp(strings[i].str,"ATTEN_LENGTH")) {
         ATTEN_LENGTH  = values[i];
         ncounter++;
       }
       if (!strcmp(strings[i].str,"C_EFFECTIVE")) {
         C_EFFECTIVE  = values[i];
         ncounter++;
       }
       if (!strcmp(strings[i].str,"SiPM_tbin_width")) {
         SiPM_tbin_width  = values[i];
         ncounter++;
       }
     }
     if (ncounter==3){
       printf("BCAL/mc_parms: ALL parameters loaded from Data Base\n");
     } else if (ncounter<3){
       printf("BCAL/mc_parms: NOT ALL necessary parameters found in Data Base %d out of 3\n",ncounter);
     } else {
       printf("BCAL/mc_parms: SOME parameters found more than once in Data Base\n");
     }
   }

	// Factors to apply effective threshold on attenuated signal.
	atten_full_length = exp(-390.0/ATTEN_LENGTH);
	THRESH_ATTENUATED_GEV = (THRESH_MEV/1000.0)*atten_full_length;

#if WRITE_OUT_BCAL_TIME_SPECTRA
	// Calculate timing spectrum histogram parameters and create
	// E_CellTruth entries for each SiPM.
	//
	// We want the histos to cover the range specified by the BGGATE
	// card in the control.in file (if it is set). However, we want
	// to honor the value of SiPM_tbin_width from the calibration DB.
	// Thus, we calculate the number of bins, but then adjust the
	// gate width to allow an integral number of bins of the specified
	// width.
	// NOTE: We change the values of BGGATE1 and BGGATE2 here, but
	// these are copies of the values used by the FORTRAN code. 
	// Changing them here doesn't change the values used by FORTRAN.
	if(BGGATE2<=BGGATE1){
	  BGGATE1 = -200.0;
	  BGGATE2 = +200.0;
	}
	NHISTO_BINS = (int)floor(0.5 + (BGGATE2 - BGGATE1)/SiPM_tbin_width);
	BGGATE2 = BGGATE1 + SiPM_tbin_width*(double)NHISTO_BINS;
	
	for(unsigned int imodule=1; imodule<=48; imodule++){
	  for(unsigned int ilayer=1; ilayer<=10; ilayer++){
	 	  for(unsigned int isector=1; isector<=4; isector++){
	 	     for(unsigned int incident_id=1; incident_id<=MAX_INCIDENT_PARTICLES; incident_id++){
				  bcal_index idxUp(imodule, ilayer, isector, incident_id, bcal_index::kUp);
				  bcal_index idxDn(imodule, ilayer, isector, incident_id, bcal_index::kDown);

				  E_CellTruth[idxUp] = 0.0;
				  E_CellTruth[idxDn] = 0.0;
		     }
		  }
	  }
	}
#endif  //WRITE_OUT_BCAL_TIME_SPECTRA

   initialized = 1;
}

//----------------------
// recordbcalentry
//----------------------
void recordbcalentry_(int *mech, int *itra, int*istak, int *ipart, float *vect, float *getot)
{
#if WRITE_OUT_BCAL_TIME_SPECTRA
	// This gets called from gustep whenever a particle is
	// "entering" one of the inner 6 BCAL layers (volumes whose name
	// starts with "BM0".) Its purpose is to record the parameters
	// for all particles entering the BCAL so they can be used in
	// mcsmear to apply the appropriate sampling fluctuations.
	//
	// A complication occurs in that particles created in the BCAL
	// will also be flagged as "entering" it. So, we check to see if
	// the particle of the current step is close to one already
	// recorded and if so, assume it is part of the same shower so
	// don't record it again. Low energy particles are also igonored.

	IncidentParticle_t mypart(vect, *getot, *ipart, *itra);
	
	bool add_to_list = true;
	float dPhi, dZ;
	for(unsigned int i=0; i<BCAL_INCIDENT_PARTICLES.size(); i++){
		
		// Only keep photons and betas 
		//if(*ipart>3)add_to_list = false;
		
		dPhi = 1000.0*mypart.dPhi(BCAL_INCIDENT_PARTICLES[i]);
		dZ = mypart.dZ(BCAL_INCIDENT_PARTICLES[i]);
		// if this is within 200 mrad and 30cm of a previously recorded
		// particle entering BCAL, assume it is part of the same shower
		// Also, ignore particles with less than 10MeV total energy.
		if(dPhi<200.0 && dZ<30.0){
			add_to_list = false;
			
			// If this particle has larger total energy than the one
			// already recorded, then replace it with this one. This 
			// would be for the case when a shower sprays from something
			// like the FDC frame so many particles enter the same area
			// but are too close together to be considered separate showers.
			if(mypart.E>BCAL_INCIDENT_PARTICLES[i].E){
				BCAL_INCIDENT_PARTICLES[i] = mypart;
			}
		}
		if(*getot<0.100)add_to_list = false;
		if(!add_to_list)break;
	}
	if(add_to_list){
		BCAL_INCIDENT_PARTICLE_COUNT++;
//_DBG_<<"*itra = "<<*itra<<"  dPhi="<<dPhi<<" dZ="<<dZ<<endl;
		if(BCAL_INCIDENT_PARTICLES.size()>=MAX_INCIDENT_PARTICLES){
			if(!SHOWED_INCIDENT_PARTICLE_LONG_WARNING){
				cerr<<endl;
				cerr<<"WARNING: The BCAL records information about certain"<<endl;
				cerr<<"particles entering it so that information can be used"<<endl;
				cerr<<"later in mcsmear to properly smear the signals. For"<<endl;
				cerr<<"this event, more than the maximum number of incident"<<endl;
				cerr<<"particles has occurred ("<<MAX_INCIDENT_PARTICLES<<") so the list is"<<endl;
				cerr<<"being truncated to the first "<<MAX_INCIDENT_PARTICLES<<"."<<endl;
				cerr<<"All of the signal in the BCAL is still being recorded,"<<endl;
				cerr<<"but the smearing may be off by a few percent for this"<<endl;
				cerr<<"event. It is probably nothing to worry about. This long"<<endl;
				cerr<<"message will only appear once and the following line just"<<endl;
				cerr<<"once per event whenever this occurs."<<endl;
				cerr<<endl;
				SHOWED_INCIDENT_PARTICLE_LONG_WARNING = true;
			}
			if(!SHOWED_INCIDENT_PARTICLE_SHORT_WARNING){
				cerr<<__FILE__<<":"<<__LINE__<<" too many particles entering BCAL! Some information will be lost."<<endl;
				SHOWED_INCIDENT_PARTICLE_SHORT_WARNING = true;
			}
		}else{
			BCAL_INCIDENT_PARTICLES.push_back(mypart);
		}
	}
#endif
}

//----------------------
// find_incident_id
//
// Find the entry in BCAL_INCIDENT_PARTICLES that is physically
// closest to the given location
//----------------------
unsigned int find_incident_id(float *x)
{
	// This should probably use the distance to the line
	// projected by the incident particle, but that would
	// be computationally expensive for something that is
	// called for every step in the BCAL shower development.
	//
	// Even better, one could search through the particle
	// stack tracing the ancestory until a track number
	// associated with an incident particle is found. That
	// too would be expensive. 

	unsigned int closest_id=0;
	float closest_dist2 = 1.0E6;
	for(unsigned int i=0; i<BCAL_INCIDENT_PARTICLES.size(); i++){
		float dx = x[0] - BCAL_INCIDENT_PARTICLES[i].x;
		float dy = x[1] - BCAL_INCIDENT_PARTICLES[i].y;
		float dz = x[2] - BCAL_INCIDENT_PARTICLES[i].z;
		float dist2 = dx*dx + dy*dy + dz*dz;
		if(dist2 < closest_dist2){
			closest_dist2 = dist2;
			closest_id = i+1;
		}
	}

	return closest_id;
}

//----------------------
// hitBarrelEMcal
//
// register hits during tracking (from gustep)
//----------------------
void hitBarrelEMcal (float xin[4], float xout[4],
                     float pin[5], float pout[5], float dEsum,
                     int track, int stack, int history, int ipart)
{
	float x[3], t;
	float xlocal[3];
	float xbcal[3];
	float xHat[] = {1,0,0};
  
	if(!initialized)initializeBarrelEMcal();
  
	x[0] = (xin[0] + xout[0])/2;
	x[1] = (xin[1] + xout[1])/2;
	x[2] = (xin[2] + xout[2])/2;
	t    = (xin[3] + xout[3])/2 * 1e9;
	transformCoord(x,"global",xlocal,"BCAL");
	transformCoord(xHat,"local",xbcal,"BCAL");

//cout<<"track:"<<track<<" stack:"<<stack<<" history:"<<history<<" ipart:"<<ipart<<endl;

	/* Under certain conditions the time in xout[3] will
	  be invalid (unusually large). Check for this and
	  use only the in time in these cases
	*/
	if(xout[3] > 1.0){
	 t = xin[3] * 1e9;
	}
  
	/* post the hit to the truth tree */

	if ((history == 0) && (pin[3] > THRESH_MEV/1e3)) {
		s_BcalTruthShowers_t* showers;
		float r = sqrt(xin[0]*xin[0]+xin[1]*xin[1]);
		float phi = atan2(xin[1],xin[0]);
		int mark = (1<<30) + showerCount;
		void** twig = getTwig(&barrelEMcalTree, mark);
		if (*twig == 0) 
		{
			s_BarrelEMcal_t* bcal = make_s_BarrelEMcal();
			*twig = bcal;
			bcal->bcalTruthShowers = showers = make_s_BcalTruthShowers(1);
			int a = thisInputEvent->physicsEvents->in[0].reactions->in[0].vertices->in[0].products->mult;
			showers->in[0].primary = (track <= a);
			showers->in[0].track = track;
			showers->in[0].z = xin[2];
			showers->in[0].r = r;
			showers->in[0].phi = phi;
			showers->in[0].t = xin[3]*1e9;
			showers->in[0].px = pin[0]*pin[4];
			showers->in[0].py = pin[1]*pin[4];
			showers->in[0].pz = pin[2]*pin[4];
			showers->in[0].E = pin[3];
			showers->in[0].ptype = ipart;
			showers->mult = 1;
			showerCount++;
		}
	}

#if WRITE_OUT_BCAL_CELL_HITS
	// ......................................................
	// This section, contains the original
	// mechanism that recorded an energy weighted time average for each
	// cell. It has since been replaced with the one right after it
	// that records full timing spectra for each hit cell. One can
	// however, re-eanble this with the WRITE_OUT_BCAL_CELL_HITS
	// flag set at the top of the file.
  
  /* post the hit to the hits tree, mark sector as hit */

   if (dEsum > 0)
   {
      int nshot;
      s_BcalHits_t* hits;
      int sector = getsector_wrapper_();
      int layer  = getlayer_wrapper_();
      int module = getmodule_wrapper_();
      float zLocal = xlocal[2];
      int mark = (module<<16)+ (layer<<9) + sector;
      
      void** twig = getTwig(&barrelEMcalTree, mark);
      if (*twig == 0)
      {
         s_BarrelEMcal_t* bcal = make_s_BarrelEMcal();
	 		*twig = bcal;
         s_BcalCells_t* cells = make_s_BcalCells(1);
         cells->mult = 1;
         cells->in[0].module = module;
         cells->in[0].layer = layer;
         cells->in[0].sector = sector;
         cells->in[0].bcalHits = hits
                               = make_s_BcalHits(MAX_HITS);
         bcal->bcalCells = cells;
         cellCount++;
      }
      else
      {
         s_BarrelEMcal_t* bcal = (s_BarrelEMcal_t*)*twig;
         hits = bcal->bcalCells->in[0].bcalHits;
      }

      for (nshot = 0; nshot < (int)hits->mult; nshot++)
      {
         if (fabs(hits->in[nshot].t - t) < TWO_HIT_RESOL)
         {
            break;
         }
      }
      if (nshot < (int)hits->mult)		/* merge with former hit */
      {
         hits->in[nshot].t =
                  (hits->in[nshot].t * hits->in[nshot].E + t * dEsum)
                / (hits->in[nshot].E + dEsum);
         hits->in[nshot].zLocal =
                  (hits->in[nshot].zLocal * hits->in[nshot].E + zLocal * dEsum)
                / (hits->in[nshot].E + dEsum);
         hits->in[nshot].E += dEsum;
      }
      else if (nshot < MAX_HITS)		/* create new hit */
      {
         hits->in[nshot].t = t;
         hits->in[nshot].E = dEsum;
         hits->in[nshot].zLocal = zLocal;
         hits->mult++;
      }
      else
      {
         fprintf(stderr,"HDGeant error in hitBarrelEMcal: ");
         fprintf(stderr,"max hit count %d exceeded, truncating!\n",MAX_HITS);
      }
   }
	// ......................................................
#endif  // WRITE_OUT_BCAL_CELL_HITS

	 
#if WRITE_OUT_BCAL_TIME_SPECTRA
	// Fill the data structures that hold the energy weighted timing
	// spectra for the cells that are hit. 
	 
	// Note that this is done differently from other hits in HDGeant in
	// that the structures filled here are not HDDM defined structures.
	// The data is copied (sparsely) into the HDDM structures in
	// pickBarrelEMcal below (as they are for other the cell hit and
	// truth structures).
	if (dEsum > 0)
	{
		// We need the incident particle id associated with the current
		// track. This is kept in an STL map so the association is made
		// only once. The way this works is a little subtle so is worth 
		// explaining:
		// When access the map via "[]" operator, an entry in the map
		// will be created if one does not already exist. For int types,
		// the new entry will have a value of "0". By using particle id
		// values starting from "1", we can test whether the association
		// has been made.
		// The map index is the track number and the value the
		// incident particle id. We make the incident particle id numbers
		// start from 1 so if we access the map

		// Guarantee a value exists at INCIDENT_ID[track]
		if(track>=(int)INCIDENT_ID.size())INCIDENT_ID.resize(track+1, 0);
		int incident_id = INCIDENT_ID[track];
		if(incident_id==0)incident_id = INCIDENT_ID[track] = find_incident_id(x);
	
		// Get map index based on layer and sector
      unsigned int sector = getsector_wrapper_();
      unsigned int layer  = getlayer_wrapper_();
      unsigned int module = getmodule_wrapper_();
		bcal_index idxUp(module, layer, sector, incident_id, bcal_index::kUp);
		bcal_index idxDn(module, layer, sector, incident_id, bcal_index::kDown);
		
		DHistogram *hup = GetSiPMHistogram(idxUp);
		DHistogram *hdn = GetSiPMHistogram(idxDn);

		double Z = xlocal[2];
		double dist_up = 390.0/2.0 + Z;
		double dist_dn = 390.0/2.0 - Z;

		double Eup = dEsum*exp(-dist_up/ATTEN_LENGTH);
		double Edn = dEsum*exp(-dist_dn/ATTEN_LENGTH);
		double t_up = t + dist_up/C_EFFECTIVE;
		double t_dn = t + dist_dn/C_EFFECTIVE;

		hup->Fill(t_up, Eup);
		hdn->Fill(t_dn, Edn);

		E_CellTruth[idxUp] += dEsum;
		E_CellTruth[idxDn] += dEsum;			
	}
#endif //WRITE_OUT_BCAL_TIME_SPECTRA
}

//----------------------
// hitbarrelemcal_
//
// entry point from fortran
//----------------------
void hitbarrelemcal_(float* xin, float* xout,
                     float* pin, float* pout, float* dEsum,
                     int* track, int* stack, int* history, int* ipart)
{
   hitBarrelEMcal(xin,xout,pin,pout,*dEsum,*track,*stack,*history, *ipart);
}


//----------------------
// pickBarrelEMcal
//
// pick and package the hits for shipping
//----------------------
s_BarrelEMcal_t* pickBarrelEMcal ()
{
   s_BarrelEMcal_t* box;  // pointer to structure we're copying into
   s_BarrelEMcal_t* item; // temporary pointer to structure we're copying from
	unsigned int spectraCount = SiPMspectra.size(); // Number non-zero spectra
#if TESTING_CAL_CONTAINMENT
  double Etotal = 0;
#endif

	// Return quickly if nothing in BCAL
   if ((cellCount == 0) && (showerCount == 0) && (spectraCount==0))
   {
      return (s_BarrelEMcal_t*)HDDM_NULL;
   }

	// Create HDDM structures to hang on output tree
   box = make_s_BarrelEMcal();
   box->bcalCells = make_s_BcalCells(cellCount);
	box->bcalSiPMSpectrums = make_s_BcalSiPMSpectrums(spectraCount);
   box->bcalTruthShowers = make_s_BcalTruthShowers(showerCount);

	// Loop over items stored in barrelEMcalTree. This includes
	// bcalCells items and bcalTruthShowers, but not bcalSiPMSpectrum
	// (Those are handled below.)
   while ( (item = (s_BarrelEMcal_t*) pickTwig(&barrelEMcalTree)) )
   {

#if WRITE_OUT_BCAL_CELL_HITS
		// ......................................................
		// This section is disabled since it corresponds to the
		// disabled section in hitBarrelEMcal above.
		
		// bcalCells (deprecated)
      s_BcalCells_t* cells = item->bcalCells;
      int cell;
      for (cell=0; cell < (int)cells->mult; ++cell)
      {
	 int m = box->bcalCells->mult;
         int mok = 0;

         s_BcalHits_t* hits = cells->in[cell].bcalHits;
          
         /* compress out the hits below threshold */
         int i,iok;
         for (iok=i=0; i < (int)hits->mult; i++)
         {
            if (hits->in[i].E >= THRESH_MEV/1e3)
            {
#if TESTING_CAL_CONTAINMENT
  Etotal += hits->in[i].E;
#endif
               if (iok < i)
               {
                  hits->in[iok] = hits->in[i];
               }
               ++iok;
               ++mok;
            }
         }
         if (hits != HDDM_NULL)
         {
            hits->mult = iok;
            if (iok == 0)
            {
               cells->in[cell].bcalHits = (s_BcalHits_t*)HDDM_NULL;
               FREE(hits);
            }
         }

         if (mok)
         {
            box->bcalCells->in[m] = cells->in[cell];
            box->bcalCells->mult++;
         }
      }
      if (cells != HDDM_NULL)
      {
         FREE(cells);
      }
		// ......................................................
#endif //WRITE_OUT_BCAL_CELL_HITS

		// bcalTruthShowers
      s_BcalTruthShowers_t* showers = item->bcalTruthShowers;
      int shower;
      for (shower=0; shower < (int)showers->mult; ++shower)
      {
         int m = box->bcalTruthShowers->mult++;
         box->bcalTruthShowers->in[m] = showers->in[shower];
      }
      if (showers != HDDM_NULL)
      {
         FREE(showers);
      }

      FREE(item);
   }

#if WRITE_OUT_BCAL_TIME_SPECTRA

	// Copy incident particle information
	box->bcalIncidentParticles = make_s_BcalIncidentParticles(BCAL_INCIDENT_PARTICLES.size());
	box->bcalIncidentParticles->mult=0;
	for(unsigned int i=0; i<BCAL_INCIDENT_PARTICLES.size(); i++){
		s_BcalIncidentParticle_t *iphddm = &box->bcalIncidentParticles->in[box->bcalIncidentParticles->mult++];
		IncidentParticle_t &ip = BCAL_INCIDENT_PARTICLES[i];
		iphddm->id = i+1;
		iphddm->ptype = ip.ptype;
		iphddm->x = ip.x;
		iphddm->y = ip.y;
		iphddm->z = ip.z;
		iphddm->px = ip.px;
		iphddm->py = ip.py;
		iphddm->pz = ip.pz;
	}

	// Sparsely copy timing spectra from local variables into HDDM structure
	map<bcal_index, DHistogram*>::iterator iter;
	for(iter=SiPMspectra.begin(); iter!=SiPMspectra.end(); iter++){

		s_BcalSiPMSpectrum_t *spectrum = &box->bcalSiPMSpectrums->in[box->bcalSiPMSpectrums->mult];
		const bcal_index &idx = iter->first;
		
		// Find index of first and last non-zero bins
		DHistogram *h = iter->second;
		
		// Fill in sparse data spectrum
		int bin_start = h->FindFirstNonZeroBin();
		int bin_end = h->FindLastNonZeroBin();
		if(bin_end>=bin_start && bin_end>0){
			spectrum->tstart = h->GetBinLowEdge(bin_start);
			spectrum->bin_width = h->GetBinWidth();
			spectrum->module = idx.module;
			spectrum->layer = idx.layer;
			spectrum->sector = idx.sector;
			spectrum->incident_id = idx.incident_id;
			spectrum->end = (idx.end==bcal_index::kUp ? 0:1);
			spectrum->Etruth = E_CellTruth[idx];
			string vals = "";
			double E_atten_sum = 0.0;
			for(int ibin=bin_start; ibin<=bin_end; ibin++){
				char str[256];
				double dE = h->GetBinContent(ibin);
				E_atten_sum += dE;
				if(dE!=0.0){
					sprintf(str, "%3.3f", dE*1000.0); // store in MeV since that cuts down on leading zeros
					vals += str;
				}else{
					vals += "0"; // lots of zeros are written and this saves space
				}
				vals += " ";
			}

			if(E_atten_sum >= THRESH_ATTENUATED_GEV){
				spectrum->vals = strdup(vals.c_str());
				box->bcalSiPMSpectrums->mult++;
			}
		}
		
		// Clear energy sum. 
		E_CellTruth[idx] = 0.0;
		
		// Reset the histo in preparation for the next event
		// and add it back to the pool
		iter->second->Reset();
		histo_pool.push_back(iter->second);
	}
	
	// Clear the sparsified list of histo objects with data
	SiPMspectra.clear();
	
	// Clear list of incident particles
	BCAL_INCIDENT_PARTICLES.clear();
	//_DBG_<<"BCAL_INCIDENT_PARTICLE_COUNT = "<<BCAL_INCIDENT_PARTICLE_COUNT<<" pool size:"<<histo_pool.size()<<endl;
	BCAL_INCIDENT_PARTICLE_COUNT = 0;
	INCIDENT_ID.assign(INCIDENT_ID.size(), 0); // don't release memory so we save reallocation in next event
	
	// Clear warning message flag so it can be shown on the next event
	SHOWED_INCIDENT_PARTICLE_SHORT_WARNING = false;
	
#endif //WRITE_OUT_BCAL_TIME_SPECTRA

	// Reset counters for next event
   cellCount = showerCount = 0;

	// Check if event is empty (after having applied thresholds above)
	// and remove branches that contain no data.
   if ((box->bcalCells != HDDM_NULL) &&
       (box->bcalCells->mult == 0))
   {
      FREE(box->bcalCells);
      box->bcalCells = (s_BcalCells_t*)HDDM_NULL;
   }
   if ((box->bcalSiPMSpectrums != HDDM_NULL) &&
       (box->bcalSiPMSpectrums->mult == 0))
   {
      FREE(box->bcalSiPMSpectrums);
      box->bcalSiPMSpectrums = (s_BcalSiPMSpectrums_t*)HDDM_NULL;
   }
   if ((box->bcalTruthShowers != HDDM_NULL) &&
       (box->bcalTruthShowers->mult == 0))
   {
      FREE(box->bcalTruthShowers);
      box->bcalTruthShowers = (s_BcalTruthShowers_t*)HDDM_NULL;
   }
   if ((box->bcalCells->mult == 0) &&
       (box->bcalSiPMSpectrums->mult == 0) &&
		 (box->bcalTruthShowers->mult == 0))
   {
      FREE(box);
      box = (s_BarrelEMcal_t*)HDDM_NULL;
   }
#if TESTING_CAL_CONTAINMENT
  printf("BCal energy sum: %f\n",Etotal);
#endif
   return box;
}
