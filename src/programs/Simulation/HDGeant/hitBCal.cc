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
#define WRITE_OUT_BCAL_TIME_SPECTRA 0


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

extern s_HDDM_t* thisInputEvent;

// Prevent name mangling so these routines keep their
// C-style names in the object
extern "C"{
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
		
		bcal_index(unsigned int module, unsigned int layer, unsigned int sector, EndType end):
			module(module),layer(layer),sector(sector),end(end){}
	
		unsigned int module;
		unsigned int layer;
		unsigned int sector;
		EndType end;
		
		bool operator<(const bcal_index &idx) const{
			if(module<idx.module)return true;
			if(module>idx.module)return false;
			if(layer<idx.layer)return true;
			if(layer>idx.layer)return false;
			if(sector<idx.sector)return true;
			if(sector>idx.sector)return false;
			if((end==kUp) && (idx.end==kDown))return true;
			return false;
		}
};

// This holds the energy-weighted time spectra for
// each SiPM in the event. It is created and initialized
// on the first call to hitBarrelEMcal
map<bcal_index, DHistogram*> SiPMspectra;

// We keep a list of the histograms that need to be 
// cleared at the end of the event so that we don't
// have to wasted time clearing histos that were never
// filled.
map<bcal_index, DHistogram*> histos_with_data;

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
	// Create time histos for each SiPM
	//
	// We want the histo to cover the range specified by the BGGATE
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
	int Nbins = (int)floor(0.5 + (BGGATE2 - BGGATE1)/SiPM_tbin_width);
	BGGATE2 = BGGATE1 + SiPM_tbin_width*(double)Nbins;
	for(unsigned int imodule=1; imodule<=48; imodule++){
	  for(unsigned int ilayer=1; ilayer<=10; ilayer++){
	 	  for(unsigned int isector=1; isector<=4; isector++){
			  bcal_index idxUp(imodule, ilayer, isector, bcal_index::kUp);
			  bcal_index idxDn(imodule, ilayer, isector, bcal_index::kDown);

			  SiPMspectra[idxUp] = new DHistogram(Nbins, BGGATE1, BGGATE2);
			  SiPMspectra[idxDn] = new DHistogram(Nbins, BGGATE1, BGGATE2);
			  
			  E_CellTruth[idxUp] = 0.0;
			  E_CellTruth[idxDn] = 0.0;
		  }
	  }
	}
	cout<<"Created "<<SiPMspectra.size()<<" histograms for BCAL timing spectra"<<endl;
#endif  //WRITE_OUT_BCAL_TIME_SPECTRA

   initialized = 1;
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
      int sector = getsector_();
      int layer  = getlayer_();
      int module = getmodule_();
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
		// Get map index based on layer and sector
      unsigned int sector = getsector_();
      unsigned int layer  = getlayer_();
      unsigned int module = getmodule_();
		bcal_index idxUp(module, layer, sector, bcal_index::kUp);
		bcal_index idxDn(module, layer, sector, bcal_index::kDown);
		map<bcal_index, DHistogram*>::iterator iterUp = SiPMspectra.find(idxUp);
		map<bcal_index, DHistogram*>::iterator iterDn = SiPMspectra.find(idxDn);
		
		if(iterUp==SiPMspectra.end() || iterDn==SiPMspectra.end()){
			cout<<"Not found: module:"<<module<<" layer:"<<layer<<" sector:"<<sector<<endl;
		}else{
			double Z = xlocal[2];
			double dist_up = 390.0/2.0 + Z;
			double dist_dn = 390.0/2.0 - Z;
			
			double Eup = dEsum*exp(-dist_up/ATTEN_LENGTH);
			double Edn = dEsum*exp(-dist_dn/ATTEN_LENGTH);
			double t_up = t + dist_up/C_EFFECTIVE;
			double t_dn = t + dist_dn/C_EFFECTIVE;
			
			DHistogram *hup = iterUp->second;
			DHistogram *hdn = iterDn->second;
			
			hup->Fill(t_up, Eup);
			hdn->Fill(t_dn, Edn);
			
			E_CellTruth[idxUp] += dEsum;
			E_CellTruth[idxDn] += dEsum;
			
			histos_with_data[idxUp] = hup;
			histos_with_data[idxDn] = hdn;
		}
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
	unsigned int spectraCount = histos_with_data.size(); // Number non-zero spectra
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
	// Sparsely copy timing spectra from local variables into HDDM structure
	map<bcal_index, DHistogram*>::iterator iter;
	for(iter=histos_with_data.begin(); iter!=histos_with_data.end(); iter++){
		
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
		iter->second->Reset();
	}
	
	// Clear the sparsified list of histo objects with data
	histos_with_data.clear();
	
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
