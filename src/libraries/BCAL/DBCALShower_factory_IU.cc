/*
 *  DBCALShower_factory_IU.cc
 *  (formerly DBCALShower_factory.cc)
 *
 *  Created by Matthew Shepherd on 3/24/11.
 *
 */

#include "DBCALShower_factory_IU.h"
#include "DBCALCluster.h"

#include "DANA/DApplication.h"

#include "units.h"

#include "TH2F.h"
#include "TROOT.h"
#include "TDirectory.h"

#include <thread>

#include <DMatrix.h>
#include <DMatrixDSym.h>


DBCALShower_factory_IU::DBCALShower_factory_IU(){
	// defaults set to minimize probles if new values are not loaded
	LOAD_CCDB_CONSTANTS = 1.;
	energy_cutoff = 100;         ///< default to use linear part for all showers
	linear_intercept = 1;        ///< default to make no change to energy
	linear_slope = 0;            ///< default to make no change to energy
	exponential_param0 = 1;      ///< default to make no change to energy
	exponential_param1 = -10000; ///< default to make no change to energy
	exponential_param2 = 0;      ///< default to make no change to energy
	VERBOSE = 0;              ///< >0 once off info ; >2 event by event ; >3 everything
	COVARIANCEFILENAME = "";  ///<  Setting the filename will take precidence over the CCDB.  Files must end in ij.txt, where i and j are integers corresponding to the element of the matrix
	
	if (gPARMS){
		gPARMS->SetDefaultParameter("BCAL:LOAD_NONLIN_CCDB", LOAD_CCDB_CONSTANTS);
		/// use to set energy corrections on command line
		gPARMS->SetDefaultParameter("BCAL:energy_cutoff", energy_cutoff);
		gPARMS->SetDefaultParameter("BCAL:linear_slope", linear_slope);
		gPARMS->SetDefaultParameter("BCAL:linear_intercept", linear_intercept);
		gPARMS->SetDefaultParameter("BCAL:exponential_param0", exponential_param0);
		gPARMS->SetDefaultParameter("BCAL:exponential_param1", exponential_param1);
		gPARMS->SetDefaultParameter("BCAL:exponential_param2", exponential_param2);
		gPARMS->SetDefaultParameter("DBCALShower:VERBOSE", VERBOSE, "Verbosity level for DBCALShower objects and factories");
		gPARMS->SetDefaultParameter("DBCALShower:COVARIANCEFILENAME", COVARIANCEFILENAME, "File name for covariance files");
	}
}

jerror_t DBCALShower_factory_IU::brun(JEventLoop *loop, int32_t runnumber) {
  DApplication* app = dynamic_cast<DApplication*>(loop->GetJApplication());
  DGeometry* geom = app->GetDGeometry(runnumber);
  geom->GetTargetZ(m_zTarget);

    //by default, energy correction parameters are obtained through ccdb
    if(LOAD_CCDB_CONSTANTS > 0.5){
		map<string, double> shower_calib;
		if (loop->GetCalib("BCAL/shower_calib", shower_calib)) {
			jerr << "Error loading from CCDB\n";
		} else {
			energy_cutoff = shower_calib["energy_cutoff"];
			linear_intercept = shower_calib["linear_intercept"];
			linear_slope = shower_calib["linear_slope"];
			exponential_param0 = shower_calib["exponential_param0"];
			exponential_param1 = shower_calib["exponential_param1"];
			exponential_param2 = shower_calib["exponential_param2"];
		}
	}
	if (VERBOSE>0) {
		printf("%20s = %f\n","energy_cutoff",energy_cutoff);
		printf("%20s = %f\n","linear_intercept",linear_intercept);
		printf("%20s = %f\n","linear_slope",linear_slope);
		printf("%20s = %f\n","exponential_param0",exponential_param0);
		printf("%20s = %f\n","exponential_param1",exponential_param1);
		printf("%20s = %f\n","exponential_param2",exponential_param2);
	}
	
	jerror_t result = LoadCovarianceLookupTables();
	if (result!=NOERROR) return result;
	
  return NOERROR;
}

jerror_t
DBCALShower_factory_IU::evnt( JEventLoop *loop, uint64_t eventnumber ){
 
  vector< const DBCALCluster* > clusters;
  loop->Get( clusters );
  
  // loop through and fill the shower structure from the cluster
  // right now just a simple 1 to 1 correspondence with 
  // an overall energy correction
  
  for( vector< const DBCALCluster* >::const_iterator clItr = clusters.begin();
       clItr != clusters.end();
      ++clItr ){
   
    float cosTh = cos( (**clItr).theta() );
    float sinTh = sin( (**clItr).theta() );
    float cosPhi = cos( (**clItr).phi() );
    float sinPhi = sin( (**clItr).phi() );
    float rho = (**clItr).rho();
	if (VERBOSE>2) printf("cluster:   E=%f  th=%f phi=%f rho=%f   t=%f\n",
						  (**clItr).E(),(**clItr).theta()/3.14159265*180,(**clItr).phi()/3.14159265*180,(**clItr).rho(),(**clItr).t());

    DBCALShower* shower = new DBCALShower();
    
    shower->E_raw = (**clItr).E();
    shower->E_preshower = (**clItr).E_preshower();
    shower->x = rho * sinTh * cosPhi;
    shower->y = rho * sinTh * sinPhi;
    shower->z = rho * cosTh + m_zTarget;

    //DBCALCluster::t() returns the time at the inner radius
    //so we need to make an adjustment so that the shower t is the time at
    //the shower location (x,y,z)
    double t = (**clItr).t();
    double inner_rad = DBCALGeometry::GetBCAL_inner_rad();
    double dist_in_BCAL = rho - inner_rad/sinTh;
    t = t + dist_in_BCAL/(30*k_cm/k_nsec);
    shower->t = t;

    shower->N_cell = (**clItr).nCells();
    
    // non-linear energy corrections can be found at https://logbooks.jlab.org/entry/3419524 
    
    if( shower->E_raw < energy_cutoff ) shower->E = shower->E_raw / (linear_intercept + linear_slope * shower->E_raw ) ;
    if( shower->E_raw >= energy_cutoff ) shower->E = shower->E_raw / (exponential_param0 - exp(exponential_param1 * shower->E_raw + exponential_param2));

	// Get covariance matrix and uncertainties
	FillCovarianceMatrix(shower);
	if (VERBOSE>2) {
		printf("shower:    E=%f   x=%f   y=%f   z=%f   t=%f\n",
			   shower->E,shower->x,shower->y,shower->z,shower->t);
		printf("shower:   dE=%f  dx=%f  dy=%f  dz=%f  dt=%f\n",
			   shower->EErr(),shower->xErr(),shower->yErr(),shower->zErr(),shower->tErr());
		printf("shower:   Ex=%f  Ey=%f  Ez=%f  Et=%f  xy=%f\n",
			   shower->EXcorr(),shower->EYcorr(),shower->EZcorr(),shower->ETcorr(),shower->XYcorr());
		printf("shower:   xz=%f  xt=%f  yz=%f  yt=%f  zt=%f\n",
			   shower->XZcorr(),shower->XTcorr(),shower->YZcorr(),shower->YTcorr(),shower->ZTcorr());
	}

    shower->AddAssociatedObject(*clItr);
    
    _data.push_back( shower );
  }
  
  return NOERROR;
}


jerror_t
DBCALShower_factory_IU::FillCovarianceMatrix(DBCALShower *shower){
	/// This function takes a BCALShower object and using the internal variables 
	/// overwrites any existing covaraince matrix using lookup tables. 

	// Get edges of lookup table histograms (assume that all histograms have the same limits.)
	TAxis *xaxis = CovarianceLookupTable[0][0]->GetXaxis();
	TAxis *yaxis = CovarianceLookupTable[0][0]->GetYaxis();
	float minElookup = xaxis->GetBinLowEdge(1);
	float maxElookup = xaxis->GetBinUpEdge(xaxis->GetNbins());
	float minthlookup = yaxis->GetBinLowEdge(1);
	float maxthlookup = yaxis->GetBinUpEdge(yaxis->GetNbins());

	float shower_r = sqrt(shower->x*shower->x + shower->y*shower->y);
	float shower_theta = atan2(shower_r,shower->z);
	float thlookup = shower_theta/3.14159265*180;
	float Elookup = shower->E;

	// Adjust values: in order to use Interpolate() must be within histogram range
	if (Elookup<minElookup) Elookup=minElookup;
	if (Elookup>maxElookup) Elookup=maxElookup-0.0001; // move below edge, on edge doesn't work.
	if (thlookup<minthlookup) thlookup=minthlookup;
	if (thlookup>maxthlookup) thlookup=maxthlookup-0.0001;
	if (VERBOSE>3) printf("(%f,%F)    limits (%f,%f)  (%f,%f)\n",Elookup,thlookup,minElookup,maxElookup,minthlookup,maxthlookup);

	DMatrixDSym ErphiztCovariance(5);
	for (int i=0; i<5; i++) {
		for (int j=0; j<=i; j++) {
			float val = CovarianceLookupTable[i][j]->Interpolate(Elookup, thlookup);
			//if (i==2) val*=shower_r; // convert phi to phihat
			//if (j==2) val*=shower_r; // convert phi to phihat
			ErphiztCovariance(i,j) = ErphiztCovariance(j,i) = val;
		}
	}
	
	float shower_phi = atan2(shower->y,shower->x);
    float cosPhi = cos(shower_phi);
    float sinPhi = sin(shower_phi);
	DMatrix rotationmatrix(5,5);
	rotationmatrix(0,0) = 1;
	rotationmatrix(3,3) = 1;
	rotationmatrix(4,4) = 1;
	rotationmatrix(1,1) = cosPhi;
	rotationmatrix(1,2) = -sinPhi;
	rotationmatrix(2,1) = sinPhi;
	rotationmatrix(2,2) = cosPhi;

	if (VERBOSE>3) {printf("(E,r,phi,z,t)  "); ErphiztCovariance.Print(); }
	DMatrixDSym &D = ErphiztCovariance.Similarity(rotationmatrix);
	shower->ExyztCovariance = D;
	if (VERBOSE>2) {printf("(E,x,y,z,t)    "); shower->ExyztCovariance.Print(); }

	return NOERROR;
}


jerror_t
DBCALShower_factory_IU::LoadCovarianceLookupTables(){
	std::thread::id this_id = std::this_thread::get_id();
	stringstream idstring;
	idstring << this_id;
	if (VERBOSE>0) printf("DBCALShower_factory_IU::LoadCovarianceLookupTables():  Thread %s\n",idstring.str().c_str());

	bool USECCDB=0;
	// if filename specified try to use filename else get info from CCDB
	if (COVARIANCEFILENAME == "") USECCDB=1;
	
	map<string,string> covariance_data;	
	if (USECCDB) {
		// load information for covariance matrix
		if (eventLoop->GetJCalibration()->GetCalib("/BCAL/shower_covariance", covariance_data)) {
			jerr << "Error loading /BCAL/shower_covariance !" << endl;
			return ERROR_OPENING_EVENT_SOURCE;
		}
		if (covariance_data.size() == 15)  {  // there are 15 elements in the covariance matrix
			// for example, print it all out
			if (VERBOSE>0) {				
				for(auto element : covariance_data) {
					cout << "\nTEST:   " << element.first << " = " << element.second << endl;
				}
			}
		} else {
			jerr << "Wrong number of elements /BCAL/shower_covariance !" << endl;
			return ERROR_OPENING_EVENT_SOURCE;	
		}
	}

	for (int i=0; i<5; i++) {
		for (int j=0; j<=i; j++) {

			japp->RootWriteLock();
			// change directory to memory so that histograms are not saved to file
			TDirectory *savedir = gDirectory;
			gROOT->cd();  

			// Read in string
			ifstream ifs;
			string line;
			stringstream ss;
			if (USECCDB) {
				stringstream matrixname; 
				matrixname << "covmatrix_" << i << j;
				if (VERBOSE>1) cout << "Using CCDB \"" << matrixname.str() << "\"  " << covariance_data[matrixname.str()] << endl;
				ss.str(covariance_data[matrixname.str()]);
			} else {
				char filename[255];
				sprintf(filename,"%s%i%i.txt",COVARIANCEFILENAME.c_str(),i,j);
				if (VERBOSE>0) cout  << filename << std::endl;
				ifs.open(filename);
				if (! ifs.is_open()) {
					jerr << " Error: Cannot open file! " << filename << std::endl;
					return ERROR_OPENING_EVENT_SOURCE;
				}
				getline(ifs, line, '\n');
				ss.str(line);
				if (VERBOSE>1) cout << "Using file " <<line<<endl;
			}

			// Parse string 
			int nxbins, nybins;
			ss>>nxbins;
			ss>>nybins;
			if (VERBOSE>1) printf("bins (%i,%i)\n",nxbins,nybins);
			Float_t xbins[nxbins+1];
			Float_t ybins[nybins+1];
			for (int count=0; count<=nxbins; count++) {
				ss>>xbins[count];
				if (VERBOSE>1) printf("(%i,%f)  ",count,xbins[count]);
			}
			if (VERBOSE>1) printf("\n");
			for (int count=0; count<=nybins; count++) {
				ss>>ybins[count];
				if (VERBOSE>1) printf("(%i,%f)  ",count,ybins[count]);
			}
			if (VERBOSE>1) printf("\n");
			int xbin=1;
			double cont;
			int ybin=1;
			// create histogram
			char histname[255];
			sprintf(histname,"covariance_%i%i_thread%s",i,j,idstring.str().c_str());
			CovarianceLookupTable[i][j] = new TH2F(histname,"Covariance histogram",nxbins,xbins,nybins,ybins);
			// fill histogram
			while(ss>>cont){
				if (VERBOSE>1) printf("(%2i,%2i)  (%2i,%2i)  %e  ",i,j,xbin,ybin,cont);
				CovarianceLookupTable[i][j]->SetBinContent(xbin,ybin,cont);
				ybin++;
				if (ybin>nybins) { xbin++; ybin=1; }
			}
			if (VERBOSE>1) printf("\n");
			// Close file
			ifs.close();

			savedir->cd();
			japp->RootUnLock(); 

		}
	}
	return NOERROR;
}

