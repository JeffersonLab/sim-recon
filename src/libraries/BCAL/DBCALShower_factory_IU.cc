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

DBCALShower_factory_IU::DBCALShower_factory_IU(){

	VERBOSE = 0; // >0 once off info ; >2 event by event ; >3 everything
	COVARIANCEFILENAME = ""; 
	// setting the filename will take precidence over the CCDB.  
	// Files must end in ij.txt, where i and j are integers corresponding to the element of the matrix
	if (gPARMS){
		gPARMS->SetDefaultParameter("DBCALShower:VERBOSE", VERBOSE, "Verbosity level for DBCALShower objects and factories");
		gPARMS->SetDefaultParameter("DBCALShower:COVARIANCEFILENAME", COVARIANCEFILENAME, "File name for covariance files");
	}

}

jerror_t DBCALShower_factory_IU::brun(JEventLoop *loop, int32_t runnumber) {
  DApplication* app = dynamic_cast<DApplication*>(loop->GetJApplication());
  DGeometry* geom = app->GetDGeometry(runnumber);
  geom->GetTargetZ(m_zTarget);

  jerror_t result = CreateCovarianceMatrix();
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
	float r = rho * sinTh;
	if (VERBOSE>2) printf("cluster: E=%f th=%f phi=%f rho=%f t=%f\n",
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
    
    shower->E = shower->E_raw;

	if (VERBOSE>2) printf("shower:  E=%f x=%f y=%f z=%f t=%f\n",
						  shower->E,shower->x,shower->y,shower->z,shower->t);
	// Get covariance matrix and uncertainties

	// Get histogram edges. Assuming that all histograms have the same limits
	TAxis *xaxis = CovarElementLookup[0][0]->GetXaxis();
	TAxis *yaxis = CovarElementLookup[0][0]->GetYaxis();
	float minElookup = xaxis->GetBinLowEdge(1);
	float maxElookup = xaxis->GetBinUpEdge(xaxis->GetNbins());
	float minthlookup = yaxis->GetBinLowEdge(1);
	float maxthlookup = yaxis->GetBinUpEdge(yaxis->GetNbins());
	float thlookup = (**clItr).theta()/3.14159265*180;
	float Elookup = shower->E;
	// Adjust values in order to use Interpolate()
	if (Elookup<minElookup) Elookup=minElookup;
	if (Elookup>maxElookup) Elookup=maxElookup-0.0001; // move below edge, on edge doesn't work.
	if (thlookup<minthlookup) thlookup=minthlookup;
	if (thlookup>maxthlookup) thlookup=maxthlookup-0.0001;
	if (VERBOSE>3) printf("(%f,%F)    limits (%f,%f)  (%f,%f)\n",Elookup,thlookup,minElookup,maxElookup,minthlookup,maxthlookup);

	DMatrix ErphiztCovariance(5,5);
	for (int i=0; i<5; i++) {
		for (int j=0; j<=i; j++) {
			float val = CovarElementLookup[i][j]->Interpolate(Elookup, thlookup);
			if (i==2) val*=r; // convert phi to phihat
			if (j==2) val*=r; // convert phi to phihat
			ErphiztCovariance(i,j) = ErphiztCovariance(j,i) = val;
		}
	}
	DMatrix rotationmatrix(5,5);
	rotationmatrix(0,0) = 1;
	rotationmatrix(3,3) = 1;
	rotationmatrix(4,4) = 1;
	rotationmatrix(1,1) = cosPhi;
	rotationmatrix(1,2) = -sinPhi;
	rotationmatrix(2,1) = sinPhi;
	rotationmatrix(2,2) = cosPhi;
	DMatrix rotationmatrixT(DMatrix::kTransposed,rotationmatrix);

	if (VERBOSE>3) {printf("(E,r,phi,z,t)  "); ErphiztCovariance.Print(); }
	shower->ExyztCovariance.ResizeTo(5,5);
	shower->ExyztCovariance = rotationmatrix*ErphiztCovariance*rotationmatrixT;  
	if (VERBOSE>2) {printf("(E,x,y,z,t)    "); shower->ExyztCovariance.Print(); }

	// Redundant duplication of uncertainties.  Should be removed
	shower->xErr = sqrt(shower->ExyztCovariance(1,1));
	shower->yErr = sqrt(shower->ExyztCovariance(2,2));
	shower->zErr = sqrt(shower->ExyztCovariance(3,3));
	shower->tErr = sqrt(shower->ExyztCovariance(4,4));

    shower->AddAssociatedObject(*clItr);
    
    _data.push_back( shower );
  }
  
  return NOERROR;
}




jerror_t
DBCALShower_factory_IU::CreateCovarianceMatrix(){
	std::thread::id this_id = std::this_thread::get_id();
	stringstream idstring;
	idstring << this_id;
	if (VERBOSE>0) printf("DBCALShower_factory_IU::CreateCovarianceMatrix():  Thread %s\n",idstring.str().c_str());

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
				if (VERBOSE>1) printf("%i  %f\n",count,xbins[count]);
			}
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
			CovarElementLookup[i][j] = new TH2F(histname,"Covariance histogram",nxbins,xbins,nybins,ybins);
			// fill histogram
			while(ss>>cont){
				if (VERBOSE>1) printf("(%2i,%2i)  (%2i,%2i)  %e\n",i,j,xbin,ybin,cont);
				CovarElementLookup[i][j]->SetBinContent(xbin,ybin,cont);
				ybin++;
				if (ybin>nybins) { xbin++; ybin=1; }
			}
			// Close file
			ifs.close();

			savedir->cd();
			japp->RootUnLock(); 

		}
	}
	return NOERROR;
}

