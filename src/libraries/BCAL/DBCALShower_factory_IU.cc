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

	VERBOSE = 0; // >0 once off info ; >2 event by event
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
	if (VERBOSE>3) printf("cluster: E=%f th=%f phi=%f rho=%f t=%f\n",
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

	if (VERBOSE>3) printf("shower:  E=%f x=%f y=%f z=%f t=%f\n",
						  shower->E,shower->x,shower->y,shower->z,shower->t);
	// Get covariance matrix and uncertainties

	// Get histogram edges. Assuming that all histograms have the same limits
	TAxis *xaxis = CovarEmntLookup[0][0]->GetXaxis();
	TAxis *yaxis = CovarEmntLookup[0][0]->GetYaxis();
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
			float val = CovarEmntLookup[i][j]->Interpolate(Elookup, thlookup);
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

	for (int i=0; i<5; i++) {
		for (int j=0; j<=i; j++) {

			japp->RootWriteLock();
			// change directory to memory so that histograms are not saved to file
			TDirectory *savedir = gDirectory;
			gROOT->cd();  


			char filename[255];
			sprintf(filename,"covmatrix/covmatrix_%i%i_bcal.txt",i,j);
			if (VERBOSE>0) cout  << filename << std::endl;
			ifstream ifs(filename);
			if (! ifs.is_open()) {
				cout << " Error: Cannot open file! " << filename << std::endl;
				return ERROR_OPENING_EVENT_SOURCE;
			}

			// Parse string and create histogram
			char histname[255];
			sprintf(histname,"covariance_%i%i_thread%s",i,j,idstring.str().c_str());
			int nxbins, nybins;
			
			string line;
			getline(ifs, line, '\n');
			if (VERBOSE>1) cout<<line<<endl;
			stringstream ss1(line);
			ss1>>nxbins;
			getline(ifs, line, '\n');
			if (VERBOSE>1) cout<<line<<endl;
			stringstream ss2(line);
			ss2>>nybins;
			if (VERBOSE>1) printf("bins %i %i\n",nxbins,nybins);

			Float_t xbins[nxbins];
			Float_t ybins[nybins];
			int count=0;
			getline(ifs, line, '\n');
			if (VERBOSE>1) cout<<line<<endl;
			stringstream ss3(line);
			while(ss3>>xbins[count]){
				if (VERBOSE>1) printf("%i  %f\n",count,xbins[count]);
				count++;
			}
			count=0;
			getline(ifs, line, '\n');
			if (VERBOSE>1) cout<<line<<endl;
			stringstream ss4(line);
			while(ss4>>ybins[count]){
				if (VERBOSE>1) printf("%i  %f\n",count,ybins[count]);
				count++;
			}
			CovarEmntLookup[i][j] = new TH2F(histname,"Covariance histogram",nxbins,xbins,nybins,ybins);
			CovarEmntLookup[i][j]->Draw();
			int xbin=0;
			while(getline(ifs, line, '\n')){
				if (VERBOSE>1) cout<<line<<endl;
				stringstream ss(line);
				double cont;
				xbin++;
				int ybin=0;
				while(ss>>cont){
					if (VERBOSE>1) cout<<cont<<endl;
					ybin++;
					if (VERBOSE>1) printf("%2i  %2i  %2i  %2i  %f\n",i,j,xbin,ybin,cont);
					CovarEmntLookup[i][j]->SetBinContent(xbin,ybin,cont);
				}
			}

			// Close file
			ifs.close();

			savedir->cd();
			japp->RootUnLock(); 

		}
	}
	return NOERROR;
}
