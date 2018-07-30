//
//    File: DFCALShower_factory.cc
// Created: Tue May 17 11:57:50 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)

#include <thread>
#include <math.h>
#include <DVector3.h>
#include "TH2F.h"
#include "TROOT.h"
#include "TDirectory.h"
using namespace std;

#include "FCAL/DFCALShower_factory.h"
#include "FCAL/DFCALGeometry.h"
#include "FCAL/DFCALCluster.h"
#include "FCAL/DFCALHit.h"
#include "TRACKING/DTrackWireBased.h"
#include <JANA/JEvent.h>
#include <JANA/JApplication.h>
using namespace jana;

//----------------
// Constructor
//----------------
DFCALShower_factory::DFCALShower_factory()
{
  // should we use CCDB constants?
  LOAD_CCDB_CONSTANTS = 1.;
  gPARMS->SetDefaultParameter("FCAL:LOAD_NONLIN_CCDB", LOAD_CCDB_CONSTANTS);

  SHOWER_ENERGY_THRESHOLD = 50*k_MeV;
  gPARMS->SetDefaultParameter("FCAL:SHOWER_ENERGY_THRESHOLD", SHOWER_ENERGY_THRESHOLD);

  // these need to come from database to ensure accuracy
  // remove default value which might be close to the right solution,
  // but not quite correct -- allow command line tuning

  cutoff_energy= 0;
  linfit_slope = 0;
  linfit_intercept = 0;
  expfit_param1 = 0;
  expfit_param2 = 0;
  expfit_param3 = 0;
	
  timeConst0 = 0;
  timeConst1 = 0; 
  timeConst2 = 0;
  timeConst3 = 0; 
  timeConst4 = 0;

  gPARMS->SetDefaultParameter("FCAL:cutoff_enegry", cutoff_energy);
  gPARMS->SetDefaultParameter("FCAL:linfit_slope", linfit_slope);
  gPARMS->SetDefaultParameter("FCAL:linfit_intercept", linfit_intercept);
  gPARMS->SetDefaultParameter("FCAL:expfit_param1", expfit_param1);
  gPARMS->SetDefaultParameter("FCAL:expfit_param2", expfit_param2);
  gPARMS->SetDefaultParameter("FCAL:expfit_param3", expfit_param3);

  gPARMS->SetDefaultParameter("FCAL:P0", timeConst0);
  gPARMS->SetDefaultParameter("FCAL:P1", timeConst1);
  gPARMS->SetDefaultParameter("FCAL:P2", timeConst2);
  gPARMS->SetDefaultParameter("FCAL:P3", timeConst3);
  gPARMS->SetDefaultParameter("FCAL:P4", timeConst4);

  // Parameters to make shower-depth correction taken from Radphi, 
  // slightly modifed to match photon-polar angle
  FCAL_RADIATION_LENGTH = 3.1;
  FCAL_CRITICAL_ENERGY = 0.035;
  FCAL_SHOWER_OFFSET = 1.0;
	
  gPARMS->SetDefaultParameter("FCAL:FCAL_RADIATION_LENGTH", FCAL_RADIATION_LENGTH);
  gPARMS->SetDefaultParameter("FCAL:FCAL_CRITICAL_ENERGY", FCAL_CRITICAL_ENERGY);
  gPARMS->SetDefaultParameter("FCAL:FCAL_SHOWER_OFFSET", FCAL_SHOWER_OFFSET);

  VERBOSE = 0;              ///< >0 once off info ; >2 event by event ; >3 everything
  COVARIANCEFILENAME = "";  ///<  Setting the filename will take precidence over the CCDB.  Files must end in ij.txt, where i and j are integers corresponding to the element of the matrix
  gPARMS->SetDefaultParameter("DFCALShower:VERBOSE", VERBOSE, "Verbosity level for DFCALShower objects and factories");
  gPARMS->SetDefaultParameter("DFCALShower:COVARIANCEFILENAME", COVARIANCEFILENAME, "File name for covariance files");

}

//------------------
// brun
//------------------
jerror_t DFCALShower_factory::brun(JEventLoop *loop, int32_t runnumber)
{
 
  // Get calibration constants
  map<string, double> fcal_parms;
  loop->GetCalib("FCAL/fcal_parms", fcal_parms);
  if (fcal_parms.find("FCAL_C_EFFECTIVE")!=fcal_parms.end()){
    FCAL_C_EFFECTIVE = fcal_parms["FCAL_C_EFFECTIVE"];
    if(debug_level>0)jout<<"FCAL_C_EFFECTIVE = "<<FCAL_C_EFFECTIVE<<endl;
  } else {
    jerr<<"Unable to get FCAL_C_EFFECTIVE from FCAL/fcal_parms in Calib database!"<<endl;
  }
  
  DApplication *dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
  const DGeometry *geom = dapp->GetDGeometry(runnumber);
    
  if (geom) {
    geom->GetTargetZ(m_zTarget);
    geom->GetFCALZ(m_FCALfront);
  }
  else{
      
    cerr << "No geometry accessbile." << endl;
    return RESOURCE_UNAVAILABLE;
  }

  // by default, load non-linear shower corrections from the CCDB
  // but allow these to be overridden by command line parameters
  if(LOAD_CCDB_CONSTANTS > 0.1) {
    map<string, double> shower_calib_piecewise;
    loop->GetCalib("FCAL/shower_calib_piecewise", shower_calib_piecewise);
    cutoff_energy = shower_calib_piecewise["cutoff_energy"];
    linfit_slope = shower_calib_piecewise["linfit_slope"];
    linfit_intercept = shower_calib_piecewise["linfit_intercept"];
    expfit_param1 = shower_calib_piecewise["expfit_param1"];
    expfit_param2 = shower_calib_piecewise["expfit_param2"];
    expfit_param3 = shower_calib_piecewise["expfit_param3"];

    if(debug_level>0) {
      jout << "cutoff_energy = " << cutoff_energy << endl;
      jout << "linfit_slope = " << linfit_slope << endl;
      jout << "linfit_intercept = " << linfit_intercept << endl;
      jout << "expfit_param1 = " << expfit_param1 << endl;
      jout << "expfit_param2 = " << expfit_param2<< endl;
      jout << "expfit_param3 = " << expfit_param3 << endl;

    }
  }

  // Get timing correction polynomial, J. Mirabelli 10/31/17
  if(LOAD_CCDB_CONSTANTS > 0.1) {
    map<string,double> timing_correction;
    loop->GetCalib("FCAL/shower_timing_correction", timing_correction); 
    timeConst0 = timing_correction["P0"];
    timeConst1 = timing_correction["P1"];     
    timeConst2 = timing_correction["P2"];
    timeConst3 = timing_correction["P3"];
    timeConst4 = timing_correction["P4"];

    if(debug_level>0) {

      jout << "timeConst0 = " << timeConst0 << endl;
      jout << "timeConst1 = " << timeConst1 << endl;
      jout << "timeConst2 = " << timeConst2 << endl;
      jout << "timeConst3 = " << timeConst3 << endl;
      jout << "timeConst4 = " << timeConst4 << endl;
    }

  }

	


  jerror_t result = LoadCovarianceLookupTables(eventLoop);
  if (result!=NOERROR) return result;

  return NOERROR;
}


jerror_t DFCALShower_factory::erun(void) {
  // delete lookup tables to prevent memory leak
  for (int i=0; i<5; i++) {
    for (int j=0; j<=i; j++) {
      delete CovarianceLookupTable[i][j];
      CovarianceLookupTable[i][j] = nullptr;
    }
  }
  return NOERROR;
}


//------------------
// evnt
//------------------
jerror_t DFCALShower_factory::evnt(JEventLoop *eventLoop, uint64_t eventnumber)
{
  vector<const DFCALCluster*> fcalClusters;
  eventLoop->Get(fcalClusters);
  if(fcalClusters.size()<1)return NOERROR;
 
  // Use the center of the target as an approximation for the vertex position
  DVector3 vertex(0.0, 0.0, m_zTarget);
  
  vector< const DTrackWireBased* > allWBTracks;
  eventLoop->Get( allWBTracks );
  vector< const DTrackWireBased* > wbTracks = filterWireBasedTracks( allWBTracks );

  // Loop over list of DFCALCluster objects and calculate the "Non-linear" corrected
  // energy and position for each. We'll use a logarithmic energy-weighting to 
  // find the final position and error. 
  for( vector< const DFCALCluster* >::const_iterator clItr = fcalClusters.begin();
       clItr != fcalClusters.end();  ++clItr ){
    const DFCALCluster* cluster=*clItr;

    // energy weighted time provides better resolution:
    double cTime = cluster->getTimeEWeight();

    double errZ;  // will be filled by call to GetCorrectedEnergyAndPosition()
		
    // Get corrected energy, position, and errZ
    double Ecorrected;
    DVector3 pos_corrected;
    GetCorrectedEnergyAndPosition( cluster , Ecorrected, pos_corrected, errZ, &vertex);

    if (Ecorrected>0.){		
      //up to this point, all times have been times at which light reaches
      //the back of the detector. Here we correct for the time that it 
      //takes the Cherenkov light to reach the back of the detector
      //so that the t reported is roughly the time of the shower at the
      //position pos_corrected	
      cTime -= ( m_FCALfront + DFCALGeometry::blockLength() - pos_corrected.Z() )/FCAL_C_EFFECTIVE;

      //Apply time-walk correction/global timing offset
      cTime += ( timeConst0  +  timeConst1 * Ecorrected  +  timeConst2 * TMath::Power( Ecorrected, 2 ) +
		 timeConst3 * TMath::Power( Ecorrected, 3 )  +  timeConst4 * TMath::Power( Ecorrected, 4 ) );

      // Make the DFCALShower object
      DFCALShower* shower = new DFCALShower;
      
      shower->setEnergy( Ecorrected );
      shower->setPosition( pos_corrected );   
      shower->setTime ( cTime );
      FillCovarianceMatrix( shower );

      if( VERBOSE > 2 ){
	printf("FCAL shower:    E=%f   x=%f   y=%f   z=%f   t=%f\n",
	       shower->getEnergy(),shower->getPosition().X(),shower->getPosition().Y(),shower->getPosition().Z(),shower->getTime());
	printf("FCAL shower:   dE=%f  dx=%f  dy=%f  dz=%f  dt=%f\n",
	       shower->EErr(),shower->xErr(),shower->yErr(),shower->zErr(),shower->tErr());
	printf("FCAL shower:   Ex=%f  Ey=%f  Ez=%f  Et=%f  xy=%f\n",
	       shower->EXcorr(),shower->EYcorr(),shower->EZcorr(),shower->ETcorr(),shower->XYcorr());
	printf("FCAL shower:   xz=%f  xt=%f  yz=%f  yt=%f  zt=%f\n",
	       shower->XZcorr(),shower->XTcorr(),shower->YZcorr(),shower->YTcorr(),shower->ZTcorr());
      }

      // now fill information related to shower shape and nearby
      // tracks -- useful for splitoff rejection later

      double docaTr = 1E6;
      double timeTr = 1E6;
      double xTr = 0;
      double yTr = 0;

      double flightTime;
      DVector3 projPos, projMom;

      // find the closest track to the shower -- here we loop over the best FOM
      // wire-based track for every track candidate not just the ones associated
      // with the topology
      for( size_t iTrk = 0; iTrk < wbTracks.size(); ++iTrk ){

	if( !wbTracks[iTrk]->GetProjection( SYS_FCAL, projPos, &projMom, &flightTime ) ) continue;
	
	// need to swim fcalPos to common z for DOCA calculation -- this really
	// shouldn't be in the loop if the z-value of projPos doesn't change
	// with each track
	
	DVector3 fcalFacePos = ( shower->getPosition() - vertex );
	fcalFacePos.SetMag( fcalFacePos.Mag() * projPos.Z() / fcalFacePos.Z() );
 
	double distance = ( fcalFacePos - projPos ).Mag();
	
	if( distance < docaTr ){

	  docaTr = distance;
	  // this is the time from the center of the target to the detector -- to compare with
	  // the FCAL time, one needs to have the t0RF at the center of the target.  That
	  // comparison happens at a later stage in the analysis.
	  timeTr = ( wbTracks[iTrk]->position().Z() - vertex.Z() ) / SPEED_OF_LIGHT + flightTime;
	  xTr = projPos.X();
	  yTr = projPos.Y();
	}
      }

      shower->setDocaTrack( docaTr );
      shower->setTimeTrack( timeTr );

      // now compute some variables at the hit level
      
      vector< const DFCALHit* > fcalHits;
      cluster->Get( fcalHits );
      shower->setNumBlocks( fcalHits.size() );
      
      double e9e25, e1e9;
      getE1925FromHits( e1e9, e9e25, fcalHits, getMaxHit( fcalHits ) );
      shower->setE1E9( e1e9 );
      shower->setE9E25( e9e25 );

      double sumU = 0;
      double sumV = 0;
      // if there is no nearest track, the defaults for xTr and yTr will result
      // in using the beam axis as the directional axis
      getUVFromHits( sumU, sumV, fcalHits,
		     DVector3( shower->getPosition().X(), shower->getPosition().Y(), 0 ),
		     DVector3( xTr, yTr, 0 ) );
      
      shower->setSumU( sumU );
      shower->setSumV( sumV );
      
      shower->AddAssociatedObject( cluster );

      _data.push_back(shower);
    }
  }

  return NOERROR;
}

//--------------------------------
// GetCorrectedEnergyAndPosition
//
// Non-linear and depth corrections should be fixed within DFCALShower member functions
//--------------------------------
void DFCALShower_factory::GetCorrectedEnergyAndPosition(const DFCALCluster* cluster, double &Ecorrected, DVector3 &pos_corrected, double &errZ, const DVector3 *vertex)
{
  // Non-linear energy correction are done here
  //int MAXITER = 1000;

  DVector3  posInCal = cluster->getCentroid();
  float x0 = posInCal.Px();
  float y0 = posInCal.Py();

  double Eclust = cluster->getEnergy();
  
  double Ecutoff = cutoff_energy;
  double A  = linfit_slope;
  double B  = linfit_intercept;
  double C  = expfit_param1;
  double D  = expfit_param2;
  double E  = expfit_param3;

	 
  double Egamma = 0.;
  
  // 06/02/2016 Shower Non-linearity Correction by Adesh. 
  
  if ( Eclust <= Ecutoff ) { 
  
    Egamma = Eclust/(A*Eclust + B); // Linear part
  
  }
  
  if ( Eclust > Ecutoff ) { 
  
    Egamma = Eclust/(C - exp(-D*Eclust+ E)); // Non-linear part
  
  }
  
  // End Correction  
  

  // then depth corrections 
  if ( Egamma > 0 ) { 
    float dxV = x0-vertex->X();
    float dyV = y0-vertex->Y();
    float zV = vertex->Z();
   
    double z0 = m_FCALfront - zV;
    double zMax = FCAL_RADIATION_LENGTH*(FCAL_SHOWER_OFFSET 
					 + log(Egamma/FCAL_CRITICAL_ENERGY));
    double zed = z0;
    double zed1 = z0 + zMax;

    double r0 = sqrt(dxV*dxV + dyV*dyV );

    int niter;
    for ( niter=0; niter<100; niter++) {
      double tt = r0/zed1;
      zed = z0 + zMax/sqrt( 1 + tt*tt );
      if ( fabs( (zed-zed1) ) < 0.001) {
	break;
      }
      zed1 = zed;
    }
    
    posInCal.SetZ( zed + zV );
    errZ = zed - zed1;
  }
  
  Ecorrected = Egamma;
  pos_corrected = posInCal;
}



jerror_t
DFCALShower_factory::FillCovarianceMatrix(DFCALShower *shower){
  /// This function takes a FCALShower object and using the internal variables
  /// overwrites any existing covaraince matrix using lookup tables.

  // Get edges of lookup table histograms (assume that all histograms have the same limits.)
  TAxis *xaxis = CovarianceLookupTable[0][0]->GetXaxis();
  TAxis *yaxis = CovarianceLookupTable[0][0]->GetYaxis();
  float minElookup = xaxis->GetBinLowEdge(1);
  float maxElookup = xaxis->GetBinUpEdge(xaxis->GetNbins());
  float minthlookup = yaxis->GetBinLowEdge(1);
  float maxthlookup = yaxis->GetBinUpEdge(yaxis->GetNbins());

  float shower_E = shower->getEnergy();
  float shower_x = shower->getPosition().X();
  float shower_y = shower->getPosition().Y();
  float shower_z = shower->getPosition().Z();
  float shower_r = sqrt(shower_x*shower_x + shower_y*shower_y);
  float shower_theta = atan2(shower_r,shower_z);
  float thlookup = shower_theta/3.14159265*180;
  float Elookup = shower_E;

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
      if (i==0 && j==0) val *= shower_E; // E variance is divided by energy in CCDB
      ErphiztCovariance(i,j) = ErphiztCovariance(j,i) = val;
    }
  }

  float shower_phi = atan2(shower_y,shower_x);
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
  for (int i=0; i<5; i++) {
    for (int j=0; j<5; j++)
      shower->ExyztCovariance(i, j) = D(i, j);
  }
  if (VERBOSE>2) {printf("(E,x,y,z,t)    "); shower->ExyztCovariance.Print(); }

  return NOERROR;
}


jerror_t
DFCALShower_factory::LoadCovarianceLookupTables(JEventLoop *eventLoop){
  std::thread::id this_id = std::this_thread::get_id();
  stringstream idstring;
  idstring << this_id;
  if (VERBOSE>0) printf("DFCALShower_factory::LoadCovarianceLookupTables():  Thread %s\n",idstring.str().c_str());

  bool USECCDB=0;
  bool DUMMYTABLES=0;
  // if filename specified try to use filename else get info from CCDB
  if (COVARIANCEFILENAME == "") USECCDB=1;

  map<string,string> covariance_data;
  if (USECCDB) {
    // load information for covariance matrix
    if (eventLoop->GetJCalibration()->GetCalib("/FCAL/shower_covariance", covariance_data)) {
      jerr << "Error loading /FCAL/shower_covariance !" << endl;
      DUMMYTABLES=1;
    }
    if (covariance_data.size() == 15)  {  // there are 15 elements in the covariance matrix
      // for example, print it all out
      if (VERBOSE>0) {
	for(auto element : covariance_data) {
	  cout << "\nTEST:   " << element.first << " = " << element.second << endl;
	}
      }
    } else {
      jerr << "Wrong number of elements /FCAL/shower_covariance !" << endl;
      DUMMYTABLES=1;
    }
  }

  for (int i=0; i<5; i++) {
    for (int j=0; j<=i; j++) {

      japp->RootWriteLock();
      // change directory to memory so that histograms are not saved to file
      TDirectory *savedir = gDirectory;

      char histname[255];
      sprintf(histname,"covariance_%i%i_thread%s",i,j,idstring.str().c_str());
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
	  DUMMYTABLES=1;
	} else {
	  getline(ifs, line, '\n');
	  ss.str(line);
	  if (VERBOSE>1) cout << filename << " dump: " <<line<<endl;
	}
      }
      if (DUMMYTABLES) {
	// create dummy histogram since something went wrong
	CovarianceLookupTable[i][j] = new TH2F(histname,"Covariance histogram",10,0,12,10,0,12);
	CovarianceLookupTable[i][j]->SetDirectory(nullptr);
      } else {
	// Parse string
	int nxbins, nybins;
	ss>>nxbins;
	ss>>nybins;
	if (VERBOSE>1) printf("parsed dump: bins (%i,%i)\n",nxbins,nybins);
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
	CovarianceLookupTable[i][j] = new TH2F(histname,"Covariance histogram",nxbins,xbins,nybins,ybins);
	CovarianceLookupTable[i][j]->SetDirectory(nullptr);
	// fill histogram
	while(ss>>cont){
	  if (VERBOSE>1) printf("(%i,%i) (%i,%i) %e  ",i,j,xbin,ybin,cont);
	  CovarianceLookupTable[i][j]->SetBinContent(xbin,ybin,cont);
	  ybin++;
	  if (ybin>nybins) { xbin++; ybin=1; }
	}
	if (VERBOSE>1) printf("\n");
	// Close file
	ifs.close();
      }
      savedir->cd();
      japp->RootUnLock(); 
    }
  }
  return NOERROR;
}

unsigned int
DFCALShower_factory::getMaxHit( const vector< const DFCALHit* >& hitVec ) const {
  
  unsigned int maxIndex = 0;
  
  double eMaxSh = 0;
  
  for( vector< const DFCALHit* >::const_iterator hit = hitVec.begin();
       hit != hitVec.end(); ++hit ){

    if( (**hit).E > eMaxSh ){

      eMaxSh = (**hit).E;
      maxIndex = hit - hitVec.begin();
    }
  }

  return maxIndex;
}

void
DFCALShower_factory::getUVFromHits( double& sumUSh, double& sumVSh, 
				    const vector< const DFCALHit* >& hits,
				    const DVector3& showerVec,
				    const DVector3& trackVec ) const {

  // This method forms an axis pointing from the shower to nearest track
  // and computes the energy-weighted second moment of the shower along
  // and perpendicular to this axis.  True photons are fairly symmetric
  // and have similar values of sumU and sumV whereas splitoffs tend
  // to be asymmetric in these variables.

  DVector3 u = ( showerVec - trackVec ).Unit();
  DVector3 z( 0, 0, 1 );
  DVector3 v = u.Cross( z );

  DVector3 hitLoc( 0, 0, 0 );

  sumUSh = 0;
  sumVSh = 0;

  double sumE = 0;
  
  for( vector< const DFCALHit* >::const_iterator hit = hits.begin();
       hit != hits.end(); ++hit ){

    hitLoc.SetX( (**hit).x - showerVec.X() );
    hitLoc.SetY( (**hit).y - showerVec.Y() );

    sumUSh += (**hit).E * pow( u.Dot( hitLoc ), 2 );
    sumVSh += (**hit).E * pow( v.Dot( hitLoc ), 2 );

    sumE += (**hit).E;
  }

  sumUSh /= sumE;
  sumVSh /= sumE;
}

void
DFCALShower_factory::getE1925FromHits( double& e1e9Sh, double& e9e25Sh, 
				       const vector< const DFCALHit* >& hits,
				       unsigned int maxIndex ) const {

  double E9 = 0;
  double E25 = 0;

  const DFCALHit* maxHit = hits[maxIndex];
  
  for( vector< const DFCALHit* >::const_iterator hit = hits.begin();
       hit != hits.end(); ++hit ){
     
    if( fabs( (**hit).x - maxHit->x ) < 4.5 && fabs( (**hit).y - maxHit->y ) < 4.5 )
      E9 += (**hit).E;

    if( fabs( (**hit).x - maxHit->x ) < 8.5 && fabs( (**hit).y - maxHit->y ) < 8.5 )
      E25 += (**hit).E;
  }

  e1e9Sh = maxHit->E/E9;
  e9e25Sh = E9/E25;
}


vector< const DTrackWireBased* >
DFCALShower_factory::filterWireBasedTracks( vector< const DTrackWireBased* >& wbTracks ) const {

  vector< const DTrackWireBased* > finalTracks;
  map< unsigned int, vector< const DTrackWireBased* > > sortedTracks;

  // first sort the wire based tracks into lists with a common candidate id
  // this means that they all come from the same track in the detector
  
  for( unsigned int i = 0; i < wbTracks.size(); ++i ){

    unsigned int id = wbTracks[i]->candidateid;

    if( sortedTracks.find( id ) == sortedTracks.end() ){
      
      sortedTracks[id] = vector< const DTrackWireBased* >();
    }

    sortedTracks[id].push_back( wbTracks[i] );
  }

  // now loop through that list of unique tracks and for each set
  // of wire based tracks, choose the one with the highest FOM
  // (this is choosing among different particle hypotheses)
  
  for( map< unsigned int, vector< const DTrackWireBased* > >::const_iterator
	 anId = sortedTracks.begin();
       anId != sortedTracks.end(); ++anId ){

    double maxFOM = 0;
    unsigned int bestIndex = 0;

    for( unsigned int i = 0; i < anId->second.size(); ++i ){

      if( anId->second[i]->Ndof < 15 ) continue;
      
      if( anId->second[i]->FOM > maxFOM ){

	maxFOM = anId->second[i]->FOM;
	bestIndex = i;
      }
    }

    finalTracks.push_back( anId->second[bestIndex] );
  }
  
  return finalTracks;
}
