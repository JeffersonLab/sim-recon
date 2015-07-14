// $Id$
//
// File: DEventProcessor_bcalfcaltof_res_tree.cc
// Created: Thu Aug 25 11:38:03 EDT 2011
// Creator: pmatt (on Darwin swire-b241.jlab.org 8.4.0 powerpc)
//

#include "DEventProcessor_bcalfcaltof_res_tree.h"

// The executable should define the ROOTfile global variable. It will
// be automatically linked when dlopen is called.
extern TFile *ROOTfile;

// Routine used to create our DEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_bcalfcaltof_res_tree());
}
} // "C"


//------------------
// init
//------------------
jerror_t DEventProcessor_bcalfcaltof_res_tree::init(void)
{

	dBCALStudyFlag = true;
	dFCALStudyFlag = true;
	dTOFStudyFlag = true;

	dBCALMCComparison = new BCALMCComparison();
	dPluginTree_BCALMCComparison = new TTree("dPluginTree_BCALMCComparison", "BCAL MC Comparison");
	dPluginTree_BCALMCComparison->Branch("dPluginBranch_BCALMCComparison", "BCALMCComparison", &dBCALMCComparison);

	dFCALMCComparison = new FCALMCComparison();
	dPluginTree_FCALMCComparison = new TTree("dPluginTree_FCALMCComparison", "FCAL MC Comparison");
	dPluginTree_FCALMCComparison->Branch("dPluginBranch_FCALMCComparison", "FCALMCComparison", &dFCALMCComparison);

	dTOFMCComparison = new TOFMCComparison();
	dPluginTree_TOFMCComparison = new TTree("dPluginTree_TOFMCComparison", "TOF MC Comparison");
	dPluginTree_TOFMCComparison->Branch("dPluginBranch_TOFMCComparison", "TOFMCComparison", &dTOFMCComparison);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_bcalfcaltof_res_tree::brun(JEventLoop *eventLoop, int runnumber)
{
	DApplication* locApplication = dynamic_cast<DApplication*>(eventLoop->GetJApplication());
	if(!locApplication){
		_DBG_<<"Cannot get DApplication from JEventLoop! (are you using a JApplication based program?)"<<endl;
		return RESOURCE_UNAVAILABLE;
	}
	dRootGeom = locApplication->GetRootGeom();


	// Get pointer to TrackFitter object that actually fits a track
	vector<const DTrackFitter *> locTrackFitters;
	eventLoop->Get(locTrackFitters);
	if(locTrackFitters.size()<1){
		_DBG_<<"Unable to get a DTrackFitter object! NO Charged track fitting will be done!"<<endl;
		return RESOURCE_UNAVAILABLE;
	}
  
	dTrackFitter = const_cast<DTrackFitter*>(locTrackFitters[0]);
	// Warn user if something happened that caused us NOT to get a fitter object pointer
	if(!dTrackFitter){
		_DBG_<<"Unable to get a DTrackFitter object! NO Charged track fitting will be done!"<<endl;
		return RESOURCE_UNAVAILABLE;
	}


	double rho_Z_over_A_LnI, radlen;
	dRootGeom->FindMat("Scintillator",dRhoZoverA,rho_Z_over_A_LnI, radlen); //defined both in hddsroot.C and hdgeant.C in the programs/Analysis/hdEventViewer/ folder
	dKRhoZoverA = 0.1535*dRhoZoverA/1000.0;
	dLnI = rho_Z_over_A_LnI/dRhoZoverA;

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_bcalfcaltof_res_tree::evnt(JEventLoop *loop, int eventnumber)
{
	unsigned int loc_i, loc_j;
	DVector3 locLabHitPosition, locLabHitPosition_Truth, locGeneratedVertex;
//	locGeneratedVertex.SetXYZ(0.0, 0.0, 65.0); //should be grabbed (whatev)
	vector<const DMCThrown*> locDMCThrownVector;
	const DMCThrown *locDMCThrown;
	loop->Get(locDMCThrownVector);

	if(dBCALStudyFlag == true){

		vector<const DBCALShower*> locBCALShowerVector;
		const DBCALShower *locBCALShower;
		float locDeltaR, locDeltaZ, locDeltaPhi, locDeltaE, locDeltaT;
		float locTrueR, locTrueZ, locTruePhi, locTrueE, locTrueT;
		float locTrueTheta, locPathLength, locPathLengthCorrection, locSurfaceR, locVelocity;
		float locActualPathLength, locActualTheta, locActualPathLengthCorrection;
		DVector3 locActualPathVector;

		loop->Get(locBCALShowerVector, "KLOE");

//can't use truth hits!:
  //project along track path until a designated distance.
    //the designated distance needs to be such that the MEAN delta-z & delta-r are zero
    //designated distance = dist to intersection with bcal + PathLength_correction, where PathLength_correction = f(theta to beamline (except at edges of chamber!), E) (thrown)

		//if more than 1 bcal hit, clustering didn't work, skip event!!
		if(locBCALShowerVector.size() > 1)
			return NOERROR;

		for(loc_i = 0; loc_i < locBCALShowerVector.size(); loc_i++){
			locBCALShower = locBCALShowerVector[loc_i];
			locLabHitPosition.SetXYZ(locBCALShower->x, locBCALShower->y, locBCALShower->z);

			locSurfaceR = 64.3; //from DBCALGeometry.cc

			// match to generated photon to get energy
			for(loc_j = 0; loc_j < locDMCThrownVector.size(); loc_j++){ //assumes generated at center of target (at least in xy)
				locDMCThrown = locDMCThrownVector[loc_j];
				if(locDMCThrown->type != int(Gamma))
					continue;
				locTrueE = locDMCThrown->energy();
				locGeneratedVertex = locDMCThrown->position();
				//project hit on bcal surface
				locPathLengthCorrection = Calc_BCALPathLengthCorrection(locTrueE);
//				locPathLengthCorrection = 0.0;

				locTrueTheta = locDMCThrown->momentum().Theta();
				locPathLength = locSurfaceR/sin(locTrueTheta) + locPathLengthCorrection;
				locTrueZ = locGeneratedVertex.Z() + locPathLength*cos(locTrueTheta);
				locPathLengthCorrection += Calc_BCALPathLengthCorrectionZPostE(locTrueZ); //must be e-corrected z!
				locPathLength = locSurfaceR/sin(locTrueTheta) + locPathLengthCorrection;
				locTrueZ = locGeneratedVertex.Z() + locPathLength*cos(locTrueTheta);

				locTrueR = locPathLength*sin(locTrueTheta);
				locVelocity = 29.9792458;
				locTrueT = 0.0 + locPathLength/locVelocity;
				locTruePhi = locDMCThrown->momentum().Phi();

				break; //just in case!
			}

			//hist to find PathLength correction!
			locActualPathVector = locLabHitPosition - locGeneratedVertex;
			locActualPathLength = locActualPathVector.Mag();
			locActualTheta = locActualPathVector.Theta();
			locActualPathLengthCorrection = locActualPathLength - locSurfaceR/sin(locTrueTheta);

			//differences
			locDeltaR = locLabHitPosition.Perp() - locTrueR;
			locDeltaPhi = locLabHitPosition.Phi() - locTruePhi;
			locDeltaZ = locLabHitPosition.Z() - locTrueZ;
			locDeltaE = locBCALShower->E - locTrueE;
			locDeltaT = locBCALShower->t - locTrueT;

			LockState();
			dBCALMCComparison->dTrueR = locTrueR;
			dBCALMCComparison->dTruePhi = locTruePhi;
			dBCALMCComparison->dTrueZ = locTrueZ;
			dBCALMCComparison->dTrueE = locTrueE;
			dBCALMCComparison->dTrueT = locTrueT;
			dBCALMCComparison->dDeltaR = locDeltaR;
			dBCALMCComparison->dDeltaPhi = locDeltaPhi;
			dBCALMCComparison->dDeltaZ = locDeltaZ;
			dBCALMCComparison->dDeltaE = locDeltaE;
			dBCALMCComparison->dDeltaT = locDeltaT;
			dBCALMCComparison->dPathLengthCorrection = locActualPathLengthCorrection;

			dBCALMCComparison->dShowerUncertaintyX = locBCALShower->xErr;
			dBCALMCComparison->dShowerUncertaintyY = locBCALShower->yErr;
			dBCALMCComparison->dShowerUncertaintyZ = locBCALShower->zErr;
			dBCALMCComparison->dShowerUncertaintyE = (locBCALShower->E >= 0.0) ? 0.0445*sqrt( locBCALShower->E ) + 0.009*locBCALShower->E : 1e-3; //from old DPhoton_factory::makeBCalPhoton() function
			dBCALMCComparison->dShowerUncertaintyT = locBCALShower->tErr;

			dPluginTree_BCALMCComparison->Fill();

			UnlockState();
		} //end DBCALShower loop

	} //end BCAL


	if(dFCALStudyFlag == true){ //energy comparison only works if one generated photon!!!

		vector<const DFCALShower*> locFCALShowerVector;
		const DFCALShower *locFCALShower;
		float locDeltaX, locDeltaY, locDeltaZ, locDeltaE, locDeltaT;
		float locTrueX, locTrueY, locTrueZ, locTrueE, locTrueT;
		float locPathLength, locPathLengthCorrection, locTrueTheta, locTruePhi, locVelocity;
		float locActualPathLength, locActualTheta, locActualPathLengthCorrection;
		DVector3 locActualPathVector, locShowerHitPositionUncertainty;

		loop->Get(locFCALShowerVector);
		float locSurfaceZ = 625.0; //560 from target center, + 65 target center

		//if more than 1 fcal hit, clustering didn't work, skip event!!
		if(locFCALShowerVector.size() > 1)
			return NOERROR;

		for(loc_i = 0; loc_i < locFCALShowerVector.size(); loc_i++){
			locFCALShower = locFCALShowerVector[loc_i];
			locLabHitPosition = locFCALShower->getPosition();

			// match to generated photon to get energy
			for(loc_j = 0; loc_j < locDMCThrownVector.size(); loc_j++){
				locDMCThrown = locDMCThrownVector[loc_j];
				if(locDMCThrown->type != int(Gamma))
					continue;
				locTrueE = locDMCThrown->energy();
				locTrueTheta = locDMCThrown->momentum().Theta();
				locGeneratedVertex = locDMCThrown->position();
				locVelocity = 29.9792458;
				locTruePhi = locDMCThrown->momentum().Phi();
				locPathLengthCorrection = Calc_FCALPathLengthCorrection(locTrueE);
//				locPathLengthCorrection = 0.0;
				locPathLength = (locSurfaceZ - locGeneratedVertex.Z())/cos(locTrueTheta) + locPathLengthCorrection;
				locTrueX = locGeneratedVertex.X() + locPathLength*sin(locTrueTheta)*cos(locTruePhi);
				locTrueY = locGeneratedVertex.Y() + locPathLength*sin(locTrueTheta)*sin(locTruePhi);
				locTrueZ = locGeneratedVertex.Z() + locPathLength*cos(locTrueTheta);
				locTrueT = 0.0 + locPathLength/locVelocity; //worst for low energy & z
				break; //just in case!
			}

			//hist to find PathLength correction!
			locActualPathVector = locLabHitPosition - locGeneratedVertex;
			locActualPathLength = locActualPathVector.Mag();
			locActualTheta = locActualPathVector.Theta();
			locActualPathLengthCorrection = locActualPathLength - (locSurfaceZ - locGeneratedVertex.Z())/cos(locTrueTheta);

			//differences
			locDeltaX = locLabHitPosition.X() - locTrueX;
			locDeltaY = locLabHitPosition.Y() - locTrueY;
			locDeltaZ = locLabHitPosition.Z() - locTrueZ;
			locDeltaE = locFCALShower->getEnergy() - locTrueE;
			locDeltaT = locFCALShower->getTime() - locTrueT;

			LockState();

			dFCALMCComparison->dTrueX = locTrueX;
			dFCALMCComparison->dTrueY = locTrueY;
			dFCALMCComparison->dTrueZ = locTrueZ;
			dFCALMCComparison->dTrueE = locTrueE;
			dFCALMCComparison->dTrueT = locTrueT;
			dFCALMCComparison->dDeltaX = locDeltaX;
			dFCALMCComparison->dDeltaY = locDeltaY;
			dFCALMCComparison->dDeltaZ = locDeltaZ;
			dFCALMCComparison->dDeltaE = locDeltaE;
			dFCALMCComparison->dDeltaT = locDeltaT;
			dFCALMCComparison->dPathLengthCorrection = locActualPathLengthCorrection;

			locShowerHitPositionUncertainty = locFCALShower->getPositionError();
			dFCALMCComparison->dShowerUncertaintyX = locShowerHitPositionUncertainty.X();
			dFCALMCComparison->dShowerUncertaintyY = locShowerHitPositionUncertainty.Y();
			dFCALMCComparison->dShowerUncertaintyZ = locShowerHitPositionUncertainty.Z();
			dFCALMCComparison->dShowerUncertaintyE = (locFCALShower->getEnergy() >= 0.0) ? 0.042*sqrt(locFCALShower->getEnergy()) + 0.0001 : 1e-3; //from old DPhoton_factory::makeFCalPhoton() function
			dFCALMCComparison->dShowerUncertaintyT = 0.0;

			dPluginTree_FCALMCComparison->Fill();

			UnlockState();
		} //end DFCALShower loop

	} //end FCAL


	if(dTOFStudyFlag == true){
		//v is forward-most TOF plane: vertical

		vector<const DTOFPoint*> locTOFPointVector;
		vector<const DTOFTruth*> locTOFTruthVector;
		vector<const DTOFHit*> locTOFHitVector;

		const DTOFPoint *locTOFPoint;
		const DTOFTruth *locTOFTruth;
		const DTOFHit *locTOFHit;
		float locDeltaX, locDeltaY, locDeltaZ, locDeltadE, locDeltaT;
		float locTheta, locPhi, locVelocity, locDeltaPathLength, locBeta;
		float locTrueX, locTrueY, locTrueZ, locTrueT, locTruedE, locTrueBetaGamma;
		int locParticleType;
		DVector3 locTrueMomentum;
		bool locVerticalPlaneFlag, locHorizontalPlaneFlag;

		loop->Get(locTOFPointVector);
		loop->Get(locTOFTruthVector);

		if((locTOFTruthVector.size() == 1) && (locTOFPointVector.size() == 1)){
			locTOFTruth = locTOFTruthVector[0];
			locLabHitPosition_Truth.SetXYZ(locTOFTruth->x, locTOFTruth->y, locTOFTruth->z);

			locTOFPoint = locTOFPointVector[0];
			locLabHitPosition = locTOFPoint->pos;
			locDeltaZ = locLabHitPosition.Z() - locTOFTruth->z;
			locTrueMomentum.SetXYZ(locTOFTruth->px, locTOFTruth->py, locTOFTruth->pz);
			locBeta = locTrueMomentum.Mag()/locTOFTruth->E;
			locTrueBetaGamma = locBeta/sqrt(1.0 - locBeta*locBeta);
			locDeltaPathLength = locDeltaZ/cos(locTheta);
			locTheta = locTrueMomentum.Theta();
			locPhi = locTrueMomentum.Phi();

			locTrueX = locTOFTruth->x + locDeltaPathLength*sin(locTheta)*cos(locPhi);
			locTrueY = locTOFTruth->y + locDeltaPathLength*sin(locTheta)*sin(locPhi);
			locTrueZ = locLabHitPosition.Z();
			locVelocity = 29.9792458*locTrueMomentum.Mag()/locTOFTruth->E;
			locTrueT = locTOFTruth->t + locDeltaPathLength/locVelocity; //v = x/t, t= x/v
         locParticleType = locTOFTruth->ptype;

			//differences
			locDeltaX = locLabHitPosition.X() - locTrueX;
			locDeltaY = locLabHitPosition.Y() - locTrueY;
			locDeltaZ = locLabHitPosition.Z() - locTrueZ;


			locTruedE = Calc_MostProbableTOFdE(locTrueMomentum, ParticleMass(Particle_t(locParticleType)));
//cout << "recon de, truth de = " << locTOFPoint->dE << ", " << locTruedE << endl;
			locDeltadE = locTOFPoint->dE - locTruedE; //locTOFPoint->dedx is actually dE!!!: rename it!!
			locDeltaT = locTOFPoint->t - locTrueT;

			//hit-based info
			locTOFPoint->GetT(locTOFHitVector);
			locVerticalPlaneFlag = false;
			locHorizontalPlaneFlag = false;
			for(loc_i = 0; loc_i < locTOFHitVector.size(); loc_i++){
				locTOFHit = locTOFHitVector[loc_i];
				if(!((locTOFHit->meantime >= 0.0) || (locTOFHit->meantime <= 0.0)))
					continue;
				if(locTOFHit->orientation)
					locHorizontalPlaneFlag = true;
				else{locVerticalPlaneFlag = true;}
			}

			LockState();

			dTOFMCComparison->dTrueX = locTrueX;
			dTOFMCComparison->dTrueY = locTrueY;
			dTOFMCComparison->dTrueZ = locTrueZ;
			dTOFMCComparison->dTruedE = locTruedE;
			dTOFMCComparison->dTrueT = locTrueT;
			dTOFMCComparison->dDeltaX = locDeltaX;
			dTOFMCComparison->dDeltaY = locDeltaY;
			dTOFMCComparison->dDeltaZ = locDeltaZ;
			dTOFMCComparison->dDeltadE = locDeltadE;
			dTOFMCComparison->dDeltaT = locDeltaT;
			dTOFMCComparison->dVerticalPlaneFlag = locVerticalPlaneFlag;
			dTOFMCComparison->dHorizontalPlaneFlag = locHorizontalPlaneFlag;
			dTOFMCComparison->dTrueBetaGamma = locTrueBetaGamma;
			dPluginTree_TOFMCComparison->Fill();


			UnlockState();
		} //end DTOFPoint loop

	} //end TOF

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_bcalfcaltof_res_tree::erun(void)
{
	// Any final calculations on histograms (like dividing them)
	// should be done here. This may get called more than once.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_bcalfcaltof_res_tree::fini(void)
{
	return NOERROR;
}

void DEventProcessor_bcalfcaltof_res_tree::Convert_Coordinates_BCALToLab(float locBCALR, float locBCALPhi, float locBCALZ, DVector3& locLabVertex){
	//check to make sure phi is in radians and not degrees!!
	locLabVertex.SetXYZ(locBCALR*cos(locBCALPhi), locBCALR*sin(locBCALPhi), locBCALZ);
}

void DEventProcessor_bcalfcaltof_res_tree::Convert_Coordinates_LabToBCAL(const DVector3& locLabVertex, float& locBCALR, float& locBCALPhi, float& locBCALZ){
	//check to make sure phi is in radians and not degrees!!
	locBCALR = locLabVertex.Perp();
	locBCALPhi = locLabVertex.Phi();
	locBCALZ = locLabVertex.Z();
}

double DEventProcessor_bcalfcaltof_res_tree::Calc_MostProbableTOFdE(const DVector3 &locMomentum, double mass){
	float locPaddleThickness = 2.54; //get from database!
   double locPathLength = locPaddleThickness/cos(locMomentum.Theta());
//dtofpoint e is averaged to a single plane, so just find the angle & use 2.54 for the thickness
	double locdE, locdEdx;
	Calc_MostProbableTOFdEdx(locMomentum.Mag(), mass, locPathLength, locdE, locdEdx);
//cout << "mean dedx, p, mass, path = " << locMeanTOFdEdx << ", " << locMomentum.Mag() << ", " << mass << ", " << locPathLength << endl;
	return locdE;
}

void DEventProcessor_bcalfcaltof_res_tree::Calc_MostProbableTOFdEdx(double p, double mass, double dx, double &dE, double &dEdx){

  double betagamma=p/mass;
  double beta2=1./(1.+1./betagamma/betagamma);
  if (beta2<1e-6) beta2=1e-6;
  
  // Electron mass 
  double Me=0.000511; //GeV

  // First (non-logarithmic) term in Bethe-Bloch formula
  double locdETerm = dKRhoZoverA*dx/beta2; //divide by 2???
 
	//density factor
	double locX0 = 0.1464;
	double locX1 = 2.4855;
	double locC = -3.1997;
	double locA = 0.16101;
	double locM = 3.2393;
	double locX = log10(betagamma);
	double locDelta; //density factor!
	if(locX < locX0)
		locDelta = 0.0;
	else if(locX < locX1)
		locDelta = 4.6052*locX + locC + locA*pow(locX1 - locX, locM);
	else
		locDelta = 4.6052*locX + locC;

//	double locDeltaNew = dTrackFitter->CalcDensityEffect(betagamma, dRhoZoverA, dLnI);
//	locDelta = locDeltaNew;
//	cout << "delta, deltanew, diff = " << locDelta << ", " << locDeltaNew << ", " << (locDelta - locDeltaNew) << endl;

  // Most probable energy loss from Landau theory (see Leo, pp. 51-52)
	//more appropriately, from RPP 2011, section 27.2.7
	dE = locdETerm*(log(locdETerm) - log((1.-beta2)/2./Me/beta2) - 2.*dLnI - beta2 + 0.200 - locDelta); //1000.0 converts MeV to GeV
	dEdx = dE/dx;

/*
TGeoMixture *mat25= new TGeoMixture("Scintillator",2,1.032); //2 elements, density = 1.032
mat25->SetUniqueID(25);
mat25->DefineElement(0,12.011,6,0.913734); //element index, A, Z, weight fraction
mat25->DefineElement(1,1.00797,1,0.0862662); //element index, A, Z, weight fraction
*/
}

float DEventProcessor_bcalfcaltof_res_tree::Calc_BCALPathLengthCorrection(float locEnergy){
	if(locEnergy >= 1.0){
		TF1 locFunction("df_BCAL_DepthCorrection", "[0] + [1]*x + [2]*exp([3]*x)", 0.0, 9.0);
		locFunction.SetParameters(9.95659, 0.142382, -2.9869, -0.56881);
		return locFunction.Eval(locEnergy);
	}
	TF1 locFunction("df_BCAL_DepthCorrection_LowE", "[0] + [1]*x + [2]*x*x", 0.0, 1.0);
	locFunction.SetParameters(4.89141, 5.98362, -2.56308);
	return locFunction.Eval(locEnergy);
}

float DEventProcessor_bcalfcaltof_res_tree::Calc_BCALPathLengthCorrectionZPostE(float locZ){
	if(locZ < 22.5){
		TF1 locFunction("df_BCAL_DepthCorrection_LowZ", "[0] + [1]*x + [2]*x*x", 0.0, 1.0);
		locFunction.SetParameters(-31.4058, 2.40814, -0.0453179);
		return locFunction.Eval(locZ);
	}else if(locZ > 389.0){
		TF1 locFunction("df_BCAL_DepthCorrection_HighZ", "[0] + [1]*x + [2]*x*x", 0.0, 1.0);
		locFunction.SetParameters(-1994.66, 10.2362, -0.0131332);
		return locFunction.Eval(locZ);
	}
	return 0.0;
}

float DEventProcessor_bcalfcaltof_res_tree::Calc_FCALPathLengthCorrection(float locEnergy){
	if(locEnergy >= 1.0){
		TF1 locFunction("df_FCAL_DepthCorrection_HighE", "[0] + [1]*x + [2]*exp([3]*x)", 1.0, 9.0);
		locFunction.SetParameters(17.1957, 0.383558, -6.9412, -0.598059);
		return locFunction.Eval(locEnergy);
	}
	TF1 locFunction("df_FCAL_DepthCorrection_LowE", "[0] + [1]*x + [2]*x*x", 0.0, 1.0);
	locFunction.SetParameters(6.83752, 12.0934, -5.24968);
	return locFunction.Eval(locEnergy);
}

