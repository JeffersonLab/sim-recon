// $Id$
//
//    File: DDIRCLut.cc
//

#include <cassert>
#include <math.h>
using namespace std;

#include "DDIRCLut.h"
#include "DANA/DApplication.h"
#include <JANA/JCalibration.h>

//---------------------------------
// DDIRCLut    (Constructor)
//---------------------------------
DDIRCLut::DDIRCLut(JEventLoop *loop) 
{
	////////////////////////////////////
	// LUT histograms and diagnostics //
	////////////////////////////////////
	DIRC_DEBUG_HISTS = true;
	//gPARMS->SetDefaultParameter("DIRC:DEBUG_HISTS",DIRC_DEBUG_HISTS);

	/*
	string locDirName = "DIRC";
	TDirectoryFile* locDirectoryFile = static_cast<TDirectoryFile*>(gDirectory->GetDirectory(locDirName.c_str()));
	if(locDirectoryFile == NULL)
		locDirectoryFile = new TDirectoryFile(locDirName.c_str(), locDirName.c_str());
	locDirectoryFile->cd();
	*/	

	// Set expected angle functions: should get parameters from CCDB for sigma_thetaC
	for(int loc_i = 0; loc_i<3; loc_i++) {
		fAngle[loc_i] = new TF1(Form("fAngle_%d",loc_i),"[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0.7,0.9);
		fAngle[loc_i]->SetParameter(0,1);         // const
		fAngle[loc_i]->SetParameter(2,0.0085);    // sigma
		// mean is set later depending on momentum of track
		
		//dapp->Lock(); //japp->RootWriteLock();
		if(DIRC_DEBUG_HISTS) {
			//gDirectory->cd("/");
			cout<<"Initializing angle histograms"<<endl;
			hAngle[loc_i] = (TH1F*)gROOT->FindObject(Form("hAngle_%d",loc_i));
			if(!hAngle[loc_i]) 
				hAngle[loc_i] = new TH1F(Form("hAngle_%d",loc_i),  "cherenkov angle;#theta_{C} [rad];entries/N_{max} [#]", 250,0.6,1);
		}
		//dapp->Unlock(); //japp->RootUnLock();
	}

	if(DIRC_DEBUG_HISTS) {
		//dapp->Lock(); //japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
        	{
			//gDirectory->cd("/");
			cout<<endl<<endl<<"Initialize histograms"<<endl<<endl;
			//hDiff = (TH1F*)gROOT->FindObject("hDiff");
			//hDiffT = (TH1F*)gROOT->FindObject("hDiffT");
			//hDiffD = (TH1F*)gROOT->FindObject("hDiffD");
			//hDiffR = (TH1F*)gROOT->FindObject("hDiffR");
			//hTime = (TH1F*)gROOT->FindObject("hTime");
			//hCalc = (TH1F*)gROOT->FindObject("hCalc");
			//hNph = (TH1F*)gROOT->FindObject("hNph");
			if(!hDiff) hDiff = new TH1F("hDiff",";t_{calc}-t_{measured} [ns];entries [#]", 400,-20,20);
			if(!hDiffT) hDiffT = new TH1F("hDiffT",";t_{calc}-t_{measured} [ns];entries [#]", 400,-20,20);
			if(!hDiffD) hDiffD = new TH1F("hDiffD",";t_{calc}-t_{measured} [ns];entries [#]", 400,-20,20);
			if(!hDiffR) hDiffR = new TH1F("hDiffR",";t_{calc}-t_{measured} [ns];entries [#]", 400,-20,20);
			if(!hTime) hTime = new TH1F("hTime",";propagation time [ns];entries [#]",   1000,0,200);
			if(!hCalc) hCalc = new TH1F("hCalc",";calculated time [ns];entries [#]",   1000,0,200);
			if(!hNph) hNph = new TH1F("hNph",";detected photons [#];entries [#]", 150,0,150);	
		
			cout<<"Found hDiff at "<<gROOT->FindObjectPathName(hDiff)<<endl;
		}
		//dapp->Unlock(); //japp->RootUnLock(); //REMOVE ROOT LOCK!!

	}
	

	DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
	JCalibration *jcalib = dapp->GetJCalibration((loop->GetJEvent()).GetRunNumber());

	/////////////////////////////////
	// retrieve from LUT from file //
	/////////////////////////////////
        const int luts = 48;
	
	TFile *fLut = new TFile("/group/halld/Users/jrsteven/2018-dirc/lut_all_flat.root");
        TTree *tLut=(TTree*) fLut->Get("lut_dirc_flat");
	//tLut->Print();
	
	vector<Double_t> *LutPixelAngleX[luts];
	vector<Double_t> *LutPixelAngleY[luts];
	vector<Double_t> *LutPixelAngleZ[luts];
	vector<Double_t> *LutPixelTime[luts];
	vector<Long64_t> *LutPixelPath[luts];

	// clear arrays to fill from TTree
	for(int l=0; l<luts; l++){
		LutPixelAngleX[l] = 0;
		LutPixelAngleY[l] = 0;
		LutPixelAngleZ[l] = 0;
		LutPixelTime[l] = 0;
		LutPixelPath[l] = 0;
	}

        for(int l=0; l<luts; l++){
		tLut->SetBranchAddress(Form("LUT_AngleX_%d",l),&LutPixelAngleX[l]); 
		tLut->SetBranchAddress(Form("LUT_AngleY_%d",l),&LutPixelAngleY[l]); 
		tLut->SetBranchAddress(Form("LUT_AngleZ_%d",l),&LutPixelAngleZ[l]); 
		tLut->SetBranchAddress(Form("LUT_Time_%d",l),&LutPixelTime[l]); 
		tLut->SetBranchAddress(Form("LUT_Path_%d",l),&LutPixelPath[l]); 
        }

        // fill nodes with LUT info for each bar/pixel combination
	for(int i=0; i<tLut->GetEntries(); i++) { // get pixels from TTree
		tLut->GetEntry(i);
		//if(i<6800 || i>7000) continue;
		//cout<<i<<endl;
		for(int l=0; l<luts; l++){ // loop over bars
			//if(l==0) cout<<"size="<<LutPixelAngleX[l]->size()<<endl;
			for(uint j=0; j<LutPixelAngleX[l]->size(); j++) { // loop over possible paths
				//if(l==31 && i==4432) 
				//	cout<<j<<" "<<LutPixelTime[l]->at(j)<<" "<<LutPixelPath[l]->at(j)<<endl;
				//cout<<"i="<<i<<" l="<<l<<" j="<<j<<endl;
				TVector3 angle(LutPixelAngleX[l]->at(j), LutPixelAngleY[l]->at(j), LutPixelAngleZ[l]->at(j));
				lutNodeAngle[l][i].push_back(angle);
				lutNodeTime[l][i].push_back(LutPixelTime[l]->at(j));
				lutNodePath[l][i].push_back(LutPixelPath[l]->at(j));
			}
		}
	}
	
	//delete LutPixelAngleX;
	//delete LutPixelAngleY;
	//delete LutPixelAngleZ;
	//delete LutPixelTime;
	//delete LutPixelPath;

}

bool DDIRCLut::CalcLUT(TVector3 locProjPos, TVector3 locProjMom, const vector<const DDIRCTruthPmtHit*> locDIRCHits, double locFlightTime, double locMass, shared_ptr<DDIRCMatchParams>& locDIRCMatchParams, shared_ptr<DDIRCLutPhotons>& locDIRCLutPhotons) const
{
	// get bar and track position/momentum from extrapolation
	int bar = 31;
	TVector3 momInBar = locProjMom;
	double tangle,luttheta,evtime;
	int64_t pathid; 
	TVector3 dir,dird;
	double criticalAngle = asin(1.00028/1.47125); // n_quarzt = 1.47125; //(1.47125 <==> 390nm)
	
	double mAngle = acos(sqrt(momInBar.Mag()*momInBar.Mag() + locMass*locMass)/momInBar.Mag()/1.473);
	
	Particle_t hypotheses[3] = {PiPlus, KPlus, Proton};
	for(int loc_i = 0; loc_i<3; loc_i++) {
		double expectedAngle = acos(sqrt(momInBar.Mag()*momInBar.Mag() + ParticleMass(hypotheses[loc_i])*ParticleMass(hypotheses[loc_i]))/momInBar.Mag()/1.473);
		fAngle[loc_i]->SetParameter(1, expectedAngle);
	}

	TVector3 fnX1 = TVector3 (1,0,0);
	TVector3 fnY1 = TVector3 (0,1,0);
	TVector3 fnZ1 = TVector3 (0,0,1);

	// timing cuts for photons
	double cut_tdiffd=2; // direct cut in ns
	double cut_tdiffr=3; // reflected cut in ns

	// loop over DIRC hits
	double logLikelihoodSum[3];
	int nPhotons = 0;
	int nPhotonsThetaC = 0;
	double meanThetaC = 0.;
	for (unsigned int loc_i = 0; loc_i < locDIRCHits.size(); loc_i++){
		const DDIRCTruthPmtHit* locDIRCHit = locDIRCHits[loc_i];
	
		//cout<<"TruthHit i="<<loc_i<<endl;
	
		// cheat and determine bar from truth info
		//bar = locDIRCHit->key_bar;

                // get channel information for LUT
		int pmt = locDIRCHit->ch/64;
		int pix = locDIRCHit->ch%64;
		int sensorId = 100*pmt + pix;

		// use hit time to determine if reflected or not
		double hitTime = locDIRCHit->t;// - locFlightTime;

		// needs to be X dependent choice for reflection cut (from CCDB?)
		bool reflected = hitTime>48;
		
		// get position along bar for calculated time
		double radiatorL = 4*1225; // get from CCDB
		double barend = 2940; // get from CCDB
		double lenz = 0;
		if(locProjPos.Y() < 0) 
			lenz = fabs(barend + locProjPos.X()*10);  
		else 
			lenz = fabs(locProjPos.X()*10 - barend);

		// get length for reflected and direct photons
		double rlenz = 2*radiatorL - lenz; // reflected
		double dlenz = lenz; // direct
		
		if(reflected) lenz = 2*radiatorL - lenz;
		
		bool isGood(false);

		// check for pixel before going through loop
		if(GetLutPixelAngleSize(bar, sensorId) == 0) continue;
		
		// loop over LUT table for this bar/pixel to calculate thetaC	     
		for(uint i = 0; i < GetLutPixelAngleSize(bar, sensorId); i++){
			
			dird   = GetLutPixelAngle(bar, sensorId, i); 
			evtime = GetLutPixelTime(bar, sensorId, i); 
			pathid = GetLutPixelPath(bar, sensorId, i); 
			
			// in MC we can check if the path of the LUT and measured photon are the same
			bool samepath(false);
			if(fabs(pathid - locDIRCHit->path)<0.0001) samepath=true;
			
			for(int r=0; r<2; r++){
				if(!reflected && r==1) continue;
				
				if(r) lenz = rlenz;
				else lenz = dlenz;
				
				for(int u = 0; u < 4; u++){
					if(u == 0) dir = dird;
					if(u == 1) dir.SetXYZ( dird.X(),-dird.Y(),  dird.Z());
					if(u == 2) dir.SetXYZ( dird.X(), dird.Y(), -dird.Z());
					if(u == 3) dir.SetXYZ( dird.X(),-dird.Y(), -dird.Z());
					if(r) dir.SetXYZ( -dir.X(), dir.Y(), dir.Z());
					if(dir.Angle(fnY1) < criticalAngle || dir.Angle(fnZ1) < criticalAngle) continue;
					
					luttheta = dir.Angle(TVector3(-1,0,0));
					if(luttheta > TMath::PiOver2()) luttheta = TMath::Pi()-luttheta;
					tangle = momInBar.Angle(dir);//-0.002; //correction
					
					double bartime = lenz/cos(luttheta)/208.0; 
					double totalTime = bartime+evtime;

					if(DIRC_DEBUG_HISTS) {

						//dapp->RootWriteLock();
						hTime->Fill(hitTime);
						hCalc->Fill(totalTime);
						//cout<<lenz<<endl;
						
						if(fabs(tangle-mAngle)<0.2){
							hDiff->Fill(totalTime-hitTime);
							if(samepath){
								hDiffT->Fill(totalTime-hitTime);
								if(r) hDiffR->Fill(totalTime-hitTime);
								else hDiffD->Fill(totalTime-hitTime);
							}
						}
						//dapp->RootUnLock();
					}

					double locDeltaT = totalTime-hitTime;
					if(!r && fabs(totalTime-hitTime)>cut_tdiffd) continue;
					if(r && fabs(totalTime-hitTime) >cut_tdiffr) continue;

					if(DIRC_DEBUG_HISTS) {
						//dapp->RootWriteLock();
						hAngle[0]->Fill(tangle);
						//cout<<"fill histograms"<<endl;
						//dapp->RootUnLock();
					}
					
					
					// remove photon candidates not used in likelihood
					if(fabs(tangle-mAngle)>0.02) continue;
					
					// save good photons to DIRCLutPhotons object
					isGood = true;
					pair<double, double> photonInfo(tangle, fabs(locDeltaT));
					//locDIRCLutPhotons->dPhoton.push_back(photonInfo);
					
					// calculate average thetaC
					nPhotonsThetaC++;
					meanThetaC += tangle;
					
					// calculate likelihood
					for(int loc_j = 0; loc_j<3; loc_j++) 
						logLikelihoodSum[loc_j] += TMath::Log(fAngle[loc_j]->Eval(tangle)+0.0001);
				}
			} // end loop over reflections

		} // end loop over nodes

		// count good photons
		if(isGood) nPhotons++;

	} // end loop over hits

	// skip tracks without enough photons
	if(nPhotons<5) 
		return false;

	// set DIRCMatchParameters contents
	/*
	
	locDIRCMatchParams->dThetaC = meanThetaC/(double)nPhotonsThetaC;
	locDIRCMatchParams->dLikelihoodPion = logLikelihoodSum[0];
	locDIRCMatchParams->dLikelihoodKaon = logLikelihoodSum[1];
	locDIRCMatchParams->dLikelihoodProton = logLikelihoodSum[2];
	locDIRCMatchParams->dNPhotons = nPhotons;
	*/

	return true;
}

uint DDIRCLut::GetLutPixelAngleSize(int bar, int pixel) const
{
	return lutNodeAngle[bar][pixel].size();
}
	
uint DDIRCLut::GetLutPixelTimeSize(int bar, int pixel) const
{
	return lutNodeTime[bar][pixel].size();
}
	
uint DDIRCLut::GetLutPixelPathSize(int bar, int pixel) const
{
	return lutNodePath[bar][pixel].size();
}

TVector3 DDIRCLut::GetLutPixelAngle(int bar, int pixel, int entry) const
{
	return lutNodeAngle[bar][pixel].at(entry);
}

Double_t DDIRCLut::GetLutPixelTime(int bar, int pixel, int entry) const
{
	return lutNodeTime[bar][pixel].at(entry);
}

Long64_t DDIRCLut::GetLutPixelPath(int bar, int pixel, int entry) const
{
	return lutNodePath[bar][pixel].at(entry);
}
