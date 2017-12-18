#include <stdlib.h>
#include <iostream>
#include <map>

#include <JANA/JApplication.h>
#include <JANA/JEvent.h>
#include "DLumi.h"

//---------------------------------
// DLumi    (Constructor)
//---------------------------------
DLumi::DLumi(JEventLoop *loop)
{
  
  compute_lumi = 1;

  // read constants for Lumi determination from calibdb
  std::vector<std::map<string,double> > result;

  loop->GetCalib("/PHOTON_BEAM/lumi/PS_accept", result);

  if ((int)result.size() != DETECTORS) {
    jerr << "Error in DLumi constructor: "
	 << "failed to read constants for PS/PSC acceptances "
	 << "from calibdb at /PHOTON_BEAM/lumi/PS_accept" << std::endl;
    m_psc_accept[0] = 0.;  m_psc_accept[1] = 0.;  m_psc_accept[2] = 0.;
    m_ps_accept[0]  = 0.;  m_ps_accept[1]  = 0.;  m_ps_accept[2]  = 0.;
  } 
  else {

    m_psc_accept[0] = (result[0])["Norm"];
    m_psc_accept[1] = (result[0])["Emin"];
    m_psc_accept[2] = (result[0])["Emax"];
    
    m_ps_accept[0] = (result[1])["Norm"];
    m_ps_accept[1] = (result[1])["Emin"];
    m_ps_accept[2] = (result[1])["Emax"]; 
  }

  loop->GetCalib("/PHOTON_BEAM/lumi/tagm_tagged", result);

  if ((int)result.size() != TAGM_CH) {
    jerr << "Error in DLumi constructor: "
	 << "failed to read constants number of TAGM hits "
	 << "from calibdb at /PHOTON_BEAM/lumi/tagm_tagged" << std::endl;
    for(int ii = 0; ii < TAGM_CH; ii++)
      tagm_tagged[ii] = 0.;

    compute_lumi = 0;

  } 
  else {
    for (int ii = 0; ii < static_cast<int>(result.size()); ii++) {
      int cnt =  (result[ii])["id"];

      if(   ((ii + 1) ==  cnt)  &&  (cnt > 0) && (cnt <= TAGM_CH))  
	tagm_tagged[ii] = (result[ii])["hit"];
      
      else {
	jerr << "Error in DLumi constructor: "
	     << "Invalid counter in the /PHOTON_BEAM/lumi/tagm_tagged table "
	     << std::endl;
	tagm_tagged[ii] = 0;
	compute_lumi = 0;
      }       

    }
  }


  loop->GetCalib("/PHOTON_BEAM/lumi/tagh_tagged", result);
  
  if ((int)result.size() != TAGH_CH) {
    jerr << "Error in DLumi constructor: "
	 << "failed to read constants number of TAGH hits "
	 << "from calibdb at /PHOTON_BEAM/lumi/tagh_tagged" << std::endl;
    for(int ii = 0; ii < TAGH_CH; ii++)
      tagh_tagged[ii] = 0.;
    
    compute_lumi = 0;

  } 
  else {
    for (int ii = 0; ii < static_cast<int>(result.size()); ii++) {
      int cnt =  (result[ii])["id"];

      if(   ((ii + 1) ==  cnt)  &&  (cnt > 0) && (cnt <= TAGH_CH))  
	tagh_tagged[ii] = (result[ii])["hit"];
      else {
	jerr << "Error in DLumi constructor: "
	     << "Invalid counter in the /PHOTON_BEAM/lumi/tagh_tagged table "
	     << std::endl;
	tagh_tagged[ii] = 0;
	compute_lumi = 0;
      }       

    }
  }

  loop->Get( taghGeomVect );
  if (taghGeomVect.size() < 1){
    jerr << "Error in DLumi constructor: "
	 << "Cannot get TAGH geometry "
	 << endl;
    compute_lumi = 0;

  }

  loop->Get( tagmGeomVect );
  if (tagmGeomVect.size() < 1){
    jerr << "Error in DLumi constructor: "
	 << "Cannot get TAGM geometry "
	 << endl;    
    compute_lumi = 0;
  }

  std::map<string,double> result1;
  loop->GetCalib("/PHOTON_BEAM/endpoint_energy", result1);
  if (result1.find("PHOTON_BEAM_ENDPOINT_ENERGY") == result1.end()) {
    std::cerr << "Error in DLumi constructor: "
	      << "failed to read photon beam endpoint energy "
	      << "from calibdb at /PHOTON_BEAM/endpoint_energy" << std::endl;
    Ebeam = 0;
  }
  else {
    Ebeam = result1["PHOTON_BEAM_ENDPOINT_ENERGY"];
  }
  
  
  if(compute_lumi)
    CalcLumi();  
  else {
    jerr << "Error in DLumi constructor: "
	 << "Cannot compute Luminosity "
	 << std::endl;
  }

}

DLumi::~DLumi() { }


void DLumi::CalcLumi(){

  double Norm = m_psc_accept[0];
  double Emin = m_psc_accept[1];
  double Emax = m_psc_accept[2];

  double Epeak = Emin + Emax;
 
  double accept = 0.;


  //  Microscope region 
  for(int ii = 0; ii <  TAGM_CH; ii++){

    double tagm_emin = tagmGeomVect[0]->getElow(1); 
    double tagm_emax = tagmGeomVect[0]->getEhigh(1);    

    double tagm_en   = (tagm_emin + tagm_emax)/2.;

    accept = 0.;
    
    if( (tagm_en < Epeak) && (tagm_en < Epeak)){
      
      accept = 1. - 2.*Emin / tagm_en;
      
      if(accept < 0.) accept = 0.;
      
    } else if( (tagm_en >= Epeak) && (tagm_en >= Epeak)){
      
      accept = 2.*Emax / tagm_en - 1.;

      if(accept < 0.) accept = 0.;
    } else
      accept = 0.;
    
    tagm_lumi[ii] = tagm_tagged[ii]*Norm*accept;

  }


  //  Hodoscope region 
  for(int ii = 0; ii <  TAGH_CH; ii++){

    double tagh_emin = taghGeomVect[0]->getElow(1); 
    double tagh_emax = taghGeomVect[0]->getEhigh(1);
    
    double tagh_en = (tagh_emin + tagh_emax) / 2.;
    
    accept = 0.;
    
    if( (tagh_en < Epeak) && (tagh_en < Epeak)){
      
      accept = 1. - 2.*Emin / tagh_en;
            
    } else if( (tagh_en >= Epeak) && (tagh_en >= Epeak)){
      
      accept = 2.*Emax / tagh_en - 1.;
      
      if(accept < 0.) accept = 0.;
      
    } else
      accept = 0.;
    
    tagh_lumi[ii] = tagh_tagged[ii]*Norm*accept;
    
  }


}

void DLumi::PrintLumi(){
  
  std::cout << std::endl;
  std::cout << " Lumi for TAGM counters (nb)  " << std::endl;
  std::cout << std::endl;
  
  for(int ii = 0; ii <  TAGM_CH; ii++)
    std::cout << " CH = " << ii + 1 <<  "    Lumi =  " <<  tagm_lumi[ii] << std::endl;
  

  std::cout << std::endl;
  std::cout << " Lumi for TAGH counters (nb)  " << std::endl;
  std::cout << std::endl;
  
  for(int ii = 0; ii <  TAGH_CH; ii++)
    std::cout << " CH = " << ii + 1 <<  "    Lumi =  " <<  tagh_lumi[ii] << std::endl;
}

void DLumi::SaveLumi(){

}

