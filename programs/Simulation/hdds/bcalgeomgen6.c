//
// Program to calculate bcal lead/scint layer parameters, and
// then generate xml code to be used in the geometry specification 
// for Hall D Geant simulation.
//
// The data input consist of:
//
// The number of phi modules (default = 108)
// The bcal inner radius (default = 65 cm)
// The bcal outer radius (default = 90 cm)
// The fibre diameter (default = 10mm)
// The lead sheet thickness (default = 5mm)
// The thickness of the inner module layer (default = 8.5 cm)
// The thickness of the middle module layer (default = 8.5 cm)
// 
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>

const char *dest1 = "./bcalgeomgen.xml";

int main(int argc , char *argv[])
{
  std::string mygetstring(int i);
  //
  // Variables from command line:
  //
  int phi_segments = 96; // number of segments in phi
  double r_inner = 65.0; // bcal inner radius in cm
  double r_outer = 90.0; // bcal outer radius in cm
  double fibre_diameter = 1.0; // nominal scint layer thickness
  double th_module = 5.00; // module layer thickness

  double pbth_a = 0.5; // nominal inner module lead thickness
  double pbth_b = 0.5; // nominal next module lead thickness
  double pbth_c = 0.5; // nominal next module lead thickness
  double pbth_d = 0.5; // nominal next module lead thickness
  double pbth_e = 0.5; // nominal outer module lead thickness

  short badopt = 0;

  std::ofstream dst(dest1); 
  if (!dst) 
    { 
      std::cout << "Could not open " << dest1 << " for writing !!" << std::endl; 
      return EXIT_FAILURE; 
    } 

  for (int next=1; next<argc; next=next+2)
    {
      //      std::cout << "Option is ..." << argv[next] << "..." << std::endl; 
      if (!strcmp(argv[next],"-f"))
	{
	  phi_segments = atoi(argv[next+1]);
	}
      else if (!strcmp(argv[next],"-i"))
	{
	  r_inner = atof(argv[next+1]);
	}
      else if (!strcmp(argv[next],"-o"))
	{
	  r_outer = atof(argv[next+1]);
	}
      else if (!strcmp(argv[next],"-a"))
	{
	  pbth_a = atof(argv[next+1]);
	}
      else if (!strcmp(argv[next],"-b"))
	{
	  pbth_b = atof(argv[next+1]);
	}
      else if (!strcmp(argv[next],"-c"))
	{
	  pbth_c = atof(argv[next+1]);
	}
      else if (!strcmp(argv[next],"-d"))
	{
	  pbth_d = atof(argv[next+1]);
	}
      else if (!strcmp(argv[next],"-e"))
	{
	  pbth_e = atof(argv[next+1]);
	}
      else if (!strcmp(argv[next],"-s"))
	{
	  fibre_diameter = atof(argv[next+1]);
	}
      else if (!strcmp(argv[next],"-m"))
	{
	  th_module = atof(argv[next+1]);
	}
      else
	{
	  std::cout << "Bad Option !!" << std::endl;
	  badopt = 1;
	}
    }

  if (badopt)
    {
      std::cout << "Usage: bcalgeomgen [-f ##] [-i ##] [-o ##] [-abcde ##] [-s ##] [-m ##]" << std::endl;
      return EXIT_FAILURE;
    }

  dst << "<?xml version='1.0' encoding='UTF-8'?>" << std::endl;
  dst << "<!--DOCTYPE HDDS SYSTEM 'HDDS_1.0.dtd'>" << std::endl;
  dst << " " << std::endl;
  dst << "  Hall D Geometry Data Base: Barrel EM calorimeter" << std::endl;
  dst << "  *************************************************" << std::endl;
  dst << " " << std::endl;
  dst << "     version 1.0: Initial version  -rtj" << std::endl;
  dst << "     version 1.1: Updated version  -ejb" << std::endl;
  dst << " " << std::endl;
  dst << "<HDDS specification='v1.0'" << std::endl;
  dst << "    xmlns='http://www.gluex.org/hdds'" << std::endl;
  dst << "    xmlns:xsi='http://www.w3.org/2001/XMLSchema-instance'" << std::endl;
  dst << "    xsi:schemaLocation='http://www.gluex.org/hdds HDDS-1_0.xsd'>" << std::endl;
  dst << "-->" << std::endl;
  dst << " " << std::endl;
  dst << "<section name        = 'BarrelEMcal'" << std::endl;
  dst << "         version     = '1.1'" << std::endl;
  dst << "         date        = '2003-10-28'" << std::endl;
  dst << "         author      = 'E.J. Brash'" << std::endl;
  dst << "         top_volume  = 'BCAL'" << std::endl;
  dst << "         specification = 'v1.0'>" << std::endl;
  dst << " " << std::endl;
  dst << "<!-- Origin of BarrelEMcal is center of upstream end. -->" << std::endl;
  dst << " " << std::endl;
  dst << "  <composition name='BarrelEMcal'>" << std::endl;
  dst << "    <posXYZ volume='barrelEMcal' X_Y_Z='0.0 0.0 220.0' />" << std::endl;
  dst << "  </composition>" << std::endl;
  dst << " " << std::endl;
  dst << "  <composition name='barrelEMcal' envelope='BCAL'>" << std::endl;
  dst << "    <mposPhi volume='BCAM' ncopy='" << phi_segments << "' Phi0='0.0' impliedRot='true'>" << std::endl;
  dst << "      <sector value='1' step='1' />" << std::endl;
  dst << "    </mposPhi>" << std::endl;
  dst << "  </composition>" << std::endl;
  dst << " " << std::endl;
  dst << "  <composition name='BCAM' envelope='BCLT'>" << std::endl;
  dst << "    <posXYZ volume='BCMI' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << std::endl;
  dst << "    <posXYZ volume='BCMJ' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << std::endl;
  dst << "    <posXYZ volume='BCMK' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << std::endl;
  dst << "    <posXYZ volume='BCML' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << std::endl;
  dst << "    <posXYZ volume='BCMM' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << std::endl;
  dst << "  </composition>" << std::endl;
  dst << " " << std::endl;


  double Pi=3.14159265;

  // Calculate some other parameters:

  double phi_seg_angle=360.0/phi_segments;
  double thickness = r_outer - r_inner;
  double th_sc = fibre_diameter*3.0/4.0;
  double th_layera = (pbth_a+th_sc);
  double th_layerb = (pbth_b+th_sc);
  double th_layerc = (pbth_c+th_sc);
  double th_layerd = (pbth_d+th_sc);
  double th_layere = (pbth_e+th_sc);
  int nla=th_module/th_layera;
  double r_a = r_inner + nla*th_layera;
  int nlb=th_module/th_layerb;
  double r_b = r_a + nlb*th_layerb;
  int nlc=th_module/th_layerc;
  double r_c = r_b + nlc*th_layerc;
  int nld=th_module/th_layerd;
  double r_d = r_c + nld*th_layerd;
  double th_outer=r_outer-r_d;
  int nle=th_outer/th_layere;
  double r_e = r_d + nle*th_layere;

  //std::cout << nla << " " << nlb << " " << nlc << " " << nld << " " << nle << std::endl;

  dst << "  <composition name='BCMI' envelope='BCLI'>" << std::endl;
  for (int i=1; i<=nla; i++) // loop over layers of inner module
    {
      std::string lstring("BQ");
      std::string sstring("BV");      
      std::string endstring = mygetstring(i);
      lstring=lstring+endstring;
      sstring=sstring+endstring;
      dst << "    <posXYZ volume='" << lstring << "' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << std::endl;
      dst << "    <posXYZ volume='" << sstring << "' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << std::endl;
    }
  dst << "  </composition>" << std::endl;
  dst << "  <composition name='BCMJ' envelope='BCLJ'>" << std::endl;
  for (int i=1; i<=nlb; i++) // loop over layers of next module
    {
      std::string lstring("BR");
      std::string sstring("BW");      
      std::string endstring = mygetstring(i);
      lstring=lstring+endstring;
      sstring=sstring+endstring;
      dst << "    <posXYZ volume='" << lstring << "' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << std::endl;
      dst << "    <posXYZ volume='" << sstring << "' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << std::endl;
    }
  dst << "  </composition>" << std::endl;
  dst << "  <composition name='BCMK' envelope='BCLK'>" << std::endl;
  for (int i=1; i<=nlc; i++) // loop over layers of next module
    {
      std::string lstring("BS");
      std::string sstring("BX");      
      std::string endstring = mygetstring(i);
      lstring=lstring+endstring;
      sstring=sstring+endstring;
      dst << "    <posXYZ volume='" << lstring << "' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << std::endl;
      dst << "    <posXYZ volume='" << sstring << "' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << std::endl;
    }
  dst << "  </composition>" << std::endl;
  dst << "  <composition name='BCML' envelope='BCLL'>" << std::endl;
  for (int i=1; i<=nld; i++) // loop over layers of next module
    {
      std::string lstring("BT");
      std::string sstring("BY");      
      std::string endstring = mygetstring(i);
      lstring=lstring+endstring;
      sstring=sstring+endstring;
      dst << "    <posXYZ volume='" << lstring << "' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << std::endl;
      dst << "    <posXYZ volume='" << sstring << "' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << std::endl;
    }
  dst << "  </composition>" << std::endl;
  dst << "  <composition name='BCMM' envelope='BCLM'>" << std::endl;
  for (int i=1; i<=nle; i++) // loop over layers of outer module
    {
      std::string lstring("BU");
      std::string sstring("BZ");      
      std::string endstring = mygetstring(i);
      lstring=lstring+endstring;
      sstring=sstring+endstring;
      dst << "    <posXYZ volume='" << lstring << "' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << std::endl;
      dst << "    <posXYZ volume='" << sstring << "' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << std::endl;
    }
  dst << "  </composition>" << std::endl;

  dst << " " << std::endl;
  dst << "  <tubs name='BCAL' Rio_Z='" << r_inner << "  " << r_e << "  440.0' unit_length='cm' material='Air'" << std::endl;
  dst << "        parameters='barrelEMcal_pars'   comment='barrel EMcal mother'       />" << std::endl;
  dst << "  <tubs name='BCLT' Rio_Z='" << r_inner << "  " << r_e << "  440.0' profile='0   " << phi_seg_angle << "'" << std::endl;
  dst << "                    unit_length='cm' material='Air'" << std::endl;
  dst << "        parameters='barrelEMcal_pars'   comment='barrel EMcal module mother' />" << std::endl;
  dst << "  <tubs name='BCLI' Rio_Z='" << r_inner << "  " << r_a << "  440.0' profile='0   " << phi_seg_angle << "'" << std::endl;
  dst << "                    unit_length='cm' material='Air'" << std::endl;
  dst << "        parameters='barrelEMcal_pars'   comment='barrel EMcal inner module mother' />" << std::endl;
  dst << "  <tubs name='BCLJ' Rio_Z='" << r_a << "  " << r_b << "  440.0' profile='0   " << phi_seg_angle << "'" << std::endl;
  dst << "                    unit_length='cm' material='Air'" << std::endl;
  dst << "        parameters='barrelEMcal_pars'   comment='barrel EMcal middle module mother' />" << std::endl;
  dst << "  <tubs name='BCLK' Rio_Z='" << r_b << "  " << r_c << "  440.0' profile='0   " << phi_seg_angle << "'" << std::endl;
  dst << "                    unit_length='cm' material='Air'" << std::endl;
  dst << "        parameters='barrelEMcal_pars'   comment='barrel EMcal inner module mother' />" << std::endl;
  dst << "  <tubs name='BCLL' Rio_Z='" << r_c << "  " << r_d << "  440.0' profile='0   " << phi_seg_angle << "'" << std::endl;
  dst << "                    unit_length='cm' material='Air'" << std::endl;
  dst << "        parameters='barrelEMcal_pars'   comment='barrel EMcal middle module mother' />" << std::endl;
  dst << "  <tubs name='BCLM' Rio_Z='" << r_d << "  " << r_e << "  440.0' profile='0   " << phi_seg_angle << "'" << std::endl;
  dst << "                    unit_length='cm' material='Air'" << std::endl;
  dst << "        parameters='barrelEMcal_pars'   comment='barrel EMcal outer module mother' />" << std::endl;
  dst << " " << std::endl;

  double r1=r_inner;
  double r2=r_inner+pbth_a;

  for (int i=1; i<=nla; i++) // loop over layers of inner module
    {
      std::string lstring("BQ");
      std::string sstring("BV");      
      std::string endstring = mygetstring(i);
      lstring=lstring+endstring;
      sstring=sstring+endstring;
      dst << "  <tubs name='" << lstring << "' Rio_Z='" << r1 << "  " << r2 << "  440.0' profile='0   " << phi_seg_angle << "'" << std::endl;
      dst << "              unit_length='cm'          unit_angle='deg'" << std::endl;
      dst << "              material='Lead'      sensitive='true'" << std::endl;
      dst << "                                  comment='bcal lead aa' />" << std::endl;
      //std::cout << "The new lead string is " << lstring << " " << r1 << " " << r2 << std::endl;
      r1=r2;
      r2=r2+th_sc;
      dst << "  <tubs name='" << sstring << "' Rio_Z='" << r1 << "  " << r2 << "  440.0' profile='0   " << phi_seg_angle << "'" << std::endl;
      dst << "              unit_length='cm'          unit_angle='deg'" << std::endl;
      dst << "              material='Scintillator'      sensitive='true'" << std::endl;
      dst << "                                  comment='bcal scint aa' />" << std::endl;
      //std::cout << "The new scint string is " << sstring << " " << r1 << " " << r2 << std::endl;
      r1=r2;
      r2=r1+pbth_a;
    }
  r2=r2-pbth_a+pbth_b;
  for (int i=1; i<=nlb; i++) // loop over layers of next module
    {
      std::string lstring("BR");
      std::string sstring("BW");      
      std::string endstring = mygetstring(i);
      lstring=lstring+endstring;
      sstring=sstring+endstring;
      dst << "  <tubs name='" << lstring << "' Rio_Z='" << r1 << "  " << r2 << "  440.0' profile='0   " << phi_seg_angle << "'" << std::endl;
      dst << "              unit_length='cm'          unit_angle='deg'" << std::endl;
      dst << "              material='Lead'      sensitive='true'" << std::endl;
      dst << "                                  comment='bcal lead aa' />" << std::endl;
      //std::cout << "The new lead string is " << lstring << " " << r1 << " " << r2 << std::endl;
      r1=r2;
      r2=r2+th_sc;
      dst << "  <tubs name='" << sstring << "' Rio_Z='" << r1 << "  " << r2 << "  440.0' profile='0   " << phi_seg_angle << "'" << std::endl;
      dst << "              unit_length='cm'          unit_angle='deg'" << std::endl;
      dst << "              material='Scintillator'      sensitive='true'" << std::endl;
      dst << "                                  comment='bcal scint aa' />" << std::endl;
      //std::cout << "The new scint string is " << sstring << " " << r1 << " " << r2 << std::endl;
      r1=r2;
      r2=r1+pbth_b;
    }
  r2=r2-pbth_b+pbth_c;
  for (int i=1; i<=nlc; i++) // loop over layers of next module
    {
      std::string lstring("BS");
      std::string sstring("BX");      
      std::string endstring = mygetstring(i);
      lstring=lstring+endstring;
      sstring=sstring+endstring;
      dst << "  <tubs name='" << lstring << "' Rio_Z='" << r1 << "  " << r2 << "  440.0' profile='0   " << phi_seg_angle << "'" << std::endl;
      dst << "              unit_length='cm'          unit_angle='deg'" << std::endl;
      dst << "              material='Lead'      sensitive='true'" << std::endl;
      dst << "                                  comment='bcal lead aa' />" << std::endl;
      //std::cout << "The new lead string is " << lstring << " " << r1 << " " << r2 << std::endl;
      r1=r2;
      r2=r2+th_sc;
      dst << "  <tubs name='" << sstring << "' Rio_Z='" << r1 << "  " << r2 << "  440.0' profile='0   " << phi_seg_angle << "'" << std::endl;
      dst << "              unit_length='cm'          unit_angle='deg'" << std::endl;
      dst << "              material='Scintillator'      sensitive='true'" << std::endl;
      dst << "                                  comment='bcal scint aa' />" << std::endl;
      //std::cout << "The new scint string is " << sstring << " " << r1 << " " << r2 << std::endl;
      r1=r2;
      r2=r1+pbth_c;
    }
  r2=r2-pbth_c+pbth_d;
  for (int i=1; i<=nld; i++) // loop over layers of next module
    {
      std::string lstring("BT");
      std::string sstring("BY");      
      std::string endstring = mygetstring(i);
      lstring=lstring+endstring;
      sstring=sstring+endstring;
      dst << "  <tubs name='" << lstring << "' Rio_Z='" << r1 << "  " << r2 << "  440.0' profile='0   " << phi_seg_angle << "'" << std::endl;
      dst << "              unit_length='cm'          unit_angle='deg'" << std::endl;
      dst << "              material='Lead'      sensitive='true'" << std::endl;
      dst << "                                  comment='bcal lead aa' />" << std::endl;
      //std::cout << "The new lead string is " << lstring << " " << r1 << " " << r2 << std::endl;
      r1=r2;
      r2=r2+th_sc;
      dst << "  <tubs name='" << sstring << "' Rio_Z='" << r1 << "  " << r2 << "  440.0' profile='0   " << phi_seg_angle << "'" << std::endl;
      dst << "              unit_length='cm'          unit_angle='deg'" << std::endl;
      dst << "              material='Scintillator'      sensitive='true'" << std::endl;
      dst << "                                  comment='bcal scint aa' />" << std::endl;
      //std::cout << "The new scint string is " << sstring << " " << r1 << " " << r2 << std::endl;
      r1=r2;
      r2=r1+pbth_d;
    }
   r2=r2-pbth_d+pbth_e;
   for (int i=1; i<=nle; i++) // loop over layers of outer module
    {
      std::string lstring("BU");
      std::string sstring("BZ");      
      std::string endstring = mygetstring(i);
      lstring=lstring+endstring;
      sstring=sstring+endstring;
      dst << "  <tubs name='" << lstring << "' Rio_Z='" << r1 << "  " << r2 << "  440.0' profile='0   " << phi_seg_angle << "'" << std::endl;
      dst << "              unit_length='cm'          unit_angle='deg'" << std::endl;
      dst << "              material='Lead'      sensitive='true'" << std::endl;
      dst << "                                  comment='bcal lead aa' />" << std::endl;
      //std::cout << "The new lead string is " << lstring << " " << r1 << " " << r2 << std::endl;
      r1=r2;
      r2=r2+th_sc;
      dst << "  <tubs name='" << sstring << "' Rio_Z='" << r1 << "  " << r2 << "  440.0' profile='0   " << phi_seg_angle << "'" << std::endl;
      dst << "              unit_length='cm'          unit_angle='deg'" << std::endl;
      dst << "              material='Scintillator'      sensitive='true'" << std::endl;
      dst << "                                  comment='bcal scint aa' />" << std::endl;
      //std::cout << "The new scint string is " << sstring << " " << r1 << " " << r2 << std::endl;
      r1=r2;
      r2=r1+pbth_e;
    }

   dst << " " << std::endl;
   dst << "  <parameters name='barrelEMcal_pars' type='mcfast'>" << std::endl;
   dst << "    <real_array name='rmin'     comment='upstream, downstream end inner radius'" << std::endl;
   dst << "                                values='   " << r_inner << "      " << r_inner << "'" << std::endl;
   dst << "                                unit='cm'       />" << std::endl;
   dst << "    <real_array name='rmax'     comment='upstream, downstream end outer radius'" << std::endl;
   dst << "                                values='   " << r_e << "      " << r_e << "'" << std::endl;
   dst << "                                unit='cm'       />" << std::endl;
   dst << "    <real       name='z0'       value='237.0'   comment='z coord. of midplane'" << std::endl;
   dst << "                                unit='cm'       />" << std::endl;
   dst << "    <real       name='zlen'     value='440.0'   comment='length of counters'" << std::endl;
   dst << "                                unit='cm'       />" << std::endl;
   dst << "    <reference  name='material' value='leadScint' />" << std::endl;
   dst << "    <reference  name='active'   value='leadScint' />" << std::endl;
   dst << "    <int        name='nphi'     value='" << phi_segments << "'      comment='phi segmentation' />" << std::endl;
   dst << "    <int        name='neta'     value='1'       comment='longitudinal seg' />" << std::endl;
   dst << "    <int        name='nlayers'  value='1'       comment='radial segmentation' />" << std::endl;
   dst << "    <real       name='siga_em'  value='0.06'    comment='rootE coefficient'" << std::endl;
   dst << "                                unit='cm'       />" << std::endl;
   dst << "    <real       name='sigb_em'  value='0.01'    comment='floor term'" << std::endl;
   dst << "                                unit='cm'       />" << std::endl;
   dst << "    <real       name='siga_had' value='0.30'    comment='rootE coefficient'" << std::endl;
   dst << "                                unit='cm'       />" << std::endl;
   dst << "    <real       name='sigb_had' value='0.03'    comment='floor term'" << std::endl;
   dst << "                                unit='cm'       />" << std::endl;
   dst << "    <real       name='em_had_ratio' value='4.0' comment='response ratio'" << std::endl;
   dst << "                                unit='cm'       />" << std::endl;
   dst << "  </parameters>" << std::endl;
   dst << " " << std::endl;
   dst << "  <mcfast model='EMCal' template='db/emcal.db' parameters='barrelEMcal_pars'>" << std::endl;
   dst << "    <string     name='name'     value='BCAL'    />" << std::endl;
   dst << "    <string     name='shape'    value='TUBE'    />" << std::endl;
   dst << "    <int        name='type'     value='1'       />" << std::endl;
   dst << "  </mcfast>" << std::endl;
   dst << " " << std::endl;
   dst << "</section>" << std::endl;
   dst << " " << std::endl;
   dst << "<!-- </HDDS> -->" << std::endl;

   return EXIT_SUCCESS;

}

std::string mygetstring(int i)
{
  int j='A'+i/26;
  int k='A'+i-26*(j-'A')-1;

  std::string fchar,schar;
  fchar=j;
  schar=k;
  std::string mchar=fchar+schar;

  //std::cout << "2nd try ..." << fchar << " ... " << schar << " ... " << mchar << "..." << std::endl; 

  return mchar;
}
