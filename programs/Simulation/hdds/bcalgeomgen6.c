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
  string mygetstring(int i);
  //
  // Variables from command line:
  //
  int phi_segments = 108; // number of segments in phi
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

  ofstream dst(dest1); 
  if (!dst) 
    { 
      cout << "Could not open " << dest1 << " for writing !!" << endl; 
      return EXIT_FAILURE; 
    } 

  for (int next=1; next<argc; next=next+2)
    {
      //      cout << "Option is ..." << argv[next] << "..." << endl; 
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
	  cout << "Bad Option !!" << endl;
	  badopt = 1;
	}
    }

  if (badopt)
    {
      cout << "Usage: bcalgeomgen [-f ##] [-i ##] [-o ##] [-abcde ##] [-s ##] [-m ##]" << endl;
      return EXIT_FAILURE;
    }

  dst << "<?xml version='1.0' encoding='UTF-8'?>" << endl;
  dst << "<!--DOCTYPE HDDS SYSTEM 'HDDS_1.0.dtd'>" << endl;
  dst << " " << endl;
  dst << "  Hall D Geometry Data Base: Barrel EM calorimeter" << endl;
  dst << "  *************************************************" << endl;
  dst << " " << endl;
  dst << "     version 1.0: Initial version  -rtj" << endl;
  dst << "     version 1.1: Updated version  -ejb" << endl;
  dst << " " << endl;
  dst << "<HDDS specification='v1.0'" << endl;
  dst << "    xmlns='http://www.gluex.org/hdds'" << endl;
  dst << "    xmlns:xsi='http://www.w3.org/2001/XMLSchema-instance'" << endl;
  dst << "    xsi:schemaLocation='http://www.gluex.org/hdds HDDS-1_0.xsd'>" << endl;
  dst << "-->" << endl;
  dst << " " << endl;
  dst << "<section name        = 'BarrelEMcal'" << endl;
  dst << "         version     = '1.1'" << endl;
  dst << "         date        = '2003-10-28'" << endl;
  dst << "         author      = 'E.J. Brash'" << endl;
  dst << "         top_volume  = 'BCAL'" << endl;
  dst << "         specification = 'v1.0'>" << endl;
  dst << " " << endl;
  dst << "<!-- Origin of BarrelEMcal is center of upstream end. -->" << endl;
  dst << " " << endl;
  dst << "  <composition name='BarrelEMcal'>" << endl;
  dst << "    <posXYZ volume='barrelEMcal' X_Y_Z='0.0 0.0 195.0' />" << endl;
  dst << "  </composition>" << endl;
  dst << " " << endl;
  dst << "  <composition name='barrelEMcal' envelope='BCAL'>" << endl;
  dst << "    <mposPhi volume='BCAM' ncopy='" << phi_segments << "' Phi0='0.0' impliedRot='true'>" << endl;
  dst << "      <sector value='1' step='1' />" << endl;
  dst << "    </mposPhi>" << endl;
  dst << "  </composition>" << endl;
  dst << " " << endl;
  dst << "  <composition name='BCAM' envelope='BCLT'>" << endl;
  dst << "    <posXYZ volume='BCMI' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << endl;
  dst << "    <posXYZ volume='BCMJ' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << endl;
  dst << "    <posXYZ volume='BCMK' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << endl;
  dst << "    <posXYZ volume='BCML' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << endl;
  dst << "    <posXYZ volume='BCMM' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << endl;
  dst << "  </composition>" << endl;
  dst << " " << endl;


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

  //cout << nla << " " << nlb << " " << nlc << " " << nld << " " << nle << endl;

  dst << "  <composition name='BCMI' envelope='BCLI'>" << endl;
  for (int i=1; i<=nla; i++) // loop over layers of inner module
    {
      string lstring("BQ");
      string sstring("BV");      
      string endstring = mygetstring(i);
      lstring=lstring+endstring;
      sstring=sstring+endstring;
      dst << "    <posXYZ volume='" << lstring << "' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << endl;
      dst << "    <posXYZ volume='" << sstring << "' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << endl;
    }
  dst << "  </composition>" << endl;
  dst << "  <composition name='BCMJ' envelope='BCLJ'>" << endl;
  for (int i=1; i<=nlb; i++) // loop over layers of next module
    {
      string lstring("BR");
      string sstring("BW");      
      string endstring = mygetstring(i);
      lstring=lstring+endstring;
      sstring=sstring+endstring;
      dst << "    <posXYZ volume='" << lstring << "' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << endl;
      dst << "    <posXYZ volume='" << sstring << "' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << endl;
    }
  dst << "  </composition>" << endl;
  dst << "  <composition name='BCMK' envelope='BCLK'>" << endl;
  for (int i=1; i<=nlc; i++) // loop over layers of next module
    {
      string lstring("BS");
      string sstring("BX");      
      string endstring = mygetstring(i);
      lstring=lstring+endstring;
      sstring=sstring+endstring;
      dst << "    <posXYZ volume='" << lstring << "' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << endl;
      dst << "    <posXYZ volume='" << sstring << "' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << endl;
    }
  dst << "  </composition>" << endl;
  dst << "  <composition name='BCML' envelope='BCLL'>" << endl;
  for (int i=1; i<=nld; i++) // loop over layers of next module
    {
      string lstring("BT");
      string sstring("BY");      
      string endstring = mygetstring(i);
      lstring=lstring+endstring;
      sstring=sstring+endstring;
      dst << "    <posXYZ volume='" << lstring << "' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << endl;
      dst << "    <posXYZ volume='" << sstring << "' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << endl;
    }
  dst << "  </composition>" << endl;
  dst << "  <composition name='BCMM' envelope='BCLM'>" << endl;
  for (int i=1; i<=nle; i++) // loop over layers of outer module
    {
      string lstring("BU");
      string sstring("BZ");      
      string endstring = mygetstring(i);
      lstring=lstring+endstring;
      sstring=sstring+endstring;
      dst << "    <posXYZ volume='" << lstring << "' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << endl;
      dst << "    <posXYZ volume='" << sstring << "' X_Y_Z='0.0 0.0 0.0'></posXYZ>" << endl;
    }
  dst << "  </composition>" << endl;

  dst << " " << endl;
  dst << "  <tubs name='BCAL' Rio_Z='" << r_inner << "  " << r_e << "  390.0' unit_length='cm' material='Air'" << endl;
  dst << "        parameters='barrelEMcal_pars'   comment='barrel EMcal mother'       />" << endl;
  dst << "  <tubs name='BCLT' Rio_Z='" << r_inner << "  " << r_e << "  390.0' profile='0   " << phi_seg_angle << "'" << endl;
  dst << "                    unit_length='cm' material='Air'" << endl;
  dst << "        parameters='barrelEMcal_pars'   comment='barrel EMcal module mother' />" << endl;
  dst << "  <tubs name='BCLI' Rio_Z='" << r_inner << "  " << r_a << "  390.0' profile='0   " << phi_seg_angle << "'" << endl;
  dst << "                    unit_length='cm' material='Air'" << endl;
  dst << "        parameters='barrelEMcal_pars'   comment='barrel EMcal inner module mother' />" << endl;
  dst << "  <tubs name='BCLJ' Rio_Z='" << r_a << "  " << r_b << "  390.0' profile='0   " << phi_seg_angle << "'" << endl;
  dst << "                    unit_length='cm' material='Air'" << endl;
  dst << "        parameters='barrelEMcal_pars'   comment='barrel EMcal middle module mother' />" << endl;
  dst << "  <tubs name='BCLK' Rio_Z='" << r_b << "  " << r_c << "  390.0' profile='0   " << phi_seg_angle << "'" << endl;
  dst << "                    unit_length='cm' material='Air'" << endl;
  dst << "        parameters='barrelEMcal_pars'   comment='barrel EMcal inner module mother' />" << endl;
  dst << "  <tubs name='BCLL' Rio_Z='" << r_c << "  " << r_d << "  390.0' profile='0   " << phi_seg_angle << "'" << endl;
  dst << "                    unit_length='cm' material='Air'" << endl;
  dst << "        parameters='barrelEMcal_pars'   comment='barrel EMcal middle module mother' />" << endl;
  dst << "  <tubs name='BCLM' Rio_Z='" << r_d << "  " << r_e << "  390.0' profile='0   " << phi_seg_angle << "'" << endl;
  dst << "                    unit_length='cm' material='Air'" << endl;
  dst << "        parameters='barrelEMcal_pars'   comment='barrel EMcal outer module mother' />" << endl;
  dst << " " << endl;

  double r1=r_inner;
  double r2=r_inner+pbth_a;

  for (int i=1; i<=nla; i++) // loop over layers of inner module
    {
      string lstring("BQ");
      string sstring("BV");      
      string endstring = mygetstring(i);
      lstring=lstring+endstring;
      sstring=sstring+endstring;
      dst << "  <tubs name='" << lstring << "' Rio_Z='" << r1 << "  " << r2 << "  390.0' profile='0   " << phi_seg_angle << "'" << endl;
      dst << "              unit_length='cm'          unit_angle='deg'" << endl;
      dst << "              material='Lead'      sensitive='true'" << endl;
      dst << "                                  comment='bcal lead aa' />" << endl;
      //cout << "The new lead string is " << lstring << " " << r1 << " " << r2 << endl;
      r1=r2;
      r2=r2+th_sc;
      dst << "  <tubs name='" << sstring << "' Rio_Z='" << r1 << "  " << r2 << "  390.0' profile='0   " << phi_seg_angle << "'" << endl;
      dst << "              unit_length='cm'          unit_angle='deg'" << endl;
      dst << "              material='Scintillator'      sensitive='true'" << endl;
      dst << "                                  comment='bcal scint aa' />" << endl;
      //cout << "The new scint string is " << sstring << " " << r1 << " " << r2 << endl;
      r1=r2;
      r2=r1+pbth_a;
    }
  r2=r2-pbth_a+pbth_b;
  for (int i=1; i<=nlb; i++) // loop over layers of next module
    {
      string lstring("BR");
      string sstring("BW");      
      string endstring = mygetstring(i);
      lstring=lstring+endstring;
      sstring=sstring+endstring;
      dst << "  <tubs name='" << lstring << "' Rio_Z='" << r1 << "  " << r2 << "  390.0' profile='0   " << phi_seg_angle << "'" << endl;
      dst << "              unit_length='cm'          unit_angle='deg'" << endl;
      dst << "              material='Lead'      sensitive='true'" << endl;
      dst << "                                  comment='bcal lead aa' />" << endl;
      //cout << "The new lead string is " << lstring << " " << r1 << " " << r2 << endl;
      r1=r2;
      r2=r2+th_sc;
      dst << "  <tubs name='" << sstring << "' Rio_Z='" << r1 << "  " << r2 << "  390.0' profile='0   " << phi_seg_angle << "'" << endl;
      dst << "              unit_length='cm'          unit_angle='deg'" << endl;
      dst << "              material='Scintillator'      sensitive='true'" << endl;
      dst << "                                  comment='bcal scint aa' />" << endl;
      //cout << "The new scint string is " << sstring << " " << r1 << " " << r2 << endl;
      r1=r2;
      r2=r1+pbth_b;
    }
  r2=r2-pbth_b+pbth_c;
  for (int i=1; i<=nlc; i++) // loop over layers of next module
    {
      string lstring("BS");
      string sstring("BX");      
      string endstring = mygetstring(i);
      lstring=lstring+endstring;
      sstring=sstring+endstring;
      dst << "  <tubs name='" << lstring << "' Rio_Z='" << r1 << "  " << r2 << "  390.0' profile='0   " << phi_seg_angle << "'" << endl;
      dst << "              unit_length='cm'          unit_angle='deg'" << endl;
      dst << "              material='Lead'      sensitive='true'" << endl;
      dst << "                                  comment='bcal lead aa' />" << endl;
      //cout << "The new lead string is " << lstring << " " << r1 << " " << r2 << endl;
      r1=r2;
      r2=r2+th_sc;
      dst << "  <tubs name='" << sstring << "' Rio_Z='" << r1 << "  " << r2 << "  390.0' profile='0   " << phi_seg_angle << "'" << endl;
      dst << "              unit_length='cm'          unit_angle='deg'" << endl;
      dst << "              material='Scintillator'      sensitive='true'" << endl;
      dst << "                                  comment='bcal scint aa' />" << endl;
      //cout << "The new scint string is " << sstring << " " << r1 << " " << r2 << endl;
      r1=r2;
      r2=r1+pbth_c;
    }
  r2=r2-pbth_c+pbth_d;
  for (int i=1; i<=nld; i++) // loop over layers of next module
    {
      string lstring("BT");
      string sstring("BY");      
      string endstring = mygetstring(i);
      lstring=lstring+endstring;
      sstring=sstring+endstring;
      dst << "  <tubs name='" << lstring << "' Rio_Z='" << r1 << "  " << r2 << "  390.0' profile='0   " << phi_seg_angle << "'" << endl;
      dst << "              unit_length='cm'          unit_angle='deg'" << endl;
      dst << "              material='Lead'      sensitive='true'" << endl;
      dst << "                                  comment='bcal lead aa' />" << endl;
      //cout << "The new lead string is " << lstring << " " << r1 << " " << r2 << endl;
      r1=r2;
      r2=r2+th_sc;
      dst << "  <tubs name='" << sstring << "' Rio_Z='" << r1 << "  " << r2 << "  390.0' profile='0   " << phi_seg_angle << "'" << endl;
      dst << "              unit_length='cm'          unit_angle='deg'" << endl;
      dst << "              material='Scintillator'      sensitive='true'" << endl;
      dst << "                                  comment='bcal scint aa' />" << endl;
      //cout << "The new scint string is " << sstring << " " << r1 << " " << r2 << endl;
      r1=r2;
      r2=r1+pbth_d;
    }
   r2=r2-pbth_d+pbth_e;
   for (int i=1; i<=nle; i++) // loop over layers of outer module
    {
      string lstring("BU");
      string sstring("BZ");      
      string endstring = mygetstring(i);
      lstring=lstring+endstring;
      sstring=sstring+endstring;
      dst << "  <tubs name='" << lstring << "' Rio_Z='" << r1 << "  " << r2 << "  390.0' profile='0   " << phi_seg_angle << "'" << endl;
      dst << "              unit_length='cm'          unit_angle='deg'" << endl;
      dst << "              material='Lead'      sensitive='true'" << endl;
      dst << "                                  comment='bcal lead aa' />" << endl;
      //cout << "The new lead string is " << lstring << " " << r1 << " " << r2 << endl;
      r1=r2;
      r2=r2+th_sc;
      dst << "  <tubs name='" << sstring << "' Rio_Z='" << r1 << "  " << r2 << "  390.0' profile='0   " << phi_seg_angle << "'" << endl;
      dst << "              unit_length='cm'          unit_angle='deg'" << endl;
      dst << "              material='Scintillator'      sensitive='true'" << endl;
      dst << "                                  comment='bcal scint aa' />" << endl;
      //cout << "The new scint string is " << sstring << " " << r1 << " " << r2 << endl;
      r1=r2;
      r2=r1+pbth_e;
    }

   dst << " " << endl;
   dst << "  <parameters name='barrelEMcal_pars' type='mcfast'>" << endl;
   dst << "    <real_array name='rmin'     comment='upstream, downstream end inner radius'" << endl;
   dst << "                                values='   " << r_inner << "      " << r_inner << "'" << endl;
   dst << "                                unit='cm'       />" << endl;
   dst << "    <real_array name='rmax'     comment='upstream, downstream end outer radius'" << endl;
   dst << "                                values='   " << r_e << "      " << r_e << "'" << endl;
   dst << "                                unit='cm'       />" << endl;
   dst << "    <real       name='z0'       value='210.0'   comment='z coord. of midplane'" << endl;
   dst << "                                unit='cm'       />" << endl;
   dst << "    <real       name='zlen'     value='390.0'   comment='length of counters'" << endl;
   dst << "                                unit='cm'       />" << endl;
   dst << "    <reference  name='material' value='leadScint' />" << endl;
   dst << "    <reference  name='active'   value='leadScint' />" << endl;
   dst << "    <int        name='nphi'     value='" << phi_segments << "'      comment='phi segmentation' />" << endl;
   dst << "    <int        name='neta'     value='1'       comment='longitudinal seg' />" << endl;
   dst << "    <int        name='nlayers'  value='1'       comment='radial segmentation' />" << endl;
   dst << "    <real       name='siga_em'  value='0.06'    comment='rootE coefficient'" << endl;
   dst << "                                unit='cm'       />" << endl;
   dst << "    <real       name='sigb_em'  value='0.01'    comment='floor term'" << endl;
   dst << "                                unit='cm'       />" << endl;
   dst << "    <real       name='siga_had' value='0.30'    comment='rootE coefficient'" << endl;
   dst << "                                unit='cm'       />" << endl;
   dst << "    <real       name='sigb_had' value='0.03'    comment='floor term'" << endl;
   dst << "                                unit='cm'       />" << endl;
   dst << "    <real       name='em_had_ratio' value='4.0' comment='response ratio'" << endl;
   dst << "                                unit='cm'       />" << endl;
   dst << "  </parameters>" << endl;
   dst << " " << endl;
   dst << "  <mcfast model='EMCal' template='db/emcal.db' parameters='barrelEMcal_pars'>" << endl;
   dst << "    <string     name='name'     value='BCAL'    />" << endl;
   dst << "    <string     name='shape'    value='TUBE'    />" << endl;
   dst << "    <int        name='type'     value='1'       />" << endl;
   dst << "  </mcfast>" << endl;
   dst << " " << endl;
   dst << "</section>" << endl;
   dst << " " << endl;
   dst << "<!-- </HDDS> -->" << endl;

   return EXIT_SUCCESS;

}

string mygetstring(int i)
{
  int j='A'+i/26;
  int k='A'+i-26*(j-'A')-1;

  string fchar,schar;
  fchar=j;
  schar=k;
  string mchar=fchar+schar;

  //cout << "2nd try ..." << fchar << " ... " << schar << " ... " << mchar << "..." << endl; 

  return mchar;
}
