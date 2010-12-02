// $Id$
//
//    File: DBCALShower_factory.h
// Created: Tue Jul  3 18:25:12 EDT 2007
// Creator: Matthew Shepherd
//
//  Author's Note: This class started as a functional clone of the default DBCALShower_factory --
//  much of the logic and data structure is as written by C. Xu.
//  Revisions were made to get approx 10x speed gain.  The intial version in the repository
//  It has undergone subsequent reision and as of 8-Aug-08 merged back in as the default
//  DBCALShower_factory
//


#ifndef _DBCALShower_factory_
#define _DBCALShower_factory_

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>
using namespace jana;

#include "BCAL/DBCALShower.h"

class DBCALShower_factory:public JFactory<DBCALShower>{
    
public:
    
    DBCALShower_factory();
    ~DBCALShower_factory(){};
    
private:
        jerror_t brun(JEventLoop *loop, int runnumber);
    jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via JEventProcessor virtual method
    void CellRecon(JEventLoop *loop);
    void CeleToArray(void);            
    void PreCluster(JEventLoop *loop);
    void Connect(int,int);
    void ClusNorm(void);
    void ClusAnalysis();
    void Clus_Break(int nclust);
    
    void clsLoop();
    // The following functions are for track fitting (using layer by layer 
    // method)
    
    void Trakfit(void);
    void Fit_ls();
    void Linefit(int ixyz,int mwt,float &a, float &b,float &siga,
                 float &sigb,float &chi2,float &q);
    float Gammq(float a,float x);
    void Gser(float &gamser,float a,float x);
    void Gcf(float &gammcf,float a,float x);
    float Gammln(float xx_gln);

#undef layermax_bcal

#define modulemax_bcal 48
#define layermax_bcal 10
#define colmax_bcal 4
#define cellmax_bcal modulemax_bcal*layermax_bcal*colmax_bcal
#define clsmax_bcal modulemax_bcal*layermax_bcal*colmax_bcal
    
    // the following data member are used bu function CellRecon()
    
    // The quantities like ecel_a,ecel_b,tcel_a,tcel_b,xcel,ycel,zcel,
    // ecel,tcel are quantities in event.inc in Kloe's code which was
    // kept in common blocks for reconstuection function to "mess up".
    // Now in our case we put them as a private members of class
    
    float  ecel_a[modulemax_bcal][layermax_bcal][colmax_bcal]; 
    float  ecel_b[modulemax_bcal][layermax_bcal][colmax_bcal]; 
    float  tcel_a[modulemax_bcal][layermax_bcal][colmax_bcal];   
    float  tcel_b[modulemax_bcal][layermax_bcal][colmax_bcal];
    float  xx[modulemax_bcal][layermax_bcal][colmax_bcal];   
    float  yy[modulemax_bcal][layermax_bcal][colmax_bcal];
    float  fiberLength;
    float  zOffset;
    
    float  xcel[modulemax_bcal][layermax_bcal][colmax_bcal];  
    float  ycel[modulemax_bcal][layermax_bcal][colmax_bcal];   
    float  zcel[modulemax_bcal][layermax_bcal][colmax_bcal];   
    float  tcel[modulemax_bcal][layermax_bcal][colmax_bcal];   
    float  ecel[modulemax_bcal][layermax_bcal][colmax_bcal];
    float  tcell_anor[modulemax_bcal][layermax_bcal][colmax_bcal];   
    float  tcell_bnor[modulemax_bcal][layermax_bcal][colmax_bcal]; 
    float  ta_offset[modulemax_bcal][layermax_bcal][colmax_bcal];
    float  tb_offset[modulemax_bcal][layermax_bcal][colmax_bcal];   
    float  ta0;
    float  tb0;   
    // The above data members are used by function CellRecon()
    
    
    // The following data are used by function CeleToArray();
    float ethr_cell;// MIN ENERGY THRESD OF cel in GeV
    int celtot;  //----- TOTAL NUMBER OF CELLS with readout --
                 //-------------------------------------------------------------------
                 // ARRAY NARR stores the module, layer, and col number of cells
                 //      NARR = Integer data for each cell where :
                 //         1 = Module number
                 //         2 = Layer  number
                 //         3 = Col    number
                 //--------------------------------------------------------------------
    int narr[4][cellmax_bcal+1]; 
    int   nclus[cellmax_bcal+1];  
    int   next[cellmax_bcal+1];
    float e_cel[cellmax_bcal+1];
    float x_cel[cellmax_bcal+1];
    float y_cel[cellmax_bcal+1];
    float z_cel[cellmax_bcal+1];
    float t_cel[cellmax_bcal+1];
    float ta_cel[cellmax_bcal+1];
    float tb_cel[cellmax_bcal+1]; 
    //-------------------------------------------------------------------
    // ARRAY CELDATA stores adc and tdc values of cells
    //      CELDATA = data for each cell where :
    //         1 = ADC Energy in side A
    //         2 = ADC Energy in side B
    //         3 = Time of Arrival in side A
    //         4 = Time of Arrival in side B
    //--------------------------------------------------------------------       
    float  celdata[4][cellmax_bcal+1];
    // The above data are used by function CeleToArray();
    
    // The following data are used by function ClusNorm() and also fitting
    
    int clstot;  //----- TOTAL NUMBER OF CLUSTERS --
    
    float  e_cls[clsmax_bcal+1];  //----- TOTAL ENERGY OF CLUSTER -
    float  x_cls[clsmax_bcal+1];  //----- X CENTROID OF CLUSTER -
    float  y_cls[clsmax_bcal+1];  //----- Y CENTROID OF CLUSTER -
    float  z_cls[clsmax_bcal+1];  //----- Z CENTROID OF CLUSTER -
    float  t_cls[clsmax_bcal+1];  //---- <T> OF CLUSTER -
    float  ea_cls[clsmax_bcal+1]; //---- E SUM OF SIDE A OF CLUSTER -
    float  eb_cls[clsmax_bcal+1]; //---- E SUM OF SIDE B OF CLUSTER - 
    
    float  ta_cls[clsmax_bcal+1]; //--- <T> OF SIDE A OF CLUSTER -
    float  tb_cls[clsmax_bcal+1]; //--- <T> OF SIDE B OF CLUSTER -
    float  tsqr_a[clsmax_bcal+1]; //--- <T^2> OF SIDE A OF CLUSTER - 
    float  tsqr_b[clsmax_bcal+1]; //--- <T^2> OF SIDE B OF CLUSTER -
    float  trms_a[clsmax_bcal+1]; //--- T RMS OF SIDE A OF CLUSTER -
    float  trms_b[clsmax_bcal+1]; //--- T RMS OF SIDE B OF CLUSTER -   
    
    float  e2_a[clsmax_bcal+1];   //--<E^2> OF SIDE A OF CLUSTER -
    float  e2_b[clsmax_bcal+1];   //--<E^2> OF SIDE B OF CLUSTER -
    
    int clspoi[clsmax_bcal+1];  //---POINTER TO FIRST CELL OF CLUSTER CHAIN 
    int ncltot[clsmax_bcal+1];    //---- TOTAL NUMBER OF CELLS INCLUSTER -
    int ntopol[clsmax_bcal+1];
    
    
    // ********************for fitting****************************************
    int  nlrtot[clsmax_bcal+1]; //----total layers of the cluster
    
    int  elyr;
    int  xlyr; 
    int  ylyr; 
    int  zlyr; 
    int  tlyr; 
    
    float clslyr[6][layermax_bcal+1][clsmax_bcal];
    
    float  apx[4][clsmax_bcal+1];
    float  eapx[4][clsmax_bcal+1];
    float  ctrk[4][clsmax_bcal+1];
    float  ectrk[4][clsmax_bcal+1];
    
    float ctrk_ix[4];
    float ectrk_ix[4];
    float apx_ix[4];
    float eapx_ix[4];   
    
    float x[layermax_bcal+1];
    float y[layermax_bcal+1];  
    float z[layermax_bcal+1];
    float e[layermax_bcal+1];
    float rt[layermax_bcal+1];
    
    float sigx[layermax_bcal+1];  // cm
    float sigy[layermax_bcal+1];  // cm
    float sigz[layermax_bcal+1];  // cm
    
    // some thresholds        
    
    float MERGE_THRESH_DIST;    // CENTROID DISTANCE THRESHOLD
    float MERGE_THRESH_TIME;    // CENTROID TIME THRESHOLD
    float MERGE_THRESH_ZDIST;   // FIBER DISTANCE THRESHOLD
    float MERGE_THRESH_XYDIST;  // CENTROID TRANSVERSE DISTANCE THRESHOLD
    float BREAK_THRESH_TRMS;    // T RMS THRESHOLD    
    float CLUST_THRESH;         // MIN ENERGY THRESD OF CLUSTER IN GEV    
    float ATTEN_LENGTH;         // Attenuation lenth and other parameters    
    float C_EFFECTIVE;          // Effective v of light in scintillator    
};

#endif // _DBCALShower_factory_

