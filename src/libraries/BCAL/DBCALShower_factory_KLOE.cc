// $Id$
//
//    File: DBCALShower_factory_KLOE.cc
// Created: Tue Jul  3 18:25:12 EDT 2007
// Creator: Matthew Shepherd
//

#include <cassert>
#include <math.h>
#include <map>

#include "BCAL/DBCALHit.h"
#include "BCAL/DBCALTDCHit.h"
#include "BCAL/DBCALPoint.h"
#include "BCAL/DBCALGeometry.h"
#include "BCAL/DBCALShower_factory_KLOE.h"

#include "DANA/DApplication.h"

#include "units.h"

using namespace std;

//------------------
// init
//------------------
jerror_t DBCALShower_factory_KLOE::init()
{
    // this should be lower than cut in mcsmear
	ethr_cell=0.0001;     // MIN ENERGY THRESD OF cell in GeV

	CLUST_THRESH = 0.03;  // MIN ENERGY THRESD OF CLUSTER IN GEV (make this match the value used in the other algorithm)
    
	elyr = 1;
	xlyr = 2; 
	ylyr = 3; 
	zlyr = 4; 
	tlyr = 5;
    
  // these four parameters are used in merging clusters
	MERGE_THRESH_DIST   = 40.0;  // CENTROID DISTANCE THRESHOLD
	MERGE_THRESH_TIME   =  2.5;  // CENTROID TIME THRESHOLD
	MERGE_THRESH_ZDIST  = 30.0;  // FIBER DISTANCE THRESHOLD
	MERGE_THRESH_XYDIST = 40.0;  // CENTROID TRANSVERSE DISTANCE THRESHOLD

  // this parameter is used to break clusters based on rms time
  BREAK_THRESH_TRMS= 5.0;   // T RMS THRESHOLD
    
  gPARMS->SetDefaultParameter( "BCALRECON:CLUST_THRESH", CLUST_THRESH );
  gPARMS->SetDefaultParameter( "BCALRECON:MERGE_THRESH_DIST", MERGE_THRESH_DIST );
  gPARMS->SetDefaultParameter( "BCALRECON:MERGE_THRESH_TIME", MERGE_THRESH_TIME );
  gPARMS->SetDefaultParameter( "BCALRECON:MERGE_THRESH_ZDIST", MERGE_THRESH_ZDIST );
  gPARMS->SetDefaultParameter( "BCALRECON:MERGE_THRESH_XYDIST", MERGE_THRESH_XYDIST );
  gPARMS->SetDefaultParameter( "BCALRECON:BREAK_THRESH_TRMS", BREAK_THRESH_TRMS );
  
  if( !DBCALGeometry::summingOn() ){

    // these are energy calibration parameters -- no summing of cells

    m_scaleZ_p0 =  0.9597;
    m_scaleZ_p1 =  0.000454875;
    m_scaleZ_p2 =  -2.29912e-06;
    m_scaleZ_p3 =  1.49757e-09;
  
    m_nonlinZ_p0 =  -0.00154122;
    m_nonlinZ_p1 =  6.73594e-05;
    m_nonlinZ_p2 =  0;
    m_nonlinZ_p3 =  0;
  
  }
  else{
    
    // these are energy calibration parameters -- 1.2.3.4 summing
  
    //last updated for svn revision 9233
    m_scaleZ_p0 = 0.99798;
    m_scaleZ_p1 = 0.000361096;
    m_scaleZ_p2 = -2.17338e-06;
    m_scaleZ_p3 = 1.32201e-09;
  
    m_nonlinZ_p0 = -0.0201272;
    m_nonlinZ_p1 = 0.000103649;
    m_nonlinZ_p2 = 0;
    m_nonlinZ_p3 = 0;
  }
  
  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DBCALShower_factory_KLOE::brun(JEventLoop *loop, int32_t runnumber)
{
    //get target position
    DApplication* app = dynamic_cast<DApplication*>(loop->GetJApplication());
    DGeometry* geom = app->GetDGeometry(runnumber);
    geom->GetTargetZ(m_z_target_center);
    
    vector<const DBCALGeometry*> bcalGeomVect;
    loop->Get( bcalGeomVect );
    const DBCALGeometry& bcalGeom = *(bcalGeomVect[0]);
    
    //////////////////////////////////////////////////////////////////
    // Calculate Cell Position
    //////////////////////////////////////////////////////////////////
    
    ATTEN_LENGTH = bcalGeom.ATTEN_LENGTH;
    C_EFFECTIVE = bcalGeom.C_EFFECTIVE;
    
    fiberLength = bcalGeom.GetBCAL_length(); // fiber length in cm
    zOffset = bcalGeom.GetBCAL_center();

    //the following uses some sad notation in which modmin=0 and modmax=48, when in fact there are 48 modules labelled either 0-47 or 1-48 depending on one's whim, although if we are using the methods from DBCALGeometry (e.g. cellId()), we must start counting from 1 and if we are accessing arrays we must of course start from 0
    int   modmin = 0;
    int   modmax = bcalGeom.NBCALMODS;
    int   rowmin1=0;
    int   rowmax1= bcalGeom.NBCALLAYSIN;
    int   rowmin2= rowmax1;
    int   rowmax2= bcalGeom.NBCALLAYSOUT+rowmin2; 
    int   colmin1=0;
    int   colmax1=bcalGeom.NBCALSECSIN;
    int   colmin2=0;
    int   colmax2=bcalGeom.NBCALSECSOUT;
    
    float r_inner= bcalGeom.GetBCAL_inner_rad();
    
    for (int i = (rowmin1+1); i < (rowmax1+1); i++){
        //this loop starts from 1, so we can use i in cellId with no adjustment
        int cellId = bcalGeom.cellId(1,i,1); //this gives us a cellId that we can use to get the radius. the module and sector numbers are irrelevant
        //rt is radius of center of layer - BCAL inner radius
        rt[i]=bcalGeom.r(cellId)-r_inner;
    }
    
    for (int i = (rowmin2+1); i < (rowmax2+1); i++){
        int cellId = bcalGeom.cellId(1,i,1);
        rt[i]=bcalGeom.r(cellId)-r_inner;
    }
    
    //these are r and phi positions of readout cells
    float r[modulemax_bcal][layermax_bcal][colmax_bcal];
    float phi[modulemax_bcal][layermax_bcal][colmax_bcal];
    
    // Now start to extract cell position information from Geometry class
    for (int k = modmin; k < modmax; k++){
        for (int i = rowmin1; i < rowmax1; i++){
            for (int j = colmin1; j < colmax1; j++){
                //in this case the loops start at 0, so we have to add 1 to the indices when calling cellId(). hooray!
                //use DBCALGeometry to get r/phi position of each cell
                int cellId = bcalGeom.cellId(k+1,i+1,j+1);
                r[k][i][j]=bcalGeom.r(cellId);
                phi[k][i][j]=bcalGeom.phi(cellId);
                //set x and y positions
                xx[k][i][j]=r[k][i][j]*cos(phi[k][i][j]);
                yy[k][i][j]=r[k][i][j]*sin(phi[k][i][j]);          
            }
        }
        
        for (int i = rowmin2; i < rowmax2; i++){
            for (int j = colmin2; j < colmax2; j++){
                int cellId = bcalGeom.cellId(k+1,i+1,j+1);
                r[k][i][j]=bcalGeom.r(cellId);
                phi[k][i][j]=bcalGeom.phi(cellId);

                xx[k][i][j]=r[k][i][j]*cos(phi[k][i][j]);
                yy[k][i][j]=r[k][i][j]*sin(phi[k][i][j]);          
            }
        }         
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // Now the cell information are already contained in xx and yy arrays.
    // xx and yy arrays are private members of this class
    ////////////////////////////////////////////////////////////////////////////
    
	 return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DBCALShower_factory_KLOE::evnt(JEventLoop *loop, uint64_t eventnumber)
{
    // Call core KLOE reconstruction routines
    CellRecon(loop);
    CeleToArray();   
    PreCluster(loop); 
    ClusNorm();
    ClusAnalysis();
    Trakfit();
    
    //Loop over reconstructed clusters and make DBCALShower objects out of them
    vector<DBCALShower*> clusters;    
    int id = 0;
    for (int i = 1; i < (clstot+1); i++){
        
        int  j=clspoi[i];
  
        if( e_cls[j] < CLUST_THRESH ) continue;
        
        // Time to cook a final shower
        DBCALShower *shower = new DBCALShower;
  
        //The algorithm has so far both clustered together hits and determined
        //the position (x,y,z,t,E) of these clusters.
        //For now, we will use only the clustering and recompute the position
        //here. The reason we do this is because the clustering process is
        //unaware of the errors on measurements of z (hits with TDC information
        //will have much better z resolution). By taking into account this
        //information we can greatly improve the z-resolution of clusters.

        /*
        //this is the old way of setting shower properties
        shower->id                  = id++;
        shower->E_raw               = e_cls[j];
        shower->x                   = x_cls[j];
        shower->y                   = y_cls[j];
        shower->z                   = z_cls[j] + zOffset;   
        shower->t                   = t_cls[j];
        shower->N_cell              = ncltot[j];
      
        shower->xErr                = eapx[1][j];
        shower->yErr                = eapx[2][j];
        shower->zErr                = eapx[3][j];

        shower->tErr                = 0.5 * sqrt( trms_a[j] * trms_a[j] +
                                                  trms_b[j] * trms_b[j] );
        */

        // Trace back to the DBCALPoint objects used in this shower and
        // add them as associated objects.
        vector<const DBCALPoint*> pointsInShower;
        FindPointsInShower(j, loop, pointsInShower);

        //Determine cluster (x,y,z,t) by averaging (x,y,z,t) of constituent
        //DBCALPoints.

        //For now just do the most naive averaging (weighting by E) to get
        //the cluster properties (x,y,t). For z, average with weight of
        //1/sig_z^2.
        //Should consider a different weighting scheme (weighting by E^2) or average different quantities (cylindrical or spherical coordinates instead of rectangular)
        double E=0,x=0,y=0,z=0,t=0;
        int N_cell=0;
        double sig_x=0,sig_y=0,sig_z=0,sig_t=0;
        double sum_z_wt=0;
        for(unsigned int j=0; j<pointsInShower.size(); j++){
            double cell_E = pointsInShower[j]->E();
            double cell_r = pointsInShower[j]->r();
            double cell_phi = pointsInShower[j]->phi();
            E += cell_E;
            x += cell_E*cell_r*cos(cell_phi);
            sig_x += cell_E*cell_r*cos(cell_phi)*cell_r*cos(cell_phi);
            y += cell_E*cell_r*sin(cell_phi);
            sig_y += cell_E*cell_r*sin(cell_phi)*cell_r*sin(cell_phi);


            double z_wt = 1/(pointsInShower[j]->sigZ()*pointsInShower[j]->sigZ());
            double cell_z = pointsInShower[j]->z();
            double cell_t = pointsInShower[j]->t();

            sum_z_wt += z_wt;
            z += z_wt*cell_z;
            sig_z += z_wt*cell_z*cell_z;

            t += cell_E*cell_t;
            sig_t += cell_E*cell_t*cell_t;
            N_cell++;

            shower->AddAssociatedObject(pointsInShower[j]);
        }

        x /= E;
        sig_x /= E;
        sig_x = sqrt(sig_x - x*x)/sqrt(N_cell); //really this should be n_effective rather than n, change this later
        y /= E;
        sig_y /= E;
        sig_y = sqrt(sig_y - y*y)/sqrt(N_cell);
        z /= sum_z_wt;
        sig_z /= sum_z_wt;
        sig_z = sqrt(sig_z - z*z)/sqrt(N_cell);
        t /= E;
        sig_t /= E;
        sig_t = sqrt(sig_t - t*t)/sqrt(N_cell);

        shower->id                  = id++;
        shower->E_raw               = E;
        shower->x                   = x;
        shower->y                   = y;
        shower->z                   = z + m_z_target_center;
        shower->t                   = t;
        shower->N_cell              = N_cell;
      
        shower->xErr                = sig_x;
        shower->yErr                = sig_y;
        shower->zErr                = sig_z;

        shower->tErr                = sig_t;
      
        // calibrate energy:
        // Energy calibration has a z dependence -- the
        // calibration comes from fitting E_rec / E_gen to scale * E_gen^nonlin
        // for slices of z.  These fit parameters (scale and nonlin) are then plotted 
        // as a function of z and fit.
  
        float r = sqrt( shower->x * shower->x + shower->y * shower->y );
      
        float zEntry = ( shower->z - m_z_target_center ) * ( DBCALGeometry::GetBCAL_inner_rad() / r );
      
        float scale = m_scaleZ_p0  + m_scaleZ_p1*zEntry + 
            m_scaleZ_p2*(zEntry*zEntry) + m_scaleZ_p3*(zEntry*zEntry*zEntry);
        float nonlin = m_nonlinZ_p0  + m_nonlinZ_p1*zEntry + 
            m_nonlinZ_p2*(zEntry*zEntry) + m_nonlinZ_p3*(zEntry*zEntry*zEntry);

        shower->E = pow( (shower->E_raw ) / scale, 1 / ( 1 + nonlin ) );

        //copy xyz errors into covariance matrix
        shower->xyzCovariance.ResizeTo(3,3);
        shower->xyzCovariance[0][0] = shower->xErr*shower->xErr;
        shower->xyzCovariance[1][1] = shower->yErr*shower->yErr;
        shower->xyzCovariance[2][2] = shower->zErr*shower->zErr;

        _data.push_back(shower);  
    }
    
    return NOERROR;
}


//------------------
// FindPointsInShower()
//------------------
void DBCALShower_factory_KLOE::FindPointsInShower(int indx, JEventLoop *loop, vector<const DBCALPoint*> &pointsInShower)
{
	/// This is called after the clusters have been completely formed. Our
	/// job is simply to find the DBCALPoint objects used to form a given cluster.
	/// This is so the DBCALPoint objects can be added to the DBCALShower object
	/// as AssociatedObjects.

	// The variable indx indexes the next[] array as the starting cell for the
	// cluster. It also indexes the narr[][] array which holds the module, layer,
	// sector(column) values for the hits.
	//
	// Here, we need to loop over next[] elements starting at next[indx] until
	// we find the one pointing to element "indx" (i.e. the start of the list of
	// cells in the cluster.) For each of these, we must find the LAST 
	// member in bcalhits to have the same module, layer, number indicating 
	// that is a DBCALHit to be added.


	vector<const DBCALPoint*> points;
	loop->Get(points);
	
	int start_indx = indx;
	do{
		int module = narr[1][indx];
		int layer  = narr[2][indx];
		int sector = narr[3][indx];
		
		// Loop over BCAL hits, trying to find this one
		for(unsigned int i=0; i<points.size(); i++){
			if(points[i]->module() !=module)continue;
			if(points[i]->layer()  !=layer)continue;
			if(points[i]->sector() !=sector)continue;
			pointsInShower.push_back(points[i]);
		}
		

		indx = next[indx];
	}while(indx != start_indx);

}


//------------------
// CellRecon()
//------------------
void DBCALShower_factory_KLOE::CellRecon(JEventLoop *loop)
{
    //**********************************************************************
    // The main purpose of this function is extracting information
    // from DBCALPoint
    // objects to form several arrays, which will be used later in the function
    // CeleToArray()
    // The four arrays ecel_a,tcel_a,ecel_b,tcel_b hold the energies and times
    // of individual hits (this it the "attenuated" energy as measured at the
    // end of the module, in GeV-equivalent units)  (however times
    // here should be *after* timewalk correction) (a=upstream, b=downstream).
    // The five arrays xcel,ycel,zcel,tcel,ecel hold information about the
    // position, time, and energy of an event in the calorimeter corresponding
    // to a DBCALPoint (one hit at each end of the detector).
    // The final two arrays tcell_anor and tcell_bnor contain the same
    // information as tcel_a and tcell_b.
    // These arrays are 3-D arrays and are rather logically are indexed by
    // [module][layer][column]
    //********************************************************************** 
    
    //First reset the arrays ecel_a,tcel_a,ecel_b,tcel_b to clear out garbage
    //information from the previous events
    memset( ecel_a, 0, modulemax_bcal * layermax_bcal *
            colmax_bcal * sizeof( float ) );
    memset( tcel_a, 0, modulemax_bcal * layermax_bcal *
            colmax_bcal * sizeof( float ) );
    memset( ecel_b, 0, modulemax_bcal * layermax_bcal *
            colmax_bcal * sizeof( float ) );
    memset( tcel_b, 0, modulemax_bcal * layermax_bcal *
            colmax_bcal * sizeof( float ) );
    //the other seven arrays will also be filled with garbage values from
    //previous events, HOWEVER
    //we don't need to zero out these arrays, as long as the ecel_a,tcel_a,ecel_b,tcel_b arrays are zeroed out
    //this is because in CeleToArray(), all cells with ecel_a=0 are skipped
    //and if ecel_a has been to set to a nonzero value for a particular cell
    //then the other arrays will also have been set properly and not full of garbage
    
    vector<const DBCALPoint*> points;
    loop->Get(points);
    if(points.size() <=0) return;

    for (vector<const DBCALPoint*>::const_iterator point_iter = points.begin();
         point_iter != points.end();
         ++point_iter) {
        const DBCALPoint &point = **point_iter;
        int module = point.module();
        int layer = point.layer();
        int sector = point.sector();
        double r = point.r();
        double phi = point.phi();
        double x = r*cos(phi);
        double y = r*sin(phi);

        xcel[module-1][layer-1][sector-1] = x;
        ycel[module-1][layer-1][sector-1] = y;
        //This factory expects z values relative to the center of the BCAL (z=212 cm)
        //but DBCALPoint_factory gives us z values relative to the center of the target
        zcel[module-1][layer-1][sector-1] = point.z()+m_z_target_center-zOffset;
        tcel[module-1][layer-1][sector-1] = point.t();
        ecel[module-1][layer-1][sector-1] = point.E();

        //we require knowledge of the times and energies of the individual upstream and downstream hits
        //to get these we need the associated objects
        double EUp=0,EDown=0,tUp=0,tDown=0;
        vector<const DBCALUnifiedHit*> assoc_hits;
        point.Get(assoc_hits);
        for (unsigned int i=0; i<assoc_hits.size(); i++) {
            if (assoc_hits[i]->end == DBCALGeometry::kUpstream) {
                EUp = assoc_hits[i]->E;
                tUp = assoc_hits[i]->t;
            }
            if (assoc_hits[i]->end == DBCALGeometry::kDownstream) {
                EDown = assoc_hits[i]->E;
                tDown = assoc_hits[i]->t;
            }

        }

        ecel_a[module-1][layer-1][sector-1] = EUp;
        //for some reason we need to record the hit time in two different arrays
        tcel_a[module-1][layer-1][sector-1] = tUp;
        tcell_anor[module-1][layer-1][sector-1] = tUp;

        ecel_b[module-1][layer-1][sector-1] = EDown;
        tcel_b[module-1][layer-1][sector-1] = tDown;
        tcell_bnor[module-1][layer-1][sector-1] = tDown;
    }
}


//------------------
// CeleToArray()
//------------------
void DBCALShower_factory_KLOE::CeleToArray(void)
{
    //  THis code is adpapted from kloe code by Chuncheng Xu on June29,2005
    //  The following part is taken from kaloe's clurec_lib.f's subroutine
    //  cele_to_arr.
    //
    //  It's original purpose is to extract the data information from array RW
    //  and IW into array NARR(j,i), CELDATA(j,i), nclus(i),next(i),
    //  and e_cel(i), x_cel(i),y_cel(i), z_cel(i), t_cel(i),ta_cel(i)
    //  and tb_cel(i)
    //  In the above explanationn, i is the index for cell, and different j
    //  is for different component for such cell.
    //  for example, NARR(1,i), NARR(2,i),NARR(3,i) are for module number,
    //  layer number and sector number in halld bcal simulation 
    
    
    //  Now  we assume that the information is stored in arrays ECEL_A(k,i,j)
    //  , ECEL_B(k,i,j), TCEL_A(k,i,j),  TCEL_B(k,i,j), XCEL(k,I,J),YCEL(k,I,J), 
    // ZCEL(K,I,J), TCEL(K,I,J), ECEL(K,I,J),TCELL_ANOR(K,I,J),TCELL_BNOR(K,I,J),
    //  which are passed into here by event.inc through common block.  
    
    //   these values e_cel(i), x_cel(i),y_cel(i), z_cel(i), t_cel(i),ta_cel(i)
    //  and tb_cel(i)) are extremly useful in later's clusterization
    //  and they are passed by common block to precluster subroutine too. (through
    //  clurec_cal.inc)
    // 
    //  we hope this is a good "bridge" from raw data to precluster. 
    //  This way we can keep good use of their strategy of clusterization
    //  as much as possible, although we know that Fortran has it's bad fame 
    //  of not so easy to be reused.


    // Essentially this function takes the cell information (e.g. E,x,y,z,t)
    // from the 3D arrays filled in CellRecon() and puts it into 1D arrays.
    // These 1D arrays are indexed by an identifier. The module/layer/sector
    // associated with this identifier can be found using the narr[][] array,
    // as described above.
    
    celtot=0;
    
    for (int k = 0; k < modulemax_bcal; k++){
        for (int i = 0; i < layermax_bcal; i++){
            for (int j = 0; j < colmax_bcal; j++){     
                
                float   ea  = ecel_a[k][i][j];
                float   eb  = ecel_b[k][i][j];
                float   ta  = tcel_a[k][i][j];
                float   tb  = tcel_b[k][i][j];
                
                if( (min(ea,eb)>ethr_cell) & (fabs(ta-tb)<35.) & (ta!=0.) & (tb!=0.)) { 
		  celtot=celtot+1;             
		} else {
		  continue;
                }
                
                
                if(celtot>cellmax_bcal) {
                    break;
                }
                
                narr[1][celtot]=k+1;    // these numbers will
                narr[2][celtot]=i+1;    //  be used by preclusters
                narr[3][celtot]=j+1;    //  which will start from index of 1
                                        // rather than from 0.
                                           
                // why 0.145? -- these variables are used as weights
                celdata[1][celtot]=ea/0.145;
                celdata[2][celtot]=eb/0.145;
                
                nclus[celtot] = celtot;
                next[celtot]  = celtot;
                
                e_cel[celtot] = ecel[k][i][j];
                x_cel[celtot] = xcel[k][i][j];
                y_cel[celtot] = ycel[k][i][j];
                z_cel[celtot] = zcel[k][i][j];
                t_cel[celtot] = tcel[k][i][j];
                
                ta_cel[celtot]=tcell_anor[k][i][j];
                tb_cel[celtot]=tcell_bnor[k][i][j];
            }
        }
    }    
}



//------------------
// PreCluster()
//------------------        
void DBCALShower_factory_KLOE::PreCluster(JEventLoop *loop)
{
  //what this function does: for each cell with a hit it finds the maximum
  //energy neighbor and Connect()'s the two. Two cells are neighbors if they
  //are within one column of each other and within one row. The situation is
  //slightly more complicated for two cells on opposite sides of the boundary
  //between inner cells and outer and is described in more detail below, but
  //essentially works out the same way. The purpose of Connect() is described
  //in that function itself.

  int k=1;     // NUMBER OF NEARBY ROWS &/OR TO LOOK FOR MAX E CELL
    
  // extract the BCAL Geometry
  vector<const DBCALGeometry*> bcalGeomVect;
  loop->Get( bcalGeomVect );
  const DBCALGeometry& bcalGeom = *(bcalGeomVect[0]);
    
  // calculate cell position
    
  int   modmin = 0;
  int   modmax = bcalGeom.NBCALMODS;


  //these values make sense as actually minima/maxima if it is implied that rowmin1=1,colmin1=1,colmin2=1
  int   rowmax1= bcalGeom.NBCALLAYSIN;
  int   rowmin2= rowmax1+1;
  //int   rowmax2= bcalGeom.NBCALLAYSOUT+rowmin2-1;
  int   colmax1=bcalGeom.NBCALSECSIN;
  int   colmax2=bcalGeom.NBCALSECSOUT;

  float r_middle= bcalGeom.BCALMIDRAD;

  //radial size of the outermost inner layer
  float thick_inner=bcalGeom.rSize(bcalGeom.cellId(1,bcalGeom.NBCALLAYSIN,1));
  //radial size of the innermost outer layer
  float thick_outer=bcalGeom.rSize(bcalGeom.cellId(1,bcalGeom.NBCALLAYSIN+1,1));

  // this is the radial distance between the center of the innermost outer layer and the outermost inner layer
  float dis_in_out=bcalGeom.r(bcalGeom.cellId(1,bcalGeom.NBCALLAYSIN+1,1))-bcalGeom.r(bcalGeom.cellId(1,bcalGeom.NBCALLAYSIN,1));
    
  float degree_permodule=360.0/(modmax-modmin);
  float half_degree_permodule=degree_permodule/2.0;

  //roughly the width of a single cell in the outermost inner layer
  float width_1=2.0*(r_middle-thick_inner/2.0)*
    sin(half_degree_permodule*3.141593/180)/colmax1;
  //roughly the width of a single cell in the innermost outer layer
  float width_2=2.0*(r_middle+thick_outer/2.0)*
    sin(half_degree_permodule*3.141593/180)/colmax2;

  //disthres is roughly the azimuthal distance between the center of an outer cell and the center of the most distant inner cell bordering an adjacent outer cell (a picture would be nice wouldn't it)
  //this value is used for determining if two cells that straddle the boundary between inner and outer layers should be considered as neighboring
  float disthres=width_2*1.5-width_1*0.5+0.0001;
    
  for (int i = 1; i < (celtot+1); i++){
        
    int maxnn=0; //cell index of the maximum energy neighbor (if one is found)
    float emin=0.; //energy of maximum energy neighbor
        
        
    for (int j = 1; j < (celtot+1); j++){
      if ( (j!=i) & (nclus[j]!=nclus[i]) & (e_cel[j]>emin)) {
                
        int k1= narr[1][i];
        int k2= narr[1][j];
        int i1= narr[2][i];
        int i2= narr[2][j];


        int  modiff = k1-k2;
        int amodif = abs(modiff);
                
        //  the following if is to check module and row distance.         
        if ( (abs(i1-i2)<=k) & ((amodif<=1) || (amodif==47)) ) { 
          //   further check col distance 
          int   j1= narr[3][i];
          int   j2= narr[3][j];

          if(amodif==0) {   // same module       
            //   further check col distance if both are inner layers
            if ( (i1<=rowmax1) & (i2<=rowmax1) & (abs(j2-j1)<=k) ) {
              emin=e_cel[j];
              maxnn=j;
            }

            //   further check col distance if both are outer layers
 
            if ( (i1>=rowmin2) & (i2>=rowmin2) & (abs(j2-j1)<=k) ) {
              emin=e_cel[j];
              maxnn=j;
            }
          }

          if(amodif>0) {  // different module          
            if( (modiff==1) || (modiff==-47) ) {      
              if ( (i1<=rowmax1) & (i2<=rowmax1) ){ 
                if(abs((j1+colmax1)-j2)<=k){
                  emin=e_cel[j];
                  maxnn=j;
                }
              }
                        
              if ( (i1>=rowmin2) & (i2>=rowmin2) ) {
                if(abs((j1+colmax2)-j2)<=k){
                  emin=e_cel[j];
                  maxnn=j;
                }
              }              
            }

            if ( (modiff==-1) || (modiff==47) ) {      

              if ( (i1<=rowmax1) & (i2<=rowmax1) ){
                if(abs((j2+colmax1)-j1)<=k){
                  emin=e_cel[j];
                  maxnn=j;
                }
              } 
                        
              if ( (i1>=rowmin2) & (i2>=rowmin2) ){
                if(abs((j2+colmax2)-j1)<=k){
                  emin=e_cel[j];
                  maxnn=j;
                }
              }              
            }
          }
                    
          // further check col distance if one is inner layer, another is outer
          // so that the two may be between the boundary of two different size
          // of cells.
          if( ( (i1 == rowmax1) & (i2 == rowmin2) ) || 
              ( (i1 == rowmin2) & (i2 == rowmax1) ) ) {

            float delta_xx=xx[k1-1][i1-1][j1-1]-xx[k2-1][i2-1][j2-1];
            float delta_yy=yy[k1-1][i1-1][j1-1]-yy[k2-1][i2-1][j2-1];

            //distance between centers of two cells
            float dis = sqrt( delta_xx * delta_xx + delta_yy * delta_yy );
            //dis_in_out is the distance in radial direction, so we now isolate distance in direction perpendicular to radius
            dis = sqrt( dis*dis - dis_in_out * dis_in_out );
            //disthres is described above
            if( dis < disthres ){
              emin = e_cel[j];
              maxnn = j;
            }
          }                    
        }             
      }
    }        // finish second loop

    if(maxnn>0){

      Connect(maxnn,i);
    }
  }       // finish first loop
}



//------------------
// Connect(n,m);
//------------------   
void DBCALShower_factory_KLOE::Connect(int n,int m)
{
    //----------------------------------------------------------------------
    //   Purpose and Methods :CONNECTS CELL M TO THE NEAREST MAX NEIGHBOR N
    //   Created  31-JUL-1992   WON KIM
    //----------------------------------------------------------------------
    
    // This little piece of code wasn't so easy to decipher, but I *think* I
    // understand what it's doing. The idea is the following:
    //
    // (Prior to entering this routine)
    // The sparsified list of hits is copied into some 1-D arrays with
    // dimension cellmax_bcal+1. This includes the nclus[] and next[]
    // arrays. The nclus[] array keeps the cluster number which is initalized
    // to the sparsified hit number. Thus, every (double-ended) hit cell
    // is it's own cluster. A loop over pairs of hits is done above to find
    // nearest neighbors that should be merged into the same cluster.
    //
    // The next[] array contains an index to the "next" cell in the cluster
    // for the given hit cell. This is a circular list such that if one follows
    // the "next[next[next[...]]]" values, the last element points back to 
    // the first element of the cluster. As such, these are also initialized
    // to index themselves as all single hits are considered a cluster of
    // one element prior to calling this Connect() routine.
    //
    // When we are called, the value of "m" will be less than than value
    // of "n". The cluster number (kept in nclus[]) is therefore updated
    // for all members of nclus[m] to be the same as nclus[n]. Furthermore,
    // the hit cell "m" is appended to the cluster nclus[n]. This also means
    // that the cluster numbers become non-sequential since all elements of
    // nclus whose value is "m" will be changed such that their values are
    // "n".
    //
    // 5/10/2011 DL
    
    if(nclus[n]!=nclus[m]){
        int j=m;
        nclus[j]=nclus[n];
        while(next[j]!=m){
            j=next[j];
            nclus[j]=nclus[n];
	    }        
        next[j]=next[n];
        next[n]=m;
    }
}


//------------------
// ClusNorm()
//------------------   
void DBCALShower_factory_KLOE::ClusNorm(void)
{    
    // fast initialization of arrays:
    memset( e_cls,  0, ( clsmax_bcal + 1 ) * sizeof( float ) );
    memset( x_cls,  0, ( clsmax_bcal + 1 ) * sizeof( float ) );
    memset( y_cls,  0, ( clsmax_bcal + 1 ) * sizeof( float ) );
    memset( z_cls,  0, ( clsmax_bcal + 1 ) * sizeof( float ) );
    memset( t_cls,  0, ( clsmax_bcal + 1 ) * sizeof( float ) );
    memset( ea_cls, 0, ( clsmax_bcal + 1 ) * sizeof( float ) );
    memset( eb_cls, 0, ( clsmax_bcal + 1 ) * sizeof( float ) );
    memset( ta_cls, 0, ( clsmax_bcal + 1 ) * sizeof( float ) );
    memset( tb_cls, 0, ( clsmax_bcal + 1 ) * sizeof( float ) );
    memset( tsqr_a, 0, ( clsmax_bcal + 1 ) * sizeof( float ) );
    memset( tsqr_b, 0, ( clsmax_bcal + 1 ) * sizeof( float ) );
    memset( trms_a, 0, ( clsmax_bcal + 1 ) * sizeof( float ) );
    memset( trms_b, 0, ( clsmax_bcal + 1 ) * sizeof( float ) );
    memset( e2_a,   0, ( clsmax_bcal + 1 ) * sizeof( float ) );
    memset( e2_b,   0, ( clsmax_bcal + 1 ) * sizeof( float ) );
    memset( clspoi, 0, ( clsmax_bcal + 1 ) * sizeof( float ) );
    memset( ncltot, 0, ( clsmax_bcal + 1 ) * sizeof( float ) );
    memset( ntopol, 0, ( clsmax_bcal + 1 ) * sizeof( float ) );
    
    // Part of what is being done here is to further sparsify the
    // data into a list of clusters. This starts to fill arrays
    // with dimension clsmax_bcal+1. One important thing is that
    // the clspoi[] array is being filled with the cluster number
    // of what should be the index of the first cell hit in the
    // cluster. 
    //
    // 5/10/2011 DL

    clstot=0;
    
    for (int ix = 1; ix < (celtot+1); ix++){
        
        //----------------------------------------------------------------------
        // KEEP TALLY OF CLUSTERS
        //----------------------------------------------------------------------
        
        int n=nclus[ix];
        int j=0;
        
        for (int i = 1; i < (clstot+1); i++){
            if(n==clspoi[i]) j=i;
        }
        
        if(j==0) {
            clstot=clstot+1;
            clspoi[clstot]=n;
        }
        
        //----------------------------------------------------------------------
        // NORMALIZED QUANTITIES OF THE TOTAL CLUSTER
        //----------------------------------------------------------------------
        if(e_cel[ix]<0.000000001)continue; // just for protection
        
        x_cls[n]=(e_cls[n]*x_cls[n]+e_cel[ix]*x_cel[ix])
            /(e_cls[n]+e_cel[ix]);
        
        y_cls[n]=(e_cls[n]*y_cls[n]+e_cel[ix]*y_cel[ix])
            /(e_cls[n]+e_cel[ix]);
        
        z_cls[n]=(e_cls[n]*z_cls[n]+e_cel[ix]*z_cel[ix])
            /(e_cls[n]+e_cel[ix]);
        
        t_cls[n]=(e_cls[n]*t_cls[n]+e_cel[ix]*t_cel[ix])
            /(e_cls[n]+e_cel[ix]);
        
        e_cls[n]=e_cls[n]+e_cel[ix];
        
        //        write(*,*)' x-y-z-T-e done'
        //----------------------------------------------------------------------
        //       NORMALIZED QUANTITIES FOR EACH SIDE
        //----------------------------------------------------------------------
        
        ta_cls[n]=(ea_cls[n]*ta_cls[n]+celdata[1][ix]*ta_cel[ix])
            /(ea_cls[n]+celdata[1][ix]);        
        tsqr_a[n]=(ea_cls[n]*tsqr_a[n]+celdata[1][ix]*ta_cel[ix]*
                   ta_cel[ix])/(ea_cls[n]+celdata[1][ix]);
        ea_cls[n]=ea_cls[n]+celdata[1][ix];
        e2_a[n]=e2_a[n]+celdata[1][ix]*celdata[1][ix];
        tb_cls[n]=(eb_cls[n]*tb_cls[n]+celdata[2][ix]*tb_cel[ix])/
            (eb_cls[n]+celdata[2][ix]);
        tsqr_b[n]=(eb_cls[n]*tsqr_b[n]+celdata[2][ix]*tb_cel[ix]*
                   tb_cel[ix])/(eb_cls[n]+celdata[2][ix]);
        eb_cls[n]=eb_cls[n]+celdata[2][ix];
        e2_b[n]=e2_b[n]+celdata[2][ix]*celdata[2][ix];
        
        //----------------------------------------------------------------------
        //       TOTAL CELLS AND CELLS IN OTHER MODULES
        //----------------------------------------------------------------------
        
        ncltot[n]++;

        if( narr[1][n] != narr[1][ix] || narr[2][n] != narr[2][ix] )
            ntopol[n]++;
        
        //----------------------------------------------------------------------
        
    }
    
    for (int n = 1; n < (clstot+1); n++){
        
        int ix=clspoi[n];
        if( ncltot[ix] > 1) {
            
            float effnum = ea_cls[ix] * ea_cls[ix] / e2_a[ix];
            trms_a[ix] = effnum / ( effnum - 1 ) * 
                ( tsqr_a[ix] - ta_cls[ix] * ta_cls[ix] );

            effnum = eb_cls[ix] * eb_cls[ix] / e2_b[ix];
            trms_b[ix] = effnum / ( effnum - 1 ) * 
                ( tsqr_b[ix] - tb_cls[ix] * tb_cls[ix] );

            if( trms_a[ix] <= 0.0 ) trms_a[ix] = 0.;
            if( trms_b[ix] <= 0.0 ) trms_b[ix] = 0.;
            trms_a[ix] = sqrt( trms_a[ix] );
            trms_b[ix] = sqrt( trms_b[ix] );
        }
        else {
            trms_a[ix] = 0.;
            trms_b[ix] = 0.;
        }
    }
}

//------------------
// ClusAnalysis()
//------------------  
void DBCALShower_factory_KLOE::ClusAnalysis()
{
    // track when clusters change to cut down
    // on excess calles to expensive ClusNorm
    bool newClust = false;

    //----------------------------------------------------------------------
    //  Check for overlapping clusters
    //----------------------------------------------------------------------
    
    for (int i = 0; i < 2; i++){
        for (int j = 1; j < (clstot+1); j++){
            int ix=clspoi[j];
            if(e_cls[ix]>0.0){
                float  dist=sqrt(trms_a[ix]*trms_a[ix]+trms_b[ix]*trms_b[ix]);
                if(dist>BREAK_THRESH_TRMS) {
                    Clus_Break(ix);
                    newClust = true;
                }
            }
        }
        
        if( newClust ){
            
            ClusNorm();
            newClust = false;
        }
    }
    
    //----------------------------------------------------------------------       
    // merge clusters likely to be from the same shower
    //----------------------------------------------------------------------       

    int icls[3];
    for (int i = 1; i < clstot; i++){
        icls[1]=0;
        icls[2]=0;
        for (int j = (i+1); j < (clstot+1); j++){
            
            int ix=clspoi[i];
            int iy=clspoi[j];
            
            
            if ( (e_cls[ix]>0.0) & (e_cls[iy]>0.0) ) {
                
                float delta_x=x_cls[ix]-x_cls[iy];
                float delta_y=y_cls[ix]-y_cls[iy];          
                float delta_z=z_cls[ix]-z_cls[iy];  
                float dist=sqrt(delta_x*delta_x+delta_y*delta_y+delta_z*delta_z);
                
                float  tdif=fabs(t_cls[ix]-t_cls[iy]);
                
                //	  float  distz=abs(z_cls[ix]-z_cls[iy]);
                //         cout<<"dist="<<dist<<" distz="<<distz<<" tdif="<<tdif<<"\n";    
                
                if ( (dist<MERGE_THRESH_DIST) & (tdif<MERGE_THRESH_TIME) ){
                    float zdif=fabs(z_cls[ix]-z_cls[iy]);
                    float distran=sqrt(delta_x*delta_x+delta_y*delta_y);
                    
                    if ( (zdif<MERGE_THRESH_ZDIST) & (distran<MERGE_THRESH_XYDIST) ){
                        if(e_cls[ix]>=e_cls[iy]) {
                            icls[1]=ix;
                            icls[2]=iy;
                        }
                        else {
                            icls[1]=iy;
                            icls[2]=ix;
                        }
                    }
                }
                
            }
            
            if(min(icls[1],icls[2])>0){
                
                Connect(icls[1],icls[2]);
                newClust = true;
            }
        }
    }
    
    if( newClust ){
        
        ClusNorm();
    }
}

//------------------
// Clus_Break(ix);
//------------------
void DBCALShower_factory_KLOE::Clus_Break(int nclust)
{
    int   nseed[5],selnum,selcel[cellmax_bcal+1];
    float tdif,tdif_a,tdif_b,tseed[5];  
    
    //----------------------------------------------------------------------
    for (int i =0; i < 5; i++){
        nseed[i]=0;
        tseed[i]=0;
    }
    //----------------------------------------------------------------------
    //   Divide cluster cells into four quadrant groups
    //----------------------------------------------------------------------
    int n=nclust;
    tdif_a=ta_cel[n]-ta_cls[nclust];
    tdif_b=tb_cel[n]-tb_cls[nclust];
    selnum=0;
    
    //----------------------------------------------------------------------
    if(tdif_a>0.0) {
        if(tdif_b>0){
            selnum=1;
        }
        else {
            selnum=2;
        }
    }
    
    else {
        if(tdif_b>0.0) {
            selnum=3;
        }
        else {
            selnum=4;
        }
    }
    
    //------------------------------------------------------------------------
    
    if(selnum>0) {
        float tdif=sqrt(tdif_a*tdif_a+tdif_b*tdif_b);
        if(tdif>tseed[selnum]){
            nseed[selnum]=n;
            tseed[selnum]=tdif;
        }
        selcel[n]=selnum;
    }
    
    //-------------------------------------------------------------------------
    
    while(next[n]!=nclust) {
        n=next[n];
        tdif_a=ta_cel[n]-ta_cls[nclust];
        tdif_b=tb_cel[n]-tb_cls[nclust];
        selnum=0;
        
        //**************************************************************************
        
        if(tdif_a>0.0) {
            
            if(tdif_b>0.0) {
                selnum=1;
            }
            else {
                selnum=2;
            }
        }
        
        else {
            if(tdif_b>0.0) {
                selnum=3;
            } 
            else {
                selnum=4;
            }
        }
        
        //***************************************************************************
        
        //...........................................................................
        
        
        if(selnum>0){
            tdif=sqrt(tdif_a*tdif_a+tdif_b*tdif_b);
            
            if(tdif>tseed[selnum]){
                nseed[selnum]=n;
                tseed[selnum]=tdif;
            }
            
            selcel[n]=selnum;
        }
        
        //...........................................................................
        
    }             //this bracket is related to "while" above. 
    
    
    //----------------------------------------------------------------------
    // If successful, divide cluster chain into the new cluster chains
    //----------------------------------------------------------------------
    
    for (int i =1; i < 5; i++){
        
        if(nseed[i]>0) {
            
            nclus[nseed[i]]=nseed[i];
            next[nseed[i]]=nseed[i];

            for (int j =1; j < (celtot+1); j++){
	      if ( (nclus[j]==nclust) & (j!=nseed[i]) ){
                    if(selcel[j]==i) {
                        nclus[j]=j;
                        next[j]=j;
                        Connect(nseed[i],j);
                    }
                }
            }
        }
    }
}


//------------------
// Trakfit()
//------------------ 

// routine modified by MRS to do expensive filling of layer information
// contains code previously in ClusNorm which is called multiple times
// per event, but there is no need to call Trakfit multiple times per event
void DBCALShower_factory_KLOE::Trakfit( void )
{
 
    float emin=0.0001;

    memset( clslyr, 0, ( clsmax_bcal + 1 ) * 
            ( layermax_bcal + 1 ) * 6 * sizeof( float ) );
    
    memset( apx,   0, ( clsmax_bcal + 1 ) * 4 * sizeof( float ) );
    memset( eapx,  0, ( clsmax_bcal + 1 ) * 4 * sizeof( float ) );
    memset( ctrk,  0, ( clsmax_bcal + 1 ) * 4 * sizeof( float ) );
    memset( ectrk, 0, ( clsmax_bcal + 1 ) * 4 * sizeof( float ) );
    
    for (int ix = 1; ix < (celtot+1); ix++){
        
        int n = nclus[ix];
        
        //----------------------------------------------------------------------
        // NORMALIZED QUANTITIES OF THE TOTAL CLUSTER PER LAYER
        //----------------------------------------------------------------------
        
        int lyr = narr[2][ix];
        
        //	write(*,*)'X, E',E_CEL(ix),CLSLYR(XLYR,LYR,N)
        clslyr[xlyr][lyr][n]=(clslyr[elyr][lyr][n]*clslyr[xlyr][lyr][n]
                              +e_cel[ix]*x_cel[ix])/(clslyr[elyr][lyr][n]+e_cel[ix]);
        
        clslyr[ylyr][lyr][n]=(clslyr[elyr][lyr][n]*clslyr[ylyr][lyr][n]
                              +e_cel[ix]*y_cel[ix])/(clslyr[elyr][lyr][n]+e_cel[ix]);
        //	write(*,*)' Y'
        
        clslyr[zlyr][lyr][n]=(clslyr[elyr][lyr][n]*clslyr[zlyr][lyr][n]
                              +e_cel[ix]*z_cel[ix])/(clslyr[elyr][lyr][n]+e_cel[ix]);
        
        //	write(*,*)' z'
        clslyr[tlyr][lyr][n]=(clslyr[elyr][lyr][n]*clslyr[tlyr][lyr][n]
                              +e_cel[ix]*t_cel[ix])/(clslyr[elyr][lyr][n]+e_cel[ix]);
        
        //	write(*,*)'E'
        clslyr[elyr][lyr][n]=clslyr[elyr][lyr][n]+e_cel[ix];
        
        //        write(*,*)' slopes done'
    }

    memset( nlrtot, 0, ( clsmax_bcal + 1 ) * sizeof( float ) );

    for (int n = 1; n < ( clstot + 1 ); n++){

        int ix=clspoi[n];
        
        for (int i = 1; i < (layermax_bcal+1); i++){
            
            if( clslyr[elyr][i][ix] > 0.0 ) nlrtot[ix]++;
        }
                    
        for (int i = 0; i < ( layermax_bcal + 1 ); i++){
        
            x[i]=0.0;
            y[i]=0.0;  
            z[i]=0.0;
            e[i]=0.0;
            sigx[i]=0.0;  // cm
            sigy[i]=0.0; // cm
            sigz[i]=0.0;
        }
        
        int nltot=0;
    
        for (int il = 1; il < (layermax_bcal+1); il++){

            if(clslyr[elyr][il][ix]>emin) {

                nltot=nltot+1;
                x[nltot]= clslyr[xlyr][il][ix];
                y[nltot]= clslyr[ylyr][il][ix];
                z[nltot]= clslyr[zlyr][il][ix];         
                e[nltot]= clslyr[elyr][il][ix]; 

                sigy[nltot] =  1.0/e[nltot];
                sigx[nltot] =  1.0/e[nltot];
                sigz[nltot] =  1.0/sqrt(e[nltot]);
            }
        }
    
        // The following error bar is the estimation of error bar
        // based on the experience of my fortran code running
        // If the structure of BCAL changes drastically, you have to make changes
        // accordingly.  Xu Chuncheng , Jan 9, 2006
        
        // **** this seems strange since it appears to overwrite what is above
    
        sigx[1]=0.5;
        sigx[2]=0.5;
        sigx[3]=0.5;
        sigx[4]=0.5;
        sigx[5]=0.5;
        sigx[6]=0.8;
        sigx[7]=0.9;
        sigx[8]=1.2;
        sigx[9]=1.3;
        
        sigy[1]=0.5;
        sigy[2]=0.5;
        sigy[3]=0.5;
        sigy[4]=0.5;
        sigy[5]=0.5;
        sigy[6]=0.8;
        sigy[7]=0.9;
        sigy[8]=1.2;
        sigy[9]=1.3;
        
        
        sigz[1]=0.5;
        sigz[2]=0.5;
        sigz[3]=0.5;
        sigz[4]=0.5;
        sigz[5]=0.5;
        sigz[6]=0.8;
        sigz[7]=0.9;
        sigz[8]=1.2;
        sigz[9]=1.3;
        
        if( nltot > 1 ){ 
            
            Fit_ls();
            
            for (int i = 1; i < 4; i++){
                
                ctrk[i][ix]=ctrk_ix[i];
                ectrk[i][ix]=ectrk_ix[i];
                apx[i][ix]=apx_ix[i];
                eapx[i][ix]=eapx_ix[i];
            }
        }    
        else{
            
            // we have 1 or less layers hit -- no fit
            
            apx[1][ix] = x[1];
            apx[2][ix] = y[1];
            apx[3][ix] = z[1]; 
            eapx[1][ix] = sigx[1];
            eapx[2][ix] = sigy[1];
            eapx[3][ix] = sigz[1];  
            ectrk[1][ix] = 0.0;
            ectrk[2][ix] = 0.0;
            ectrk[3][ix] = 0.0;
            ctrk[1][ix] = 999.0;
            ctrk[2][ix] = 999.0;
            ctrk[3][ix] = 999.0;
        }            
    }
}

//------------------
// Fit_ls()
//------------------  
void DBCALShower_factory_KLOE::Fit_ls()
{      
    float a,b,c;
    float d,e,f,chi2,q,norm;
    float siga,sigb,sigc,sigd,sige,sigf;
    float sigb2,sigd2,sigf2;
    
    //    fitting for X=a+bt
    Linefit(1,1,a,b,siga,sigb,chi2,q);
    //    fitting for Y=c+dt
    Linefit(2,1,c,d,sigc,sigd,chi2,q);
    //    fitting for Z=e+ft
    Linefit(3,1,e,f,sige,sigf,chi2,q);
    sigb2=sigb*sigb;
    sigd2=sigd*sigd;
    sigf2=sigf*sigf;
    
    apx_ix[1]=a;
    apx_ix[2]=c;
    apx_ix[3]=e;
    eapx_ix[1]=siga;
    eapx_ix[2]=sigc;
    eapx_ix[3]=sige;
    
    //      write(*,*) "b,d,f = ",b,d,f
    //		     cout<<"b= "<<b<<" d="<<d<<" f="<<f<<"\n";
    
    //      write(*,*) "sigb,sigd,sigf = ",sigb,sigd,sigf
    
    norm=sqrt(b*b+d*d+f*f);
    
    ctrk_ix[1]=b/norm;
    ctrk_ix[2]=d/norm;
    ctrk_ix[3]=f/norm;
    
    float norm3=norm*norm*norm;
    
    ectrk_ix[1]=sqrt((d*d+f*f)*(d*d+f*f)*sigb2+b*b*d*d*sigd2+b*b*f*f*sigf2)/norm3;
    ectrk_ix[2]=sqrt((b*b+f*f)*(b*b+f*f)*sigd2+d*d*b*b*sigb2+d*d*f*f*sigf2)/norm3;
    ectrk_ix[3]=sqrt((b*b+d*d)*(b*b+d*d)*sigf2+f*f*b*b*sigb2+f*f*d*d*sigd2)/norm3;
    
    return;
}


//------------------
// Linefit(x,y,ndata,sig,mwt,a,b,siga,sigb,chi2,q)
//------------------  
void DBCALShower_factory_KLOE::Linefit(int ixyz,int mwt,float &a,
                                  float &b,float &siga,float &sigb,float &chi2,float &q)
{
    
    // This programme is taken from Garth's book
    //  "Numerical Recipes in Fortran"
    //  and we tested it. It works well.  Chuncheng Xu,2005
    
    
    float sig[layermax_bcal+1],etemp;
    float xtemp[layermax_bcal+1],ytemp[layermax_bcal+1];
    // uses gammq
    
    //  Given a set of data points X(1:ndata),Y(1:ndata) with individual standard
    //  deviations sig(1:ndata), fit them to a strait line y=a+bx by minimum 
    // chisq. Returned are a,b and their respective probable uncertainties
    // siga and sigb, the chisq, and the goodness-of-fit probability q(that the
    // fit would have chisq this large or larger). if mwt=0 on input, then the 
    // the standard deviations are assumed to be unavalable:q is returned as 
    // 1.0 and the normalization of chi2 is to the unit standard deviation on all
    // points.
    
    
    float sigdat,ss,st2,sx,sxoss,sy,t,wt;
    sx=0.0;             // Initialize sums to zero
    sy=0.0;
    st2=0.0;
    b=0.0;
    
    int ndata=0;
    
    
    
    if(ixyz==1) {
        for (int i = 1; i < (layermax_bcal+1); i++){
            xtemp[i]=rt[i];
            ytemp[i]=x[i];
            sig[i]=sigx[i];
            etemp=e[i];
            if(etemp>0.0001)ndata=ndata+1;
            
        }   
    }
    else if(ixyz==2) {
        for (int i = 1; i < (layermax_bcal+1); i++){
            xtemp[i]=rt[i];
            ytemp[i]=y[i];
            sig[i]=sigy[i];
            etemp=e[i];
            if(etemp>0.000001)ndata=ndata+1;   
        }
    }
    else if(ixyz==3) {
        for (unsigned int i = 1; i < (layermax_bcal+1); i++){
            xtemp[i]=rt[i];
            ytemp[i]=z[i]; 
            sig[i]=sigz[i];  
            etemp=e[i];
            if(etemp>0.000001)ndata=ndata+1;
        }
    }
    
    if(mwt!=0) {   // Accumulate sums
        ss=0.0;
        for (int i = 1; i < (ndata+1); i++){
            wt=1.0/(sig[i]*sig[i]);
            ss=ss+wt;
            sx=sx+xtemp[i]*wt;
            sy=sy+ytemp[i]*wt;
        }
    }
    
    else  {
        for (int i = 1; i < (ndata+1); i++){
            sx=sx+xtemp[i];
            sy=sy+ytemp[i];
        } 
        ss=float(ndata);
    }
    
    sxoss=sx/ss;
    
    if(mwt!=0) {
        for (int i = 1; i < (ndata+1); i++){
            t=(xtemp[i]-sxoss)/sig[i];
            st2=st2+t*t;
            b=b+t*ytemp[i]/sig[i];
        }
        
    }
    
    else  {
        
        for (int i = 1; i < (ndata+1); i++){
            t=xtemp[i]-sxoss;
            st2=st2+t*t;
            b=b+t*ytemp[i];
        }
    }
    
    b=b/st2;                // Solve for a,b,siga and sigb
    a=(sy-sx*b)/ss;
    siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
    sigb=sqrt(1.0/st2);
    chi2=0.0;        //   Calculate chisq
    
    if(mwt==0) {
        for (int i = 1; i < (ndata+1); i++){
            chi2=chi2+(ytemp[i]-a-b*xtemp[i])*(ytemp[i]-a-b*xtemp[i]);
        }
        q=1.0;
        sigdat=sqrt(chi2/(ndata-2));
        siga=siga*sigdat;
        sigb=sigb*sigdat;
    }
    else  {
        for (int i = 1; i < (ndata+1); i++){
            chi2=chi2+((ytemp[i]-a-b*xtemp[i])/
                       sig[i])*((ytemp[i]-a-b*xtemp[i])/sig[i]);
        }
        q=Gammq(0.5*(ndata-2),0.5*chi2); 
    }
    
}



//------------------
// Gammq(a_gammq,x_gammq)
//------------------  
float DBCALShower_factory_KLOE::Gammq(float a_gammq,float x_gammq)
{
    
    //============================================================================
    
    
    float gammq;
    
    //    uses gcf,gser
    //Returns the incomplete gamma function Q(a_gammq,x_gammq)=1-P(a_gammq,x_gammq)
    
    
    float gammcf,gamser;
    
    if(a_gammq==0.0) { 
        gammq=999.0;
        return gammq;
    }
    
    
    
    if(x_gammq<0. || a_gammq<= 0.0) {
        //  cout<<"bad arguments in gammq"<<"\n";
        return 999.0;// pause 'bad arguments in gammq'
    }
    
    if(x_gammq<(a_gammq+1.)) {     // Use the series representation
        Gser(gamser,a_gammq,x_gammq);
        gammq=1.0-gamser;    // and take its complement
    }
    else  {               // Use continued fraction representation
        Gcf(gammcf,a_gammq,x_gammq);
        gammq=gammcf;
    }
    return gammq;
}


//------------------
// Gser(gamser,a_gser,x_gser,gln)
//------------------  
void DBCALShower_factory_KLOE::Gser(float &gamser,float a_gser,float x_gser)
{
    //============================================================================
    int itmax=100;
    float eps=3.0e-7;
    float gln;            
    
    //    uses gammln
    //    Returns the incomplete gamma function P(a,x) evaluated by its
    //    series representation as gamser. Also returns ln(Gamma(a)) as gln
    
	   float ap, del,sum;
       
       gln=Gammln(a_gser);
       
       if(x_gser<=0.0) {
           if(x_gser<0.0) cout<<"x_gser<0 in gser"<<"\n";
           gamser=0.0;
           return;
       }
       
       ap=a_gser;
       sum=1.0/a_gser;
       del=sum;
       
       
       for (int n = 1; n < (itmax+1); n++){
           ap=ap+1.0;
           del=del*x_gser/ap;
           sum=sum+del;
           
           if(fabs(del)<fabs(sum)*eps) {
               gamser=sum*exp(-x_gser+a_gser*log(x_gser)-gln);
               return;
           }
           
       }
       
       // cout<< "a too large, ITMAX too small in gser"<<"\n";
       return;
}


//------------------
//  Gcf(gammcf,a,x,gln)
//------------------  
void DBCALShower_factory_KLOE::Gcf(float &gammcf,float a_gcf,float x_gcf)
{
    //========================================================================
    
    int itmax=100;
    float eps=3.0e-7;
    float fpmin=1.0e-30;
    
    float gln;
    
    //   uses gammln
    //        Returns the incomplete gamma function Q(a,x) evaluated by 
    //   its continued fraction representation as gammcf. Also returns
    //   ln(Gamma(a)) as gln
    
    //   Parameters: ITMAX is the maximum allowed number of iterations;
    //   EPS is the relative accuracy; FPMIN is a number near the smallest
    //   representable floating-point number.
    
    float an,b,c,d,del,h;
    
    gln=Gammln(a_gcf);
    b=x_gcf+1.0-a_gcf;
    c=1.0/fpmin;
    d=1.0/b;
    h=d;
    
    
    for (int i = 1; i < (itmax+1); i++){
        an=-i*(i-a_gcf);
        b=b+2.0;
        d=an*d+b;
        if(fabs(d)<fpmin)d=fpmin;
        c=b+an/c;
        if(fabs(c)<fpmin)c=fpmin;
        d=1.0/d;
        del=d*c;
        h=h*del;
        if(fabs(del-1.0)<eps) {	 
            gammcf=exp(-x_gcf+a_gcf*log(x_gcf)-gln)*h; // Put factors in front 
            return;
        }
        // cout<< "a too large,ITMAX too small in gcf"<<"\n";    
        return;
    }
} //???

//------------------
// Gammln(xx)
//------------------  
float DBCALShower_factory_KLOE::Gammln(float xx_gln)
{
    //    Returns the value ln[Gamma(xx_gln)] for xx_gln>0.0
    
	   float ser,stp,tmp,x_gln,y_gln;
    float cof[7];
    float gammln; 
    
    //    Internal arithmetic will be done in double precision, a nicety that
    //    that you can omit if five-figure accuracy is good enoug
    //         save cof,stp
    //         data cof,stp/76.18009172947146d0,-86.50532032941677d0,
    //     * 24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
    //     * -.5395239384953d-5,2.5066282746310005d0/
    
    stp=2.5066282746310005;
    cof[1]=76.18009172947146;
    cof[2]=-86.50532032941677;
    cof[3]=24.01409824083091;
    cof[4]=-1.231739572450155;
    cof[5]=.1208650973866179e-2;
    cof[6]=-.5395239384953e-5;
    
    x_gln=xx_gln;
    y_gln=x_gln;
    tmp=x_gln+5.5;
    tmp=(x_gln+0.5)*log(tmp)-tmp;
    
    ser=1.000000000190015;
    
    for (int j = 1; j < 7; j++){
        y_gln=y_gln+1.0;
        ser=ser+cof[j]/y_gln;
    }
    
    
    gammln=tmp+log(stp*ser/x_gln);
    return gammln;
}

