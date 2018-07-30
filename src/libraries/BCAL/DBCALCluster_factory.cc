/*
 *  DBCALCluster_factory.cc
 *
 *  Created by Matthew Shepherd on 3/12/11.
 *
 */

#include <iostream>

using namespace std;

#include "DANA/DApplication.h"
#include "BCAL/DBCALGeometry.h"
#include "BCAL/DBCALHit.h"
#include "BCAL/DBCALUnifiedHit.h"

#include "BCAL/DBCALCluster_factory.h"

#include "units.h"
#include <TMath.h>

bool PointSort( const DBCALPoint* p1, const DBCALPoint* p2 ){

	return ( p1->E() > p2->E() );
}

bool ClusterSort( const DBCALCluster* c1, const DBCALCluster* c2 ){

	return ( c1->E() > c2->E() );
}

DBCALCluster_factory::DBCALCluster_factory() : 
	m_mergeSig( 5 ), 
	m_moliereRadius( 3.7*k_cm ),
	m_clust_hit_timecut ( 20.0*k_nsec ),
	m_timeCut( 8.0*k_nsec ){

	// The phi and theta direction inclusion curves are described in: 
	// http://argus.phys.uregina.ca/gluex/DocDB/0029/002998/003/CAL_meeting_may5.pdf.
	// The theta direction inclusion curve needs to be a function of theta. C1_parm and
	// C2_parm are parameters [0] and [1] in dtheta_inclusion_curve. 
	}

jerror_t
DBCALCluster_factory::init(void){

	m_BCALGeom = NULL;
	return NOERROR;

}

jerror_t
DBCALCluster_factory::fini( void ){

	return NOERROR;
}

jerror_t DBCALCluster_factory::brun(JEventLoop *loop, int32_t runnumber) {
	DApplication* app = dynamic_cast<DApplication*>(loop->GetJApplication());
	DGeometry* geom = app->GetDGeometry(runnumber);
	geom->GetTargetZ(m_z_target_center);

	// load BCAL Geometry
	vector<const DBCALGeometry *> BCALGeomVec;
        loop->Get(BCALGeomVec);
        if(BCALGeomVec.size() == 0)
                throw JException("Could not load DBCALGeometry object!");
        m_BCALGeom = BCALGeomVec[0];


	loop->GetCalib("/BCAL/effective_velocities", effective_velocities);

	loop->GetCalib("/BCAL/attenuation_parameters",attenuation_parameters);

	BCALCLUSTERVERBOSE = 0;
	gPARMS->SetDefaultParameter("BCALCLUSTERVERBOSE", BCALCLUSTERVERBOSE, "VERBOSE level for BCAL Cluster overlap success and conditions");
	//command line parameter to investigate what points are being added to clusters and what clusters are being merged together. // Track fitterer helper class


  vector<const DTrackFitter *> fitters;
  loop->Get(fitters);

  if(fitters.size()<1){
    _DBG_<<"Unable to get a DTrackFinder object!"<<endl;
    return RESOURCE_UNAVAILABLE;
  }

  fitter = fitters[0];

	return NOERROR;
}

jerror_t
DBCALCluster_factory::evnt( JEventLoop *loop, uint64_t eventnumber ){

	vector< const DBCALPoint* > twoEndPoint;
	vector< const DBCALPoint* > usedPoints;
	loop->Get(twoEndPoint);

	// Want to add singled-ended hits to the Clusters. 

	// Looking for hits that are single-ended.

	vector< const DBCALUnifiedHit* > hits;
	loop->Get(hits);

	vector< const DTrackWireBased* > tracks;
	loop->Get(tracks);

	// first arrange the list of hits so they are grouped by cell
	map< int, vector< const DBCALUnifiedHit* > > cellHitMap;
	for( vector< const DBCALUnifiedHit* >::const_iterator hitPtr = hits.begin();
			hitPtr != hits.end();
			++hitPtr ){

		const DBCALUnifiedHit& hit = (**hitPtr);

		int id = m_BCALGeom->cellId( hit.module, hit.layer, hit.sector );

		if( cellHitMap.find( id ) == cellHitMap.end() ){

			cellHitMap[id] = vector< const DBCALUnifiedHit* >();
		}

		cellHitMap[id].push_back( *hitPtr );  
	}

	// now we should try to add on single-ended hits ... 
	vector< const DBCALUnifiedHit* > single_ended_hits; 

	for( map< int, vector< const DBCALUnifiedHit* > >::iterator mapItr = cellHitMap.begin();
			mapItr != cellHitMap.end();
			++mapItr ){

		if( mapItr->second.size() == 1 ){      
			// only one hit in the cell

			const DBCALUnifiedHit* hit = mapItr->second[0];

			single_ended_hits.push_back(hit);

		}
	}

	vector<DBCALCluster*> clusters = clusterize( twoEndPoint, usedPoints,  single_ended_hits, tracks );

	// load our vector of clusters into the factory member data
	for( vector<DBCALCluster*>::iterator clust = clusters.begin();
			clust != clusters.end();
			++clust ){
		
		if( isnan((**clust).t()) == 1 || isnan((**clust).phi()) == 1 || isnan((**clust).theta()) == 1 ) continue;
		// put in an energy threshold for clusters
		if( (**clust).E() < 5*k_MeV ) {
			delete *clust;
			continue;
		}
		vector<const DBCALPoint*>points=(**clust).points();
		for (unsigned int i=0;i<points.size();i++){
		  (**clust).AddAssociatedObject(points[i]);
		}
		_data.push_back(*clust);
	}
	return NOERROR;
}

vector<DBCALCluster*>
DBCALCluster_factory::clusterize( vector< const DBCALPoint* > points , vector< const DBCALPoint* > usedPoints ,  vector< const DBCALUnifiedHit* > hits, vector< const DTrackWireBased* > tracks ) const {

  // first sort the points by energy
  sort( points.begin(), points.end(), PointSort );
  
  vector<DBCALCluster*> clusters(0);
  
  // ahh.. more hard coded numbers that should probably
  // come from a database or something
  float seedThresh = 1.*k_GeV;
  float minSeed = 10*k_MeV;
  //We have a big problem with noise in the outer layer of the detector
  //(the noise is the greatest in the outer layer, since the number of SiPMs
  //being summed is also the greatest here).
  //Thus there are a lot of DBCALPoint's in this layer that are pure noise hits.
  //The simplest way to deal with this is to prevent outer layer points
  //from seeding clusters. So hits in the outer layer can be associated
  //with existing clusters, but cannot create their own cluster.
  //This is okay since since isolated hits in the outer layer
  //is not really a signature we expect for many physical showers.
  //However, if a hit is sufficiently energetic, it is unlikely to be a noise
  //hit. For this reason, we allow 4th layer hits to seed clusters,
  //but we need a different (higher) minimum seed energy.
  float layer4_minSeed = 50*k_MeV;
  float tracked_phi = 0.;
  float matched_dphi = .175;
  float matched_dtheta = .087; 
  
  while( seedThresh > minSeed ) {
    
    bool usedPoint = false;

    for( vector< const DBCALPoint* >::iterator pt = points.begin();
	 pt != points.end();
	 ++pt ){

      // first see if point should be added to an existing
      // cluster
      
      int q = 0;

      // Check if a point is matched to a track      
      for( vector< const DTrackWireBased* >::iterator trk = tracks.begin();
	   trk != tracks.end();
	   ++trk ){
	DVector3 track_pos(0.0, 0.0, 0.0);
	double point_r = (**pt).r();
	double point_z = (**pt).z();
	vector<DTrackFitter::Extrapolation_t>extrapolations=(*trk)->extrapolations.at(SYS_BCAL);
	if (fitter->ExtrapolateToRadius(point_r,extrapolations,track_pos)){
	  double dPhi=track_pos.Phi()-(**pt).phi();
	  if (dPhi<-M_PI) dPhi+=2.*M_PI;
	  if (dPhi>M_PI) dPhi-=2.*M_PI;
	  double point_theta_global = fabs(atan2(point_r,point_z + m_z_target_center ));  // convert point z-position origin to global frame to match tracks origin
	  double dTheta = fabs(point_theta_global - track_pos.Theta());
	  matched_dphi=0.175+0.175*exp(-0.8*extrapolations[0].momentum.Mag());
	  if(fabs(dPhi) < matched_dphi && dTheta < matched_dtheta){
	    q = 1; // if point and track are matched then set q = 1
	    tracked_phi = extrapolations[0].position.Phi();
	    break;
	  }
	}
      }

      for( vector<DBCALCluster*>::iterator clust = clusters.begin();
	   clust != clusters.end();
	   ++clust ){
	
	if((**clust).Q()==1){
	  if(overlap_charged( **clust,*pt, tracked_phi ) ){
	    usedPoints.push_back( *pt );
	    int point_q = 1;
	    (**clust).addPoint( *pt, point_q );
	    points.erase( pt );
	    usedPoint = true;
	    break;
	  }
	}
	if( overlap( **clust, *pt ) ){
	  if (q==1 && (**pt).layer()!=1) q=0;
	  // assign point q=1 if it's in layer 1 because track matching tends to be improved in layer 1 than later layers where the cluster seed is. This would allow us to jump into the charged clustering routines on the fly.
	  usedPoints.push_back( *pt );  
	  (**clust).addPoint( *pt , q);
	  points.erase( pt );
	  usedPoint = true;
	  break;
	}
	// once we erase a point the iterator is no longer useful
	// and we start the loop over, so that a point doesn't get added to
	// multiple clusters. We will recycle through points later to 
	// check if a point was added to its closest cluster.
      }
    
      if( usedPoint ) break;
      
      // if the point doesn't overlap with a cluster see if it can become a 
      // new seed
      if( (**pt).E() > seedThresh && ((**pt).layer() != 4 || (**pt).E() > layer4_minSeed) ){
	clusters.push_back(new DBCALCluster( *pt, m_z_target_center, q, m_BCALGeom  ) );
	usedPoints.push_back( *pt );
	points.erase( pt );
	usedPoint = true;
	break;
      }
    }

    recycle_points( usedPoints, clusters);	
    // recycle through points that were added to a cluster and check if they
    // were added to their closest cluster. If they weren't then we remove 
    // the point from its original cluster and add it to its closest cluster.
  
    double point_reatten_E = 0.;  
    merge( clusters, point_reatten_E );
    // lower the threshold to look for new seeds if none of 
    // the existing points were used as new clusters or assigned
    // to existing clusters
    if( !usedPoint ) seedThresh /= 2;
  }

  // add the single-ended hits that overlap with a cluster that was made from points
  for( vector< const DBCALUnifiedHit* >::iterator ht = hits.begin();
       ht != hits.end();
       ++ht){
    bool usedHit = false;	 

    for( vector<DBCALCluster*>::iterator clust = clusters.begin();
	 clust != clusters.end();
	 ++clust ){
      
      if( overlap( **clust, *ht ) ){

	int channel_calib = 16*((**ht).module-1)+4*((**ht).layer-1)+(**ht).sector-1; // need to use cellID for objects in DBCALGeometry but the CCDB uses a different channel numbering scheme, so use channel_calib when accessing CCDB tables.

	// given the location of the cluster, we need the best guess
	// for z with respect to target at this radius
	
	double z = (**clust).rho()*cos((**clust).theta()) + m_z_target_center;
	double d = ( ((**ht).end == 0) ? (z  - m_BCALGeom->GetBCAL_center() + m_BCALGeom->GetBCAL_length()/2.0) : (m_BCALGeom->GetBCAL_center() + m_BCALGeom->GetBCAL_length()/2.0 - z));  // d gives the distance to upstream or downstream end of BCAL depending on where the hit was with respect to the cluster z position.
	double lambda = attenuation_parameters[channel_calib][0];
	double hit_E = (**ht).E;
	double hit_E_unattenuated = hit_E/exp(-d/lambda);  // hit energy unattenuated wrt the cluster z position
	
	(**clust).addHit( *ht, hit_E_unattenuated );
	usedHit = true;
      }
      if( usedHit ) break;
    }
  }     
  return clusters;
}

void
DBCALCluster_factory::recycle_points( vector<const DBCALPoint*> usedPoints, vector<DBCALCluster*>& clusters) const{

  if ( clusters.size() <= 1 ) return;

  int q = 2;
	
  sort( clusters.begin(), clusters.end(), ClusterSort );

  for( vector<const DBCALPoint*>::const_iterator usedpt = usedPoints.begin();
       usedpt != usedPoints.end();
       ++usedpt ){		
    
    bool got_overlap=false;				
    double min_phi=1e6;
		
    for( vector<DBCALCluster*>::iterator clust = clusters.begin();
	 clust != clusters.end();
	 ++clust ){
      
      if( overlap( **clust, *usedpt ) ){
	got_overlap=true;

	float deltaPhi = (**clust).phi() - (*usedpt)->phi();
	if (deltaPhi<-M_PI) deltaPhi+=2.*M_PI;
	if (deltaPhi>M_PI) deltaPhi-=2.*M_PI;
	if (fabs(deltaPhi)<min_phi){
	  min_phi=fabs(deltaPhi);
	}
      }
    }
    
    if(got_overlap==false) break;

    // Find the points closest cluster in distance along the sphere and in phi
    for( vector<DBCALCluster*>::iterator clust = clusters.begin();
	 clust != clusters.end();
	 ++clust ){
      bool best_clust = false;
      vector<const DBCALPoint*>associated_points=(**clust).points();

      float deltaPhi = (**clust).phi() - (*usedpt)->phi();
      if (deltaPhi<-M_PI) deltaPhi+=2.*M_PI;
      if (deltaPhi>M_PI) deltaPhi-=2.*M_PI;
      deltaPhi=fabs(deltaPhi);
    
      for(unsigned int j = 0 ; j < associated_points.size(); j++){
	// Check to see if the point we are comparing to the cluster
	// is already in that cluster.
	if (fabs((*usedpt)->E()-associated_points[j]->E())<1e-4
	    && fabs(deltaPhi-min_phi)<1e-4) best_clust=true;
	if(BCALCLUSTERVERBOSE>1)cout << " clust E = " << (**clust).E() <<" assoc point E = " << associated_points[j]->E() << " points E = " << (*usedpt)->E() <<  " clust match = " << best_clust <<  endl;
      }
      if(best_clust==true) break;
      // if the point was originally placed in its "best" cluster then we don't want to touch it.
      if(best_clust==0){
	int added_point = 0;
	int removed_point = 0;
	for(unsigned int i = 0 ; i < associated_points.size(); i++){
	  bool point_match = (fabs((*usedpt)->E()-associated_points[i]->E())<1e-4);
	  if( point_match==0 && added_point==0 && fabs(deltaPhi-min_phi)<1e-4){
	    (**clust).addPoint( *usedpt , q );
	    // if the point found a closer cluster then we add it to the closer cluster.
	    // The point is now an associated object of the closer cluster.
	    added_point=1;
	  }
	  if( point_match==1 && removed_point==0 && fabs(deltaPhi-min_phi)>1e-4){
	    (**clust).removePoint( *usedpt );
	    // Now we remove the point from its original cluster since it has been added
	    // to its closest cluster. The point is no longer an associated object of
	    // the original cluster.
	    removed_point=1;
	  }
	}
      }   
    }
  }
}	

void
DBCALCluster_factory::merge( vector<DBCALCluster*>& clusters, double point_reatten_E ) const {

	if( clusters.size() <= 1 ) return;

	sort( clusters.begin(), clusters.end(), ClusterSort );

	bool stillMerging = true;

	float low_z_lim = -100.;
	float high_z_lim = 500.;

	while( stillMerging ){

		stillMerging = false;
		for( vector<DBCALCluster*>::iterator hClust = clusters.begin();
			hClust != clusters.end() - 1;
			++hClust ){

			vector<const DBCALPoint*>hClust_points=(**hClust).points();

			for( vector<DBCALCluster*>::iterator lClust = hClust + 1;
				lClust != clusters.end();
				++lClust ){

				vector<const DBCALPoint*>lClust_points=(**lClust).points();
				vector<const DBCALPoint*>hClust_points=(**hClust).points();
			
				if( overlap( **hClust, **lClust ) ){

					point_reatten_E = 0.;

					if (hClust_points.size() == 1) {

						for( unsigned int i = 0 ; i < hClust_points.size() ; i++){

							if (hClust_points[i]->z() > low_z_lim && hClust_points[i]->z() < high_z_lim) point_reatten_E = 0.;
							else {
							      int channel_calib = 16*(hClust_points[i]->module()-1)+4*(hClust_points[i]->layer()-1)+hClust_points[i]->sector()-1;

							      double fibLen = m_BCALGeom->GetBCAL_length();
	
							      double point_z = hClust_points[i]->z();
							      double zLocal = point_z + m_z_target_center - m_BCALGeom->GetBCAL_center();

							      double dUp = 0.5 * fibLen + zLocal;
							      double dDown = 0.5 * fibLen - zLocal;
							      if (dUp>fibLen)   dUp=fibLen;
							      if (dUp<0)        dUp=0;
							      if (dDown>fibLen) dDown=fibLen;
							      if (dDown<0)      dDown=0;

							      double lambda = attenuation_parameters[channel_calib][0];
							      double attUp = exp( -dUp / lambda );
							      double attDown = exp( -dDown / lambda );

							      double US_unatten_E = hClust_points[i]->E_US()*attUp;
							      double DS_unatten_E = hClust_points[i]->E_DS()*attDown;
	
							      double zLocal_clust = m_BCALGeom->GetBCAL_inner_rad()/tan((**lClust).theta()) + m_z_target_center - m_BCALGeom->GetBCAL_center();
							      double dUp_clust = 0.5 * fibLen + zLocal_clust;
							      double dDown_clust = 0.5 * fibLen - zLocal_clust;

							     double attUp_clust = exp( -dUp_clust / lambda );
							     double attDown_clust = exp( -dDown_clust / lambda );

							     double US_reattn_E = US_unatten_E/attUp_clust;
							     double DS_reattn_E = DS_unatten_E/attDown_clust;
							     point_reatten_E = 0.5 * ( US_reattn_E + DS_reattn_E);

							}
						}
					}

	                                if (lClust_points.size() == 1) {

                                                for( unsigned int i = 0 ; i < lClust_points.size() ; i++){

                                                        if (lClust_points[i]->z() > low_z_lim && lClust_points[i]->z() < high_z_lim) point_reatten_E = 0.;
                                                        else{
                                                              int channel_calib = 16*(lClust_points[i]->module()-1)+4*(lClust_points[i]->layer()-1)+lClust_points[i]->sector()-1;

                                                              double fibLen = m_BCALGeom->GetBCAL_length();

                                                              double point_z = lClust_points[i]->z();
                                                              double zLocal = point_z + m_z_target_center - m_BCALGeom->GetBCAL_center();

                                                              double dUp = 0.5 * fibLen + zLocal;
                                                              double dDown = 0.5 * fibLen - zLocal;
                                                              if (dUp>fibLen)   dUp=fibLen;
                                                              if (dUp<0)        dUp=0;
                                                              if (dDown>fibLen) dDown=fibLen;
                                                              if (dDown<0)      dDown=0;

                                                              double lambda = attenuation_parameters[channel_calib][0];
                                                              double attUp = exp( -dUp / lambda );
                                                              double attDown = exp( -dDown / lambda );

                                                              double US_unatten_E = lClust_points[i]->E_US()*attUp;
                                                              double DS_unatten_E = lClust_points[i]->E_DS()*attDown;

                                                              double zLocal_clust = m_BCALGeom->GetBCAL_inner_rad()/tan((**hClust).theta()) + m_z_target_center - m_BCALGeom->GetBCAL_center();
                                                              double dUp_clust = 0.5 * fibLen + zLocal_clust;
                                                              double dDown_clust = 0.5 * fibLen - zLocal_clust;

                                                             double attUp_clust = exp( -dUp_clust / lambda );
                                                             double attDown_clust = exp( -dDown_clust / lambda );

                                                             double US_reattn_E = US_unatten_E/attUp_clust;
                                                             double DS_reattn_E = DS_unatten_E/attDown_clust;
                                                             point_reatten_E = 0.5 * ( US_reattn_E + DS_reattn_E);

                                                        }
                                                }
                                       } 
                                        
						if( (**lClust).Q() == 1 && (**hClust).Q() == 0) {
	                                                (**lClust).mergeClust(**hClust, point_reatten_E);
        	                                        delete *hClust;
                	                                clusters.erase( hClust );
                        	                }
        
						else {
							(**hClust).mergeClust(**lClust, point_reatten_E);
                                                	delete *lClust;
                                                	clusters.erase( lClust );
						}
                                        
					// now iterators are invalid and we need to bail out of loops
					stillMerging = true;
					break;
				}
			}
			if( stillMerging ) break;
		}
	}
}

bool
DBCALCluster_factory::overlap( const DBCALCluster& highEClust,
		const DBCALCluster& lowEClust ) const {

	float sigTheta = fabs( highEClust.theta() - lowEClust.theta() ) / 
		sqrt( highEClust.sigTheta() * highEClust.sigTheta() +
				lowEClust.sigTheta()  * lowEClust.sigTheta() );

	// difference in phi is tricky due to overlap at 0/2pi
	// order based on phi and then take the minimum of the difference
	// and the difference with 2pi added to the smallest

	float deltaPhi = highEClust.phi() - lowEClust.phi();
	float deltaPhiAlt = ( highEClust.phi() > lowEClust.phi() ? 
			highEClust.phi() - lowEClust.phi() - 2*TMath::Pi() :
			lowEClust.phi() - highEClust.phi() - 2*TMath::Pi() );

	deltaPhi = min( fabs( deltaPhi ), fabs( deltaPhiAlt ) );

	float sigPhi = deltaPhi / 
		sqrt( highEClust.sigPhi() * highEClust.sigPhi() +
				lowEClust.sigPhi()  * lowEClust.sigPhi() );

	//We can't rely entirely on sigTheta and sigPhi as defined above.
	//For high-energy clusters, the position uncertainties will be very small,
	//so sigTheta/sigPhi will be large, and clusters may not merge properly.
	//To fix this, force clusters to merge if delta_z and delta_phi are both
	//very small. This is hopefully only a temporary fix.

	//deltaPhi_force_merge and delta_z_force_merge were determined by looking
	//at the separation of decay photons from pi0's from a pythia sample.
	//There are no events where the decay photons have separation
	//(delta_phi < 0.2 && delta_z < 25 cm), so in most cases it should be safe
	//to merge clusters together if they are so close.
	const double deltaPhi_force_merge = 0.1; //radians
	const double delta_z_force_merge = 15.0*k_cm;

	//A major cause of extra clusters are lower energy hits, which have poor
	//z-resolution and so are not properly merged. Treat low energy
	//clusters (< 40 MeV) as a special case. Again, hopefully this is only
	//a temporary fix until we have a more comprehensive solution.
	const double delta_z_force_merge_low_E = 40.0*k_cm;
	const double low_E = .04*k_GeV;

	double z1 = m_BCALGeom->GetBCAL_inner_rad()/tan(highEClust.theta());
	double z2 = m_BCALGeom->GetBCAL_inner_rad()/tan(lowEClust.theta());
	double delta_z = fabs(z1-z2);

	bool theta_match = (sigTheta < m_mergeSig) || (delta_z < delta_z_force_merge) || (delta_z < delta_z_force_merge_low_E && lowEClust.E() < low_E);

	bool phi_match = (sigPhi < m_mergeSig) || (deltaPhi < deltaPhi_force_merge);

	//very loose cut to check that the two clusters are in time
	bool time_match = (highEClust.t() - lowEClust.t()) < m_timeCut;

	if(BCALCLUSTERVERBOSE>1) cout << " clust merge: " << " theta match success = " << theta_match << " phi match = " << phi_match << " time match = " << time_match << " high E = " << highEClust.E() << " low E = " << lowEClust.E() << " highE z = " << z1 << " lowE z = " << z2 << " deltaTheta = " << fabs(highEClust.theta()-lowEClust.theta()) << " sigTheta = " << sigTheta << " highE sigTheta = " << highEClust.sigTheta() << " lowE sigTheta = " << lowEClust.sigTheta() << endl;

	vector<const DBCALPoint*> highE_points;
        highE_points = (highEClust).points();

	vector<const DBCALPoint*> lowE_points;
	lowE_points = (lowEClust).points();

	double highE_summed_z = 0.;
        double highE_summed_phi = 0.;
        double highE_summed_zphi = 0.;
        double highE_summed_z_sq = 0.;
        double highE_slope = 0.;
        double highE_y_intercept = 0.;

	double lowE_summed_z = 0.;
        double lowE_summed_phi = 0.;
        double lowE_summed_zphi = 0.;
        double lowE_summed_z_sq = 0.;
        double lowE_slope = 0.;
        double lowE_y_intercept = 0.;

	int connected = 0;
//	double z_match = 50.;
	double slope_match = 0.01;
	double intercept_match = 1.8;
	double deltaPhi_match = 0.2;

	int lowE_global_sector = 0;
	int highE_global_sector = 0;
	int lowE_point_layer = 0;

        for(unsigned int i = 0 ; i < lowE_points.size() ; i ++){
		// adjust the points phi position to be close to the cluster phi position at the 0/2pi phi boundary
		if(lowEClust.phi() > lowE_points[i]->phi() ){
	   		if( fabs( lowEClust.phi() - lowE_points[i]->phi() - 2*TMath::Pi() ) < TMath::Pi() ) lowE_points[i]->add2Pi();
	  	}
	  	else{
   	 		if( fabs( lowE_points[i]->phi() - lowEClust.phi() - 2*TMath::Pi() ) < TMath::Pi() ) lowE_points[i]->sub2Pi();
	 	 }

          	// compute quantities to be used to calculate the direction of the lower energy cluster if we need it for merging. 
	        lowE_summed_z += lowE_points[i]->z();
                lowE_summed_phi += lowE_points[i]->phi();;
                lowE_summed_zphi += lowE_points[i]->z()*lowE_points[i]->phi();
                lowE_summed_z_sq += lowE_points[i]->z()*lowE_points[i]->z();
                if(lowE_points.size()==1) {
			lowE_global_sector = 4*(lowE_points[i]->module()-1) + lowE_points[i]->sector();
			lowE_point_layer = lowE_points[i]->layer();  
  		}
	 }
       
	 for(unsigned int i = 0 ; i < highE_points.size() ; i ++){
        	// adjust the points phi position to be close to the cluster phi position at the 0/2pi phi boundary	  
	        if(highEClust.phi() > highE_points[i]->phi() ){
                        if( fabs( highEClust.phi() - highE_points[i]->phi() - 2*TMath::Pi() ) < TMath::Pi() ) highE_points[i]->add2Pi();
                }
                else{
                        if( fabs( highE_points[i]->phi() - highEClust.phi() - 2*TMath::Pi() ) < TMath::Pi() ) highE_points[i]->sub2Pi();
                }
		// compute quantities to be used to calculate the direction of the higher energy cluster if we need it for merging.
		highE_summed_z += highE_points[i]->z();
                highE_summed_phi += highE_points[i]->phi();;
                highE_summed_zphi += highE_points[i]->z()*highE_points[i]->phi();
                highE_summed_z_sq += highE_points[i]->z()*highE_points[i]->z();
        	highE_global_sector = 4*(highE_points[i]->module()-1) + highE_points[i]->sector();
		if(lowE_points.size()==1 && lowE_point_layer == highE_points[i]->layer() && ( lowE_global_sector+1 == highE_global_sector || lowE_global_sector-1 == highE_global_sector ) ) connected = 1; // clustesr that contain only a single point won't have any fit parameters and will make it hard for them to merge, this connected int will force a merge if a single point cluster is connected to a cluster without any points adjacent to it.
	}

	// calculate slopes and intercepts of the 2 clusters direction and if one of the clusters is matched to a track then we will require their fit parameter quantities
	// to match for their merging criteria. This allows us to relax their phi position proximity merging criteria since split clusters matched to tracks tend to be
	// further distributed in the azimuthal direction than neutral clusters. 

	highE_slope = (highE_summed_z*highE_summed_phi - highE_points.size()*highE_summed_zphi)/(highE_summed_z*highE_summed_z - highE_points.size()*highE_summed_z_sq);
        highE_y_intercept = (highE_summed_zphi*highE_summed_z - highE_summed_phi*highE_summed_z_sq)/(highE_summed_z*highE_summed_z - highE_points.size()*highE_summed_z_sq);	

	lowE_slope = (lowE_summed_z*lowE_summed_phi - lowE_points.size()*lowE_summed_zphi)/(lowE_summed_z*lowE_summed_z - lowE_points.size()*lowE_summed_z_sq);
        lowE_y_intercept = (lowE_summed_zphi*lowE_summed_z - lowE_summed_phi*lowE_summed_z_sq)/(lowE_summed_z*lowE_summed_z - lowE_points.size()*lowE_summed_z_sq);

	double delta_slope = fabs(highE_slope - lowE_slope) ;
	double delta_intercept = fabs(highE_y_intercept - lowE_y_intercept) ;

	highE_points.clear();
	lowE_points.clear();


	// If both clusters trying to merge together were NOT matched to a track then use neutral clusterizer merging critera of theta and phi matching.
	// If EITHER of the 2 clusters trying to merge together were amtched to a track then use the information about the direction of the cluster for merging.

	if (highEClust.Q() == 0 && lowEClust.Q() == 0 ) return theta_match && phi_match && time_match;

	else return ( ( delta_slope < slope_match && delta_intercept < intercept_match && deltaPhi < deltaPhi_match ) || connected == 1 ) ;


}

bool
DBCALCluster_factory::overlap( const DBCALCluster& clust,
		const DBCALPoint* point ) const {


	float deltaTheta = fabs( clust.theta() - point->theta() );
	/* sigTheta not used
	   float sigTheta = deltaTheta / sqrt( clust.sigTheta() * clust.sigTheta() +
	   point->sigTheta()  * point->sigTheta() );
	   */

	// difference in phi is tricky due to overlap at 0/2pi
	// order based on phi and then take the minimum of the difference
	// and the difference with 2pi added to the smallest

	float deltaPhi = clust.phi() - point->phi();
	float deltaPhiAlt = ( clust.phi() > point->phi() ? 
			clust.phi() - point->phi() - 2*TMath::Pi() :
			point->phi() - clust.phi() - 2*TMath::Pi() );

	deltaPhi = min( fabs( deltaPhi ), fabs( deltaPhiAlt ) );

	/* sigPhi not used
	   float sigPhi = deltaPhi / 
	   sqrt( clust.sigPhi() * clust.sigPhi() +
	   point->sigPhi()  * point->sigPhi() );
	   */

	float rho = ( clust.rho() + point->rho() ) / 2;
	float theta = ( clust.theta() + point->theta() ) / 2;

	float sep = sqrt( ( rho * deltaTheta ) * ( rho * deltaTheta ) +
			( rho * sin( theta ) * deltaPhi ) * ( rho * sin( theta ) * deltaPhi ) );

	float sep_term1 = rho*deltaTheta;
	float sep_term2 = rho*sin(theta)*deltaPhi;

	//very loose cuts to make sure the two hits are in time
	bool time_match = fabs(clust.t() - point->t()) < m_timeCut;

	double clust_z = clust.rho()*cos(clust.theta());

	//double c1 = C1_parm->Eval(clust_z);
	double c1=23.389+19.093*tanh(-0.0104*(clust_z-201.722));

	//double c2 = C2_parm->Eval(clust_z);
	double c2=0.151+0.149*tanh(-0.016*(clust_z-275.194));

	//dtheta_inclusion_curve->SetParameter(0,c1);
	//dtheta_inclusion_curve->SetParameter(1,c2); 
	
	//double inclusion_val = sep_inclusion_curve->Eval(sep);
	double inclusion_val=exp(-sep/30.)-0.1;

        //double inclusion_val1 = dtheta_inclusion_curve->Eval(sep_term1);
	double inclusion_val1=exp(-(sep_term1-0.1)/c1)-c2+.15;
	
        //double inclusion_val2 = dphi_inclusion_curve->Eval(sep_term2);	
	double inclusion_val2=exp(-(sep_term2-2.)/2.5)-sep_term2*0.002+0.07;
	
	// We consider fractional energy (point.E/Clust.E) as a function of spatial separation between
	// a point and cluster to determine if the point should be included in the cluster.
	// These distributions are tighter in the phihat direction than along thetahat. For more details
	// on how the selection criteria for cluster,point overlap function go to logbook entry 3396018.	

	if(BCALCLUSTERVERBOSE>0) cout << "(m,l,s) = (" <<point->module()<<","<<point->layer()<<","<<point->sector()<<")" <<  " sep = " << sep << "sep1 = " << sep_term1 << " sep2 = " << sep_term2 << " inclusion value = " << inclusion_val << " inclusion val1= " << inclusion_val1 << " inclusion val2= " << inclusion_val2<< " time match = " << time_match << " clust E = " << clust.E() << " point E = " << point->E() << " energy ratio = " << point->E()/(point->E()+clust.E()) <<  " clust theta = " << clust.theta()*180./3.14159 << " point theta = " << point->theta()*180./3.14159 << " sep rho*deltaTheta = " << ( rho * deltaTheta ) << endl;

	if(sep>m_moliereRadius && sep<7.*m_moliereRadius &&sep_term2>=2.*m_moliereRadius){
                return ((point->E()/(point->E()+clust.E())) < inclusion_val1 ) && ((point->E()/(point->E()+clust.E())) < inclusion_val2 ) && time_match && deltaPhi*180./3.14159<10.;
        }

	else{
		return ((point->E()/(point->E()+clust.E())) < (inclusion_val1+.2)) && sep < 10.*m_moliereRadius && time_match && sep_term2<2.*m_moliereRadius;
	}

}


bool
DBCALCluster_factory::overlap_charged( const DBCALCluster& clust,
		const DBCALPoint* point, float tracked_phi) const {


	// difference in phi is tricky due to overlap at 0/2pi
	// order based on phi and then take the minimum of the difference
	// and the difference with 2pi added to the smallest

	float phiCut = 0.65417;

	vector<const DBCALPoint*> assoc_points;
	assoc_points = (clust).points();

	double summed_r = 0.;
	double summed_phi = 0.;
	double summed_rphi = 0.;
	double summed_r_sq = 0.;

	double summed_z = 0.;
	double summed_zphi = 0.;
	double summed_z_sq = 0.;
	
	double slope = 0.;
	double y_intercept = 0.;

	int point_global_sector = 4*(point->module()-1) + point->sector();
	int point_layer = point->layer();
	int connected = 0;

	for(unsigned int i = 0 ; i < assoc_points.size() ; i ++){
		int assoc_point_global_sector = 4*(assoc_points[i]->module() - 1) + assoc_points[i]->sector();
		if( point_layer == assoc_points[i]->layer() && ( point_global_sector + 1 == assoc_point_global_sector || point_global_sector - 1 == assoc_point_global_sector) ) connected = 1;
		summed_r += assoc_points[i]->r();
		summed_z += assoc_points[i]->z();
		if( tracked_phi > assoc_points[i]->phi() ){
                        if( fabs( tracked_phi - assoc_points[i]->phi() - 2*TMath::Pi() ) < TMath::Pi() ) assoc_points[i]->add2Pi();
                }
                else{
                        if( fabs( assoc_points[i]->phi() - tracked_phi - 2*TMath::Pi() ) < TMath::Pi() ) assoc_points[i]->sub2Pi();
                 }
	
		summed_phi += assoc_points[i]->phi();
                summed_rphi += assoc_points[i]->r()*assoc_points[i]->phi();
                summed_r_sq += assoc_points[i]->r()*assoc_points[i]->r();
                summed_zphi += assoc_points[i]->z()*assoc_points[i]->phi();
                summed_z_sq += assoc_points[i]->z()*assoc_points[i]->z();

	}

	if(assoc_points.size()<2){
		slope = (tracked_phi - summed_phi)/(m_BCALGeom->GetBCAL_inner_rad() - summed_r);
		y_intercept = tracked_phi - slope*m_BCALGeom->GetBCAL_inner_rad();
	}

        else{
		slope = (summed_z*summed_phi - assoc_points.size()*summed_zphi)/(summed_z*summed_z - assoc_points.size()*summed_z_sq);
                y_intercept = (summed_zphi*summed_z - summed_phi*summed_z_sq)/(summed_z*summed_z - assoc_points.size()*summed_z_sq);
	}

	float fit_phi = 0.;

	if(assoc_points.size() < 2) fit_phi = slope*point->r() + y_intercept;
        else fit_phi = slope*point->z() + y_intercept;

	assoc_points.clear();

	float deltaPhi = fit_phi-point->phi();
	float deltaPhiAlt = ( fit_phi  > point->phi() ? 
                        fit_phi  - point->phi() - 2*TMath::Pi() :
                        point->phi() - fit_phi - 2*TMath::Pi() );

	deltaPhi = min( fabs( deltaPhi ), fabs( deltaPhiAlt ) );

	float rho = point->rho();
	float theta = point->theta();

	float deltaTheta = fabs( clust.theta() - point->theta() );

	float sep = sqrt( ( rho * deltaTheta ) * ( rho * deltaTheta ) +
			( rho * sin( theta ) * deltaPhi ) * ( rho * sin( theta ) * deltaPhi ) );

	float sep_term1 = rho*deltaTheta;
	float sep_term2 = rho*sin(theta)*deltaPhi;

	//very loose cuts to make sure the two hits are in time
	bool time_match = fabs(clust.t() - point->t()) < m_timeCut;

        bool phi_match = fabs( clust.phi() - point->phi() ) < phiCut;

	double clust_z = clust.rho()*cos(clust.theta());

	//double c1 = C1_parm->Eval(clust_z);
	double c1=23.389+19.093*tanh(-0.0104*(clust_z-201.722));

	//double c2 = C2_parm->Eval(clust_z);
	double c2=0.151+0.149*tanh(-0.016*(clust_z-275.194));

	//dtheta_inclusion_curve->SetParameter(0,c1);
	//dtheta_inclusion_curve->SetParameter(1,c2); 
	
	//double inclusion_val = sep_inclusion_curve->Eval(sep);
	double inclusion_val=exp(-sep/30.)-0.1;

        //double inclusion_val1 = dtheta_inclusion_curve->Eval(sep_term1);
	double inclusion_val1=exp(-(sep_term1-0.1)/c1)-c2+.15;

	double inclusion_val2 = exp(-(sep_term2-2.)/1.5) - sep_term2*.007 + .15;
	
	// We consider fractional energy (point.E/Clust.E) as a function of spatial separation between
	// a point and cluster to determine if the point should be included in the cluster.
	// These distributions are tighter in the phihat direction than along thetahat. For more details
	// on how the selection criteria for cluster,point overlap function go to logbook entry 3396018.	

	if(BCALCLUSTERVERBOSE>1) cout << "(m,l,s) = (" <<point->module()<<","<<point->layer()<<","<<point->sector()<<")" <<  " sep = " << sep << "sep1 = " << sep_term1 << " sep2 = " << sep_term2 << " inclusion value = " << inclusion_val << " inclusion val1= " << inclusion_val1 << " inclusion val2= " << inclusion_val2<< " time match = " << time_match << " clust E = " << clust.E() << " point E = " << point->E() << " energy ratio = " << point->E()/(point->E()+clust.E()) <<  " clust theta = " << clust.theta()*180./3.14159 << " point theta = " << point->theta()*180./3.14159 << " sep rho*deltaTheta = " << ( rho * deltaTheta ) << endl;

	if(sep>m_moliereRadius && sep<7.*m_moliereRadius &&sep_term2>=2.*m_moliereRadius){
                return ((point->E()/(point->E()+clust.E())) < (inclusion_val1) ) && ((point->E()/(point->E()+clust.E())) < (inclusion_val2) ) && time_match && phi_match;
        }

        else{
                return ((point->E()/(point->E()+clust.E())) < (inclusion_val1 + .2)) && sep < 10.*m_moliereRadius && time_match && sep_term2<2.*m_moliereRadius;
        }

	return connected == 1;

}


bool
DBCALCluster_factory::overlap( const DBCALCluster& clust,
		const DBCALUnifiedHit* hit ) const {

	int cellId = m_BCALGeom->cellId( hit->module, hit->layer, hit->sector );

	float cellPhi = m_BCALGeom->phi( cellId );
	float cellSigPhi = m_BCALGeom->phiSize( cellId );

	// annoying +- 2pi business to try to find the best delta phi

	float deltaPhi = clust.phi() - cellPhi;
	float deltaPhiAlt = ( clust.phi() > cellPhi ? 
			clust.phi() - cellPhi - 2*TMath::Pi() :
			cellPhi - clust.phi() - 2*TMath::Pi() );
	deltaPhi = min( fabs( deltaPhi ), fabs( deltaPhiAlt ) );  

	float sigPhi = deltaPhi / 
		sqrt( clust.sigPhi() * clust.sigPhi() + cellSigPhi  * cellSigPhi );

	int channel_calib = 16*(hit->module-1)+4*(hit->layer-1)+hit->sector-1; // need to use cellID for objects in DBCALGeometry but the CCDB uses a different channel numbering scheme, so use channel_calib when accessing CCDB tables.
	// given the location of the cluster, we need the best guess
	// for z with respect to target at this radius
	double z = clust.rho()*cos(clust.theta()) + m_z_target_center;        
	double d = ( (hit->end == 0) ? (z - m_BCALGeom->GetBCAL_center() + m_BCALGeom->GetBCAL_length()/2.0) : (m_BCALGeom->GetBCAL_center() + m_BCALGeom->GetBCAL_length()/2.0 - z));  // d gives the distance to upstream or downstream end of BCAL depending on where the hit was with respect to the cluster z position.
	double time_corr = hit->t - d/effective_velocities[channel_calib];  // hit time corrected to the interaction point in the bar.        
	double time_diff = TMath::Abs(clust.t() - time_corr); // time cut between cluster time and hit time - 20 ns is a very loose time cut.

	return( sigPhi < m_mergeSig && time_diff < m_clust_hit_timecut ); 

}
