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
}

#ifdef BCAL_CLUSTER_DIAGNOSTIC

jerror_t
DBCALCluster_factory::init(void){

  m_rootFile = new TFile( "bcal_clust_diag.root", "RECREATE" );
  
  m_twoEndPtTr = new TTree( "twoEndPtTr", "BCAL two-ended points" );
  m_twoEndPtTr->Branch( "n2EPt", &m_n2EPt, "n2EPt/I" );
  m_twoEndPtTr->Branch( "rhoPt", m_rhoPt, "rhoPt[n2EPt]/F" );
  m_twoEndPtTr->Branch( "phiPt", m_phiPt, "phiPt[n2EPt]/F" );
  m_twoEndPtTr->Branch( "thePt", m_thePt, "thePt[n2EPt]/F" );
  m_twoEndPtTr->Branch( "rhoSPt", m_rhoSPt, "rhoSPt[n2EPt]/F" );
  m_twoEndPtTr->Branch( "phiSPt", m_phiSPt, "phiSPt[n2EPt]/F" );
  m_twoEndPtTr->Branch( "theSPt", m_theSPt, "theSPt[n2EPt]/F" );
  m_twoEndPtTr->Branch( "ePt", m_ePt, "ePt[n2EPt]/F" );
  m_twoEndPtTr->Branch( "tPt", m_tPt, "tPt[n2EPt]/F" );
  m_twoEndPtTr->Branch( "t0Pt", m_t0Pt, "t0Pt[n2EPt]/F" );
  
  m_firstClustTr = new TTree( "firstClustTr", "First-pass clusters" );
  m_firstClustTr->Branch( "nCl", &m_nCl, "nCl/I" );
  m_firstClustTr->Branch( "nPts", m_nPts, "nPts[nCl]/I" );
  m_firstClustTr->Branch( "rhoCl", m_rhoCl, "rhoCl[nCl]/F" );
  m_firstClustTr->Branch( "phiCl", m_phiCl, "phiCl[nCl]/F" );
  m_firstClustTr->Branch( "theCl", m_theCl, "theCl[nCl]/F" );
  m_firstClustTr->Branch( "rhoSCl", m_rhoSCl, "rhoSCl[nCl]/F" );
  m_firstClustTr->Branch( "phiSCl", m_phiSCl, "phiSCl[nCl]/F" );
  m_firstClustTr->Branch( "theSCl", m_theSCl, "theSCl[nCl]/F" );
  m_firstClustTr->Branch( "eCl", m_eCl, "eCl[nCl]/F" );
  m_firstClustTr->Branch( "tCl", m_tCl, "tCl[nCl]/F" );
  
  m_ovrlpTr = new TTree( "ovrlpTr", "Point Cluster Overlap" );
  m_ovrlpTr->Branch( "dPhi", &m_dPhi, "dPhi/F" );
  m_ovrlpTr->Branch( "dThe", &m_dThe, "dThe/F" );
  m_ovrlpTr->Branch( "sep", &m_sep, "sep/F" );
  m_ovrlpTr->Branch( "sigPhi", &m_sigPhi, "sigPhi/F" );
  m_ovrlpTr->Branch( "sigThe", &m_sigThe, "sigThe/F" );
  m_ovrlpTr->Branch( "eClus", &m_eClus, "eClus/F" );
  m_ovrlpTr->Branch( "rhoClus", &m_rhoClus, "rhoClus/F" );
  m_ovrlpTr->Branch( "phiClus", &m_phiClus, "phiClus/F" );
  m_ovrlpTr->Branch( "theClus", &m_theClus, "theClus/F" );
  m_ovrlpTr->Branch( "nClClus", &m_nClClus, "nClClus/I" );
  
  return NOERROR;

}

jerror_t
DBCALCluster_factory::fini( void ){
  
  m_rootFile->cd();
  
  m_twoEndPtTr->Write();
  m_firstClustTr->Write();
  m_ovrlpTr->Write();
  m_rootFile->Write();
  m_rootFile->Close();
  
  return NOERROR;
}

#endif

jerror_t DBCALCluster_factory::brun(JEventLoop *loop, int runnumber) {
  DApplication* app = dynamic_cast<DApplication*>(loop->GetJApplication());
  DGeometry* geom = app->GetDGeometry(runnumber);
  geom->GetTargetZ(m_z_target_center);

  loop->GetCalib("/BCAL/effective_velocities", effective_velocities);

  loop->GetCalib("/BCAL/attenuation_parameters",attenuation_parameters);

  return NOERROR;
}

jerror_t
DBCALCluster_factory::evnt( JEventLoop *loop, int eventnumber ){

  vector< const DBCALPoint* > twoEndPoint;
  loop->Get(twoEndPoint);

#ifdef BCAL_CLUSTER_DIAGNOSTIC
  
  m_rootFile->cd();
  
  m_n2EPt = twoEndPoint.size();
  assert( m_n2EPt <= MAX_POINT );
  
  for( int i = 0; i < m_n2EPt; ++i ){
   
    m_rhoPt[i] = twoEndPoint[i]->rho();
    m_phiPt[i] = twoEndPoint[i]->phi();
    m_thePt[i] = twoEndPoint[i]->theta();
    m_rhoSPt[i] = twoEndPoint[i]->sigRho();
    m_phiSPt[i] = twoEndPoint[i]->sigPhi();
    m_theSPt[i] = twoEndPoint[i]->sigTheta();
    m_ePt[i] = twoEndPoint[i]->E();
    m_tPt[i] = twoEndPoint[i]->t();
    m_t0Pt[i] = twoEndPoint[i]->tInnerRadius();
  }
  
  m_twoEndPtTr->Fill();
  
#endif // BCAL_CLUSTER_DIAGNOSTIC
  
#ifdef BCAL_CLUSTER_DIAGNOSTIC
  
  m_nCl = clusters.size();
  assert( m_nCl <= MAX_CLUST );
  
  for( int i = 0; i < m_nCl; ++i ){
    
    m_nPts[i] = clusters[i]->nCells();
    m_rhoCl[i] = clusters[i]->rho();
    m_phiCl[i] = clusters[i]->phi();
    m_theCl[i] = clusters[i]->theta();
    m_rhoSCl[i] = clusters[i]->sigRho();
    m_phiSCl[i] = clusters[i]->sigPhi();
    m_theSCl[i] = clusters[i]->.sigTheta();
    m_eCl[i] = clusters[i]->E();
    m_tCl[i] = clusters[i]->t();
  }
  
  m_firstClustTr->Fill();
  
#endif // BCAL_CLUSTER_DIAGNOSTIC
  

  // Want to add singled-ended hits to the Clusters. 

  // Looking for hits that are single-ended.

  vector< const DBCALUnifiedHit* > hits;
  loop->Get(hits);

  // first arrange the list of hits so they are grouped by cell
  map< int, vector< const DBCALUnifiedHit* > > cellHitMap;
  for( vector< const DBCALUnifiedHit* >::const_iterator hitPtr = hits.begin();
      hitPtr != hits.end();
      ++hitPtr ){

    const DBCALUnifiedHit& hit = (**hitPtr);
    
    int id = DBCALGeometry::cellId( hit.module, hit.layer, hit.sector );
    
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
 vector<DBCALCluster*> clusters = clusterize( twoEndPoint, single_ended_hits );
  
  // load our vector of clusters into the factory member data
  for( vector<DBCALCluster*>::iterator clust = clusters.begin();
      clust != clusters.end();
      ++clust ){
   
    // put in an energy threshold for clusters
    if( (**clust).E() < 5*k_MeV ) {
      delete *clust;
      continue;
    }

    _data.push_back(*clust);
  }

  return NOERROR;
}

vector<DBCALCluster*>
DBCALCluster_factory::clusterize( vector< const DBCALPoint* > points , vector< const DBCALUnifiedHit* > hits ) const {

  // first sort the points by energy
  sort( points.begin(), points.end(), PointSort );

  vector<DBCALCluster*> clusters(0);
  
  // ahh.. more hard coded numbers that should probably
  // come from a database or something
  float seedThresh = 1*k_GeV;
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

  while( seedThresh > minSeed ) {

    bool usedPoint = false;
 
    for( vector< const DBCALPoint* >::iterator pt = points.begin();
        pt != points.end();
        ++pt ){
       
      // first see if point should be added to an existing
      // cluster
      
      for( vector<DBCALCluster*>::iterator clust = clusters.begin();
           clust != clusters.end();
          ++clust ){
        
        if( overlap( **clust, *pt ) ){
                    
          (**clust).addPoint( *pt );
          points.erase( pt );
          usedPoint = true;
        }
        
        // once we erase a point the iterator is no longer useful
        // and we start the loop over
        if( usedPoint ) break;
      }
      
      if( usedPoint ) break;
        
      // if the point doesn't overlap with a cluster
      // see if it can become a new seed
      if( (**pt).E() > seedThresh && ((**pt).layer() != 4 || (**pt).E() > layer4_minSeed) ){

        clusters.push_back(new DBCALCluster( *pt, m_z_target_center ) );
        points.erase( pt );
        usedPoint = true;
      }
      
      if( usedPoint ) break;
    }

    merge( clusters );

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
        double d = ( ((**ht).end == 0) ? (z  - DBCALGeometry::GetBCAL_center() + DBCALGeometry::GetBCAL_length()/2.0) : (DBCALGeometry::GetBCAL_center() + DBCALGeometry::GetBCAL_length()/2.0 - z));  // d gives the distance to upstream or downstream end of BCAL depending on where the hit was with respect to the cluster z position.
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
DBCALCluster_factory::merge( vector<DBCALCluster*>& clusters ) const {
  
  if( clusters.size() <= 1 ) return;
  
  sort( clusters.begin(), clusters.end(), ClusterSort );
  
  bool stillMerging = true;
  
  while( stillMerging ){
  
    stillMerging = false;
    for( vector<DBCALCluster*>::iterator hClust = clusters.begin();
        hClust != clusters.end() - 1;
        ++hClust ){
    
      for( vector<DBCALCluster*>::iterator lClust = hClust + 1;
          lClust != clusters.end();
          ++lClust ){
      
        if( overlap( **hClust, **lClust ) ){
                  
          (**hClust).mergeClust(**lClust);
          delete *lClust;
          clusters.erase( lClust );
          
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
                       highEClust.phi() - lowEClust.phi() - 2*PI :
                       lowEClust.phi() - highEClust.phi() - 2*PI );

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
  
  double z1 = DBCALGeometry::GetBCAL_inner_rad()/tan(highEClust.theta());
  double z2 = DBCALGeometry::GetBCAL_inner_rad()/tan(lowEClust.theta());
  double delta_z = fabs(z1-z2);

  bool theta_match = (sigTheta < m_mergeSig) || (delta_z < delta_z_force_merge) || (delta_z < delta_z_force_merge_low_E && lowEClust.E() < low_E);

  bool phi_match = (sigPhi < m_mergeSig) || (deltaPhi < deltaPhi_force_merge);

  //very loose cut to check that the two clusters are in time
  bool time_match = (highEClust.t() - lowEClust.t()) < m_timeCut;

  return theta_match && phi_match && time_match;

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
                        clust.phi() - point->phi() - 2*PI :
                        point->phi() - clust.phi() - 2*PI );

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
  
#ifdef BCAL_CLUSTER_DIAGNOSTIC
  
  m_dPhi = deltaPhi;
  m_dThe = fabs( clust.theta() - point->theta() );
  m_sep  = sep;
  m_sigPhi = sigPhi;
  m_sigThe = sigTheta;
  m_eClus = clust.E();
  m_rhoClus = clust.rho();
  m_phiClus = clust.phi();
  m_theClus = clust.theta();
  m_nClClus = clust.nCells();
  
  m_ovrlpTr->Fill();
  
#endif // BCAL_CLUSTER_DIAGNOSTIC

  //very loose cuts to make sure the two hits are in time
  bool time_match = fabs(clust.t() - point->t()) < m_timeCut;

  if( point->E() / clust.E() < 0.1 ){
    return (sep < 5*m_moliereRadius) && time_match;
  }
  else{
    return (sep < 2*m_moliereRadius) && time_match;
  }

}

bool
DBCALCluster_factory::overlap( const DBCALCluster& clust,
                               const DBCALUnifiedHit* hit ) const {
             
  int cellId = DBCALGeometry::cellId( hit->module, hit->layer, hit->sector );

  float cellPhi = DBCALGeometry::phi( cellId );
  float cellSigPhi = DBCALGeometry::phiSize( cellId );
             
  // annoying +- 2pi business to try to find the best delta phi
               
  float deltaPhi = clust.phi() - cellPhi;
  float deltaPhiAlt = ( clust.phi() > cellPhi ? 
                        clust.phi() - cellPhi - 2*PI :
                        cellPhi - clust.phi() - 2*PI );
  deltaPhi = min( fabs( deltaPhi ), fabs( deltaPhiAlt ) );  
             
  float sigPhi = deltaPhi / 
       sqrt( clust.sigPhi() * clust.sigPhi() + cellSigPhi  * cellSigPhi );
             
  int channel_calib = 16*(hit->module-1)+4*(hit->layer-1)+hit->sector-1; // need to use cellID for objects in DBCALGeometry but the CCDB uses a different channel numbering scheme, so use channel_calib when accessing CCDB tables.
  // given the location of the cluster, we need the best guess
  // for z with respect to target at this radius
  double z = clust.rho()*cos(clust.theta()) + m_z_target_center;        
  double d = ( (hit->end == 0) ? (z - DBCALGeometry::GetBCAL_center() + DBCALGeometry::GetBCAL_length()/2.0) : (DBCALGeometry::GetBCAL_center() + DBCALGeometry::GetBCAL_length()/2.0 - z));  // d gives the distance to upstream or downstream end of BCAL depending on where the hit was with respect to the cluster z position.
  double time_corr = hit->t - d/effective_velocities[channel_calib];  // hit time corrected to the interaction point in the bar.        
  double time_diff = TMath::Abs(clust.t() - time_corr); // time cut between cluster time and hit time - 20 ns is a very loose time cut.
  
  return( sigPhi < m_mergeSig && time_diff < m_clust_hit_timecut ); 

}
