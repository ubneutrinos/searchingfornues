#ifndef PFPHITDISTANCE_H
#define PFPHITDISTANCE_H

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

namespace searchingfornues
{

  /**
   * @brief given two PFP return 2D min distance between hits on the three planes
   * @input pfp_pxy1 : proxy for first pfparticle
   * @input pfp_pxy2 : proxy for second pfparticle
   * @input hitcoll  : proxy connecting clusters to hits
   * @return vector of distances on the three planes [U,V,Y]
   */
  std::vector<float> GetTrackShowerScore(const ProxyPfpElem_t &pfp_pxy1,
					 const ProxyPfpElem_t &pfp_pxy2,
					 const ProxyClusColl_t &hitcoll)
    {

    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    float w2cm = geom->WirePitch(0,0,0);
    float t2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );

      std::vector<float> dist_v = {-1,-1,-1};

      auto clus_pxy1_v = pfp_pxy1.get<recob::Cluster>();
      auto clus_pxy2_v = pfp_pxy2.get<recob::Cluster>();

      // store hits for each plane
      std::vector< std::vector< art::Ptr<recob::Hit> > > cluster1_hits_v;
      std::vector< std::vector< art::Ptr<recob::Hit> > > cluster2_hits_v;
      
      cluster1_hits_v.resize(3);
      cluster2_hits_v.resize(3);

      for (auto ass_clus : clus_pxy1_v) {
	// get cluster proxy
	const auto &clus = hitcoll[ass_clus.key()];
	auto clus_hit_v = clus.get<recob::Hit>();
	auto plane = clus->Plane().Plane;
	if ( (plane >=0) && (plane < 3) ) {
	  cluster1_hits_v[plane].clear();
	  for (size_t h=0; h < clus_hit_v.size(); h++) {
	    cluster1_hits_v[plane].push_back( clus_hit_v[h] );
	  }// for all hits in cluster
	}// if plane is ok
      }// for all clusters for PFP 1

      for (auto ass_clus : clus_pxy2_v) {
	// get cluster proxy
	const auto &clus = hitcoll[ass_clus.key()];
	auto clus_hit_v = clus.get<recob::Hit>();
	auto plane = clus->Plane().Plane;
	if ( (plane >=0) && (plane < 3) ) {
	  cluster2_hits_v[plane].clear();
	  for (size_t h=0; h < clus_hit_v.size(); h++) {
	    cluster2_hits_v[plane].push_back( clus_hit_v[h] );
	  }// for all hits in cluster
	}// if plane is ok
      }// for all clusters for PFP 1

      for (size_t plane=0; plane < 3; plane++) {

	auto hits1 = cluster1_hits_v[plane];
	auto hits2 = cluster2_hits_v[plane];

	if ( hits1.size() && hits2.size() ) {

	  float dmin = 1e6;

	  for (size_t h1 = 0; h1 < hits1.size(); h1++) {
	    for (size_t h2 = 0; h2 < hits2.size(); h2++) {

	      auto hit1 = hits1[h1];
	      auto hit2 = hits2[h2];

	      float dt = (hit1->PeakTime()    - hit2->PeakTime()) * t2cm;
	      float dw = (hit1->WireID().Wire - hit2->WireID().Wire) * w2cm;
	      
	      float dhit = sqrt ( (dt*dt) + (dw*dw) );

	      if (dhit < dmin) { dmin = dhit; }
	      
	    }// second hit loop
	  }// first hit loop

	  dist_v[plane] = dmin;
	  
	}// if both clusters have hits
	
      }// for all three planes

      return dist_v;
    }
  
} // namespace searchingfornues

#endif
