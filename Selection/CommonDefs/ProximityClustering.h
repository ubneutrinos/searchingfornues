/**
 * \file ProximityClusterer.h
 *
 * \ingroup Clusterer
 * 
 * \brief Class def header for a class ProximityClusterer
 *
 * @author david caratelli
 */

/** \addtogroup Clusterer

    @{*/

#ifndef PROXIMITYCLUSTERING_H
#define PROXIMITYCLUSTERING_H

#include <map>

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

namespace searchingfornues {
  /**
     \class ProximityClustering
     User custom analysis class made by david caratelli
   */
  /// Map making function
  void MakeHitMap(const std::vector< art::Ptr<recob::Hit> >& hitlist,
		  int plane,
		  const float& _time2cm, const float& _wire2cm,
		  const float& _cellSize,
		  std::map<std::pair<int,int>, std::vector<size_t> >& _hitMap) {

    _hitMap.clear();
    // temporary pair
    std::pair<int,int> tmpPair;
    
    
    for (size_t h=0; h < hitlist.size(); h++){
      
      auto const& hit = hitlist.at(h);
      // skip if not of plane we want
      if (hit->View() != plane)
	continue;

      double t = hit->PeakTime()*_time2cm;
      double w = hit->WireID().Wire*_wire2cm;

      // map is (i,j) -> hit list
      // i : ith bin in wire of some width
      // j : jth bin in time of some width
      int i = int(w/_cellSize);
      int j = int(t/_cellSize);
      tmpPair = std::make_pair(i,j);
      // does this entry exist in the map?
      // if yes -> append to vector
      // if no create new vector and add to map
      if (_hitMap.find(tmpPair) == _hitMap.end()){
	std::vector<size_t> aaa = {h};
	_hitMap[tmpPair] = aaa;
      }
      else
	_hitMap[tmpPair].push_back(h);
    }// for all hits

    return;
  }

  /// check if time overlaps
  bool TimeOverlap(const art::Ptr<recob::Hit>& h1, const art::Ptr<recob::Hit>& h2,
		   const float& _time2cm, double& dmin) {

    auto T1 = h1->PeakTime() * _time2cm; // time of first hit
    auto T2 = h2->PeakTime() * _time2cm;
    auto W1 = h1->RMS() * _time2cm;
    auto W2 = h2->RMS() * _time2cm;
    
    double d = dmin;
    
    if (T1 > T2) {
      
      if ( (T2+W2) > (T1-W1) ) return true;
      
      d = (T1-W1) - (T2+W2);
      if (d < dmin) dmin = d;
      
    }
    
    else {
      
      if ( (T1+W1) > (T2-W2) ) return true;
      
      d = (T2-W2) - (T1+W1);
      if (d < dmin) dmin = d;
      
    }
    
    return false;
  }



  
  /// Functions to decide if two hits should belong to the same cluster or not
  bool HitsCompatible(const art::Ptr<recob::Hit>& h1, const art::Ptr<recob::Hit>& h2,
		      const float& _time2cm, const float& _wire2cm, const float& _radius) {

    if (h1->WireID().Plane != h2->WireID().Plane)
      return false;
    
    double dt = ( h1->PeakTime() - h2->PeakTime() ) * _time2cm;
    //  if the hit time-ranges overlap, this distnce should be 0
    if (TimeOverlap(h1,h2,_time2cm,dt) == true)
      dt = 0;
    double dw = fabs(((double)h1->Channel()-(double)h2->Channel())*_wire2cm);
    if (dw >  0.3) dw -= 0.3;
    //if (dw < -0.3) dw  = 0.3;
    double d = dt*dt + dw*dw;

    if (d > (_radius*_radius))
      return false;

    return true;
  }

  
  /// Function to get neighboring hits (from self + neighoring cells)
  void getNeighboringHits(const std::pair<int,int>& pair, std::vector<size_t>& hitIndices,
			  std::map<std::pair<int,int>, std::vector<size_t> >& _hitMap) {

    auto const& i       = pair.first;
    // time-space cell index
    auto const& j       = pair.second;

    // _________
    // |__|__|__|
    // |__|XX|__|
    // |__|__|__|
    if (_hitMap.find(std::make_pair(i,j)) != _hitMap.end()){
      for (auto &h : _hitMap[std::make_pair(i,j)])
	hitIndices.push_back(h);
    }

    // now look at neighboring cells, if they exist
    // _________
    // |__|__|__|
    // |XX|__|__|
    // |__|__|__|
    if (_hitMap.find(std::make_pair(i-1,j)) != _hitMap.end()){
      for (auto &h : _hitMap[std::make_pair(i-1,j)])
	hitIndices.push_back(h);
    }
    // _________
    // |__|__|__|
    // |__|__|__|
    // |__|XX|__|
    if (_hitMap.find(std::make_pair(i,j-1)) != _hitMap.end()){
      for (auto &h : _hitMap[std::make_pair(i,j-1)])
	hitIndices.push_back(h);
    }
    // _________
    // |__|__|__|
    // |__|__|__|
    // |XX|__|__|
    if ( _hitMap.find(std::make_pair(i-1,j-1)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i-1,j-1)])
	hitIndices.push_back(h);
    }
    // _________
    // |__|XX|__|
    // |__|__|__|
    // |__|__|__|
    if ( _hitMap.find(std::make_pair(i,j+1)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i,j+1)])
	hitIndices.push_back(h);
    }
    // _________
    // |__|__|__|
    // |__|__|XX|
    // |__|__|__|
    if ( _hitMap.find(std::make_pair(i+1,j)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i+1,j)])
	hitIndices.push_back(h);
    }
    // _________
    // |__|__|XX|
    // |__|__|__|
    // |__|__|__|
    if ( _hitMap.find(std::make_pair(i+1,j+1)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i+1,j+1)])
	hitIndices.push_back(h);
    }
    // _________
    // |XX|__|__|
    // |__|__|__|
    // |__|__|__|
    if ( _hitMap.find(std::make_pair(i-1,j+1)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i-1,j+1)])
	hitIndices.push_back(h);
    }
    // _________
    // |__|__|__|
    // |__|__|__|
    // |__|__|XX|
    if ( _hitMap.find(std::make_pair(i+1,j-1)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i+1,j-1)])
	hitIndices.push_back(h);
    }
  }

  
  
  
  /**
       Cluster function:
       @brief cluster hits based on proximity
       Input: pointer to event hit record.
       Output: vector of vector of hit indices which make up clusters
  */
  bool cluster(const std::vector< art::Ptr<recob::Hit> >& hit_ptr_v,
	       std::vector<std::vector<unsigned int> >& _out_cluster_vector,
	       const float& cellSize, const float& radius) {

    if (hit_ptr_v.size() == 0)
      return false;
    
    /// size of each cell [cm]
    double _cellSize = cellSize;
    
    /// radius to count charge around [cm]
    double _radius = radius;
    
    /// conversion factors for hits
    double _wire2cm, _time2cm;
    
    // get detector specific properties
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    _wire2cm = geom->WirePitch(0,0,0);
    _time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );
    
    
    /// map connecting coordinate index (i,j) to [h1,h2,h3] (hit index list)
    std::map<std::pair<int,int>, std::vector<size_t> > _hitMap;

    // a map to connect hit index wih a cluster index
    // each hit gets a cluster index
    // _clusterMap[hit_index] -> cluster_index
    std::map<size_t, size_t> _clusterMap;
    // a map to connect the cluster index with the vector of hit indices for that cluster
    // _clusters[index] -> vector of hit indices for that cluster
    std::map<size_t,std::vector<size_t> > _clusters;

    // keep track of largest cluster ID created
    size_t maxClusterID = 0;

    for (int pl=0; pl < 3; pl++){
      
      // hit map will only contain hits we want to use for clustering
      MakeHitMap(hit_ptr_v,pl,_time2cm, _wire2cm, _cellSize, _hitMap);
      
      // iterator for hit cell map
      std::map<std::pair<int,int>, std::vector<size_t> >::iterator it;
      
      // loop through hits in each cell to find matches
      for (it = _hitMap.begin(); it != _hitMap.end(); it++){
	
	// pair = (i,j) indices of this cell in the _hitMap
	auto const& pair = it->first;
	
	// wire-space cell index
	// prepare a hit list of all neighboring cells
	// _________
	// |__|__|__|
	// |__|__|__|
	// |__|__|__|
	std::vector<size_t> cellhits = it->second;

	std::vector<size_t> neighborhits;
	getNeighboringHits(pair,neighborhits, _hitMap);

	for (size_t h1=0; h1 < cellhits.size(); h1++){

	  // has this hit been added to a cluster?
	  // if so not necessary to look at
	  auto const& hit1 = cellhits[h1];
	  // keep track if the hit will ever be matched to another
	  bool matched = false;
	  // if not find hits it should be clustered with and add it to the appropriate cluster
	  for (size_t h2=0; h2 < neighborhits.size(); h2++){
	    auto const& hit2 = neighborhits[h2];
	    if (hit1 == hit2) continue;
	    // are the hits compatible?
	    bool compat = HitsCompatible(hit_ptr_v.at(hit1),
					 hit_ptr_v.at(hit2),
					 _time2cm, _wire2cm, _radius);
	    // should the hits go in the same cluster?
	    if (compat){
	      matched = true;
	      // if both hits have already been assigned to a cluster then we can merge the cluster indices!
	      if ( (_clusterMap.find(hit1) != _clusterMap.end()) and
		   (_clusterMap.find(hit2) != _clusterMap.end()) ){
		// if in the same cluster -> do nothing
		// if they are in different clusters:
		if (_clusterMap[hit1] != _clusterMap[hit2]){
		  auto idx1 = _clusterMap[hit1];
		  auto idx2 = _clusterMap[hit2];
		  // hit indices for 1st cluster:
		  auto hits1 = _clusters[idx1];
		  auto hits2 = _clusters[idx2];
		  // append hits2 to hits1
		  for (auto h : hits2){
		    hits1.push_back(h);
		    // also change the index that the hit goes to (idx1 instead of idx2)
		    _clusterMap[h] = idx1;
		  }
		  _clusters[idx1] = hits1;
		  // erase cluster @ index2
		  _clusters.erase(idx2);
		}// if they are in different clusters
	      }
	      // if compatible and the 2nd hit has been added to a cluster
	      // add hit1 to the same cluster
	      else if ( (_clusterMap.find(hit2) != _clusterMap.end()) and
			(_clusterMap.find(hit1) == _clusterMap.end()) ){
		auto clusIdx = _clusterMap[hit2];
		_clusterMap[hit1] = clusIdx;
		_clusters[clusIdx].push_back(hit1);
	      }
	      // otherwise, add both to a new cluster
	      else if ( (_clusterMap.find(hit1) != _clusterMap.end()) and
			(_clusterMap.find(hit2) == _clusterMap.end()) ){
		auto clusIdx = _clusterMap[hit1];
		_clusterMap[hit2] = clusIdx;
		_clusters[clusIdx].push_back(hit2);
	      }
	      // if neither has a cluster yet
	      else{
		// create a new cluster for this match
		_clusterMap[hit1] = maxClusterID;
		_clusterMap[hit2] = maxClusterID;
		std::vector<size_t> cl = {hit1,hit2};
		_clusters[maxClusterID] = cl;
		maxClusterID += 1;
	      }
	    }// if the two hits are compatible
	  }// 2nd loop through hits in the cell
	  // has this hit been matched? if not we still need to add it as its own cluster
	  if (matched == false){
	    _clusterMap[hit1] = maxClusterID;
	    _clusters[maxClusterID] = {hit1};
	    maxClusterID += 1;
	  }
	}// 1st loop through hits in the cell
      }// loop through all cells

    }// loop through all planes

    // make a vector for the clusters
    for (auto it = _clusters.begin(); it != _clusters.end(); it++){
      auto indices = it->second;
      // if there are enough indices, make a cluster
      if (indices.size() >= 1){
	std::vector<unsigned int> clus;
	for (auto idx : indices)
	  clus.push_back(idx);
	_out_cluster_vector.push_back(clus);
      }// if there are 2 hits in cluster
    }
    
    return true;
  }

      
}
#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
