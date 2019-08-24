#ifndef PFPHITDISTANCE_H
#define PFPHITDISTANCE_H

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "ubreco/ShowerReco/ProximityClustering/Algorithms/ProximityClusterer.h"

namespace searchingfornues
{

  /**
   * @brief given two PFP return 2D min distance between hits on the three planes
   * @input pfp_pxy1 : proxy for first pfparticle
   * @input pfp_pxy2 : proxy for second pfparticle
   * @input hitcoll  : proxy connecting clusters to hits
   * @return vector of distances on the three planes [U,V,Y]
   */
  int GetPFParticleSegments(const ProxyPfpElem_t &pfp_pxy,
			    const ProxyClusColl_t &hitcoll) {
    
    // load proximity clustering algorithm
    ::gammacatcher::ProximityClusterer* _ProximityClusterer = new gammacatcher::ProximityClusterer();
    _ProximityClusterer->initialize();
    _ProximityClusterer->setRadius(2.0);
    _ProximityClusterer->setCellSize(2.0);
    
    return 1;
    
  }
  
} // namespace searchingfornues

#endif
