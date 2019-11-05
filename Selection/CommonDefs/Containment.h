#ifndef TRUTHCONTAINMENT_H
#define CONTAINMENT_H

// services for detector properties
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larcore/Geometry/Geometry.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCStep.h"

#include "lardata/Utilities/FindManyInChainP.h"

namespace searchingfornues
{

  /**
   * @brief Determine if the specified point is in the fiducial volume
   *        Not recommended, no array size checking is done.
   *
   * @param x array of 3D location
   * @return True if the point is inside the fiducial volume
   */  
  bool isFiducial(const double x[3], 
		  const double fFidvolXstart, const double fFidvolYstart, const double fFidvolZstart,
		  const double fFidvolXend, const double fFidvolYend, const double fFidvolZend)
  {
    
    art::ServiceHandle<geo::Geometry> geo;
    geo::TPCGeo const &thisTPC = geo->TPC();
    geo::BoxBoundedGeo theTpcGeo = thisTPC.ActiveBoundingBox();
    std::vector<double> bnd = {theTpcGeo.MinX(), theTpcGeo.MaxX(), theTpcGeo.MinY(), theTpcGeo.MaxY(), theTpcGeo.MinZ(), theTpcGeo.MaxZ()};
    bool is_x =
      x[0] > (bnd[0] + fFidvolXstart) && x[0] < (bnd[1] - fFidvolXend);
    bool is_y =
      x[1] > (bnd[2] + fFidvolYstart) && x[1] < (bnd[3] - fFidvolYend);
    bool is_z =
      x[2] > (bnd[4] + fFidvolZstart) && x[2] < (bnd[5] - fFidvolZend);
    
    return is_x && is_y && is_z;
  }


  bool TruthContained(const float& FVxS, const float& FVyS, const float& FVzS,
		      const float& FVxE, const float& FVyE, const float& FVzE,
		      const std::vector<sim::MCShower> &inputMCShower,
		      const std::vector<sim::MCTrack> &inputMCTrack ) {
    
    // require truth-containment by
    // (1) requiring the vertex is in the FV
    // (2) require all MCTracks are contained within the FV
    // (3) require all MCShowers to deposit > some fraction of energy in the FV
    
    bool contained = true;
    
    for (auto mcs : inputMCShower) {
      if (mcs.Process() == "primary" || (mcs.MotherPdgCode() == 111 && mcs.Process() == "Decay" && mcs.MotherProcess() == "primary") ) {
	if (mcs.DetProfile().E() / mcs.Start().E() < 0.8) {
	  return false;
	}
      }// if primary
    }

    return contained;
    
  }
  
} // namespace searchingfornues

#endif
