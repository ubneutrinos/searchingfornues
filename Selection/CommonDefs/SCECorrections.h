#ifndef SCECORRECTIONSFUNCS_H
#define SCECORRECTIONSFUNCS_H

#include "larevt/SpaceChargeServices/SpaceChargeService.h"

namespace searchingfornues
{

  // apply the mapping of XYZ true -> XYZ position after SCE-induced shift.
  // to be applied to truth xyz in order to compare to reconstructed variables
  // e.g. used for resolution plots
  void ApplySCEMappingXYZ(float& x, float& y, float& z) {

    auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

    if (SCE->EnableSpatialSCE() == true) {
    
      auto offset = SCE->GetPosOffsets(geo::Point_t(vtx_x,vtx_y,vtx_z));
      x -= offset.X();
      y += offset.Y();
      z += offset.Z();

    }
    
  }

  // apply the SCE corrections to a reconstructed XYZ to see where the
  // XYZ position associated to the actual energy deposition should be
  // to be applied to reconstructed quantities to get a better XYZ coordinate.
  void ApplySCECorrectionXYZ(float& x, float& y, float& z) {

    auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
    
    if (SCE->EnableCalSpatialSCE() == true) {
      
      auto offset = SCE->GetCalPosOffsets(geo::Point_t(x,y,z));
      x -= offset.X();
      y += offset.Y();
      z += offset.Z();
      
    }// if spatial offset calibrations are enabled

  }
  
  return;
}

} // namespace searchingfornues

#endif
