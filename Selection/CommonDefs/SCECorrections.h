#ifndef SCECORRECTIONSFUNCS_H
#define SCECORRECTIONSFUNCS_H

#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

namespace searchingfornues
{

  // apply the mapping of XYZ true -> XYZ position after SCE-induced shift.
  // to be applied to truth xyz in order to compare to reconstructed variables
  // e.g. used for resolution plots
  void ApplySCEMappingXYZ(float& x, float& y, float& z)
  {
    auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

    if (SCE->EnableSimSpatialSCE() == true)
    {
      auto offset = SCE->GetPosOffsets(geo::Point_t(x, y, z));
      x -= offset.X();
      y += offset.Y();
      z += offset.Z();
    }
  }

  // apply the SCE corrections to a reconstructed XYZ to see where the
  // XYZ position associated to the actual energy deposition should be
  // to be applied to reconstructed quantities to get a better XYZ coordinate.
  void ApplySCECorrectionXYZ(float& x, float& y, float& z)
  {
    auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

    if (SCE->EnableCalSpatialSCE() == true)
    {

      auto offset = SCE->GetCalPosOffsets(geo::Point_t(x, y, z));
      x -= offset.X();
      y += offset.Y();
      z += offset.Z();
    }// if spatial offset calibrations are enabled
  }

  // apply the mapping of XYZ true -> XYZ position as it would be recosntructed.
  // takes into account SCE, trigger time offset, and wirecell-pandora offset.
  // to be applied to truth xyz in order to compare to reconstructed variables
  // e.g. used for resolution plots
  void True2RecoMappingXYZ(float& t, float& x, float& y, float& z)
  {
    ApplySCEMappingXYZ(x, y, z);

    auto const &detProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    auto const &detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
    double g4Ticks = detClocks->TPCG4Time2Tick(t) + detProperties->GetXTicksOffset(0, 0, 0) - detProperties->TriggerOffset();
    float _xtimeoffset = detProperties->ConvertTicksToX(g4Ticks, 0, 0, 0);

    x += _xtimeoffset;
    x += 0.6;
  }
} // namespace searchingfornues

#endif
