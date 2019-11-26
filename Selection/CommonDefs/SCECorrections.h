#ifndef SCECORRECTIONSFUNCS_H
#define SCECORRECTIONSFUNCS_H

#include "larevt/SpaceCharge/SpaceCharge.h"
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


  // apply the mapping of XYZ true -> XYZ position after SCE-induced shift.
  // to be applied to truth xyz in order to compare to reconstructed variables
  // e.g. used for resolution plots
  void ApplySCEMappingXYZ(float x, float y, float z, float out[3])
  {
    auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

    if (SCE->EnableSimSpatialSCE() == true)
    {
      auto offset = SCE->GetPosOffsets(geo::Point_t(x, y, z));
      out[0] = (x - offset.X());
      out[1] = (y + offset.Y());
      out[2] = (z + offset.Z());
    }
  }

  // apply the SCE corrections to a reconstructed XYZ to see where the
  // XYZ position associated to the actual energy deposition should be
  // to be applied to reconstructed quantities to get a better XYZ coordinate.
  void ApplySCECorrectionXYZ(float x, float y, float z, float out[3])
  {
    auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

    if (SCE->EnableCalSpatialSCE() == true)
    {
      auto offset = SCE->GetCalPosOffsets(geo::Point_t(x, y, z));
      out[0] = (x - offset.X());
      out[1] = (y + offset.Y());
      out[2] = (z + offset.Z());
    }// if spatial offset calibrations are enabled
  }

  float x_offset(float t)
  {
    auto const &detProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    auto const &detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
    double g4Ticks = detClocks->TPCG4Time2Tick(t) + detProperties->GetXTicksOffset(0, 0, 0) - detProperties->TriggerOffset();
    float xoffset = detProperties->ConvertTicksToX(g4Ticks, 0, 0, 0);
    xoffset += 0.6;
    return xoffset;
  }

  // apply the mapping of XYZ true -> XYZ position as it would be recosntructed.
  // takes into account SCE, trigger time offset, and wirecell-pandora offset.
  // to be applied to truth xyz in order to compare to reconstructed variables
  // e.g. used for resolution plots
  void True2RecoMappingXYZ(float t, float x, float y, float z, float out[3])
  {
    ApplySCEMappingXYZ(x, y, z, out);
    float _xoffset = x_offset(t);
    out[0] += _xoffset;
  }

  float GetLocalEFieldMag(const float x, const float y, const float z) {
    
    const detinfo::DetectorProperties* detprop = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
    auto const *sce = lar::providerFrom<spacecharge::SpaceChargeService>();
    
    double E_field_nominal = detprop->Efield();        // Electric Field in the drift region in KV/cm

    //correct Efield for SCE
    geo::Vector_t E_field_offsets = {0.,0.,0.};
    E_field_offsets = sce->GetCalEfieldOffsets(geo::Point_t{x,y, z});
    TVector3 E_field_vector = {E_field_nominal*(1 + E_field_offsets.X()), E_field_nominal*E_field_offsets.Y(), E_field_nominal*E_field_offsets.Z()};
    float E_field = E_field_vector.Mag();

    return E_field;
  }

} // namespace searchingfornues

#endif
