#ifndef GEOMETRYFUNCS_H
#define GEOMETRYFUNCS_H

#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

namespace searchingfornues
{
  float distance2d(const float& x1, const float& y1,
                  const float& x2, const float& y2)
  {
    return sqrt((x1-x2)*(x1-x2) +
                (y1-y2)*(y1-y2));
  }

  float distance3d(const float& x1, const float& y1, const float& z1,
                  const float& x2, const float& y2, const float& z2)
  {
    return sqrt((x1-x2)*(x1-x2) +
                (y1-y2)*(y1-y2) +
                (z1-z2)*(z1-z2));
  }

  float YZtoUcoordinate(const float& y, const float& z)
  {
    return -sqrt(3.)/2. * y + 1./2. * z;
  }

  float YZtoVcoordinate(const float& y, const float& z)
  {
    return +sqrt(3.)/2. * y + 1./2. * z;
  }

  float YZtoYcoordinate(const float& y, const float& z)
  {
    return z;
  }

  float YZtoPlanecoordinate(const float& y, const float& z, int plane)
  {
    if (plane == 0) return YZtoUcoordinate(y, z);
    if (plane == 1) return YZtoVcoordinate(y, z);
    if (plane == 2) return YZtoYcoordinate(y, z);
  }

} // namespace searchingfornues

#endif
