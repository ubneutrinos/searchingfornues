#ifndef GEOMETRYFUNCS_H
#define GEOMETRYFUNCS_H

#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

namespace searchingfornues
{
  float distance2d(float& x1, float& y1,
                  float& x2, float& y2)
  {
    return sqrt((x1-x2)*(x1-x2) +
                (y1-y2)*(y1-y2));
  }

  float distance3d(float& x1, float& y1, float& z1,
                  float& x2, float& y2, float& z2)
  {
    return sqrt((x1-x2)*(x1-x2) +
                (y1-y2)*(y1-y2) +
                (z1-z2)*(z1-z2));
  }

  float YZtoUcoordinate(float& y, float& z)
  {
    return -sqrt(3.)/2. * y + 1./2. * z;
  }

  float YZtoVcoordinate(float& y, float& z)
  {
    return +sqrt(3.)/2. * y + 1./2. * z;
  }

  float YZtoYcoordinate(float& y, float& z)
  {
    return z;
  }

} // namespace searchingfornues

#endif
