#ifndef TRACKFITTERFUNCTIONS_H
#define TRACKFITTERFUNCTIONS_H

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

namespace searchingfornues
{

  /**
   * @brief Return dE/dx in MeV/cm for shower track-fit
   * @input tkcalo     : calorimetry object associated to fitted track
   * @input cmskip     : centimeters to skip from vertex for dE/dx calculation
   * @input cmlen      : centimeters over which dE/dx should be calculated (starting from cmskip)
   * @input localdEdx  : use local dE/dx from calorimetry? True -> yes. False: constant Q -> MeV conversion
   * @output dE/dx in MeV/cm
   * @output number of hits used for dE/dx calculation
   */
  void GetTrackFitdEdx(const art::Ptr<anab::Calorimetry>& tkcalo,
		       const float& cmskip, const float& cmlen, const bool& localdEdx,
		       float& dedx, int& nhits) {

    
    if (tkcalo->ResidualRange().size() == 0) {
      dedx  = std::numeric_limits<float>::lowest();
      nhits = std::numeric_limits<int>::lowest();
      return;
    }// if no points with which to calculate dE/dx

    // vector where to store dE/dx
    std::vector<float> dedxNcm;
    
    for (size_t ic = 0; ic < tkcalo->ResidualRange().size(); ++ic) {
      // check that at least cmskip cm away from vertex
      if ((tkcalo->ResidualRange().back() - tkcalo->ResidualRange()[ic]) > cmskip) {
	// check that no more then cmskip + cmlen away from vertex
	if ((tkcalo->ResidualRange().back() - tkcalo->ResidualRange()[ic]) < (cmskip + cmlen) ) {
	  if (localdEdx) 
	    dedxNcm.push_back(tkcalo->dEdx()[ic]); // if we want to use local dEdx from calo
	  else
	    dedxNcm.push_back(tkcalo->dQdx()[ic]); // if we should use Q and convert to MeV
	}
      }
    }

    float dedxNcm_med = -1.;
    
    if (dedxNcm.size() > 0) {
      
      std::sort(dedxNcm.begin(), dedxNcm.end());
      
      if (dedxNcm.size() % 2 == 1)
	dedxNcm_med = dedxNcm[dedxNcm.size() / 2];
      else
	dedxNcm_med = 0.5 * (dedxNcm[dedxNcm.size() / 2] + dedxNcm[dedxNcm.size() / 2 - 1]);
    }// if dedx vector has at least one netry

    dedx  = dedxNcm_med;
    nhits = dedxNcm.size();

    return;
  }
  
} // namespace searchingfornues

#endif
