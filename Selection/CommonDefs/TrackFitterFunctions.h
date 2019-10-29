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
      if ((tkcalo->ResidualRange().back() - tkcalo->ResidualRange()[ic]) >= cmskip) {
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


  /**
   * @brief Return median angular deviation of a track in a given cm interval
   * @input trk  : track being provided as input
   * @input dmax : distance from track start point over which to calculate the angular deviations
   * @output median deflection angle [radins]
   */
  float GetTrackMedianDeflection(const searchingfornues::ProxyCaloElem_t& trk,
				 const float dmax) {
    
    float medangle = 0.;
    std::vector<float> dir_v; // vector of all directions
    //int npoints = 0; // how many direction vectors have we looked at?

    // store vertex of track
    // needed to understand if we are within distance dmax

    auto vtx = trk->Vertex();

    TVector3 trkdir(0,0,0);
    
    for(size_t i=0; i < trk->NumberTrajectoryPoints(); i++) {
      
      if (trk->HasValidPoint(i)) { // check this point is valid
	auto pt = trk->LocationAtPoint(i);
	
	auto dvtx = sqrt( ( (pt.X() - vtx.X()) * (pt.X() - vtx.X()) ) +
			  ( (pt.Y() - vtx.Y()) * (pt.Y() - vtx.Y()) ) +
			  ( (pt.Z() - vtx.Z()) * (pt.Z() - vtx.Z()) ) );
	
	// quit if we've passed maximum distance
	if (dvtx > dmax) break;
	
	auto dir = trk->DirectionAtPoint(i);

	TVector3 thisdir( dir.X(), dir.Y(), dir.Z() );

	// if we are at the first point, nothing to compare to
	if ( (trkdir.X() == 0) && (trkdir.X() == 0) && (trkdir.X() == 0) ) {
	  trkdir = thisdir;
	}
	// if we found a new point and can compare to the previous!
	else {
	  float dot = trkdir.Dot(thisdir) / ( trkdir.Mag() * thisdir.Mag() );
	  trkdir = thisdir;
	  float angle = acos(dot);
	  if (angle > 1e-5) 
	    dir_v.push_back(angle);
	}
      }// if point is valid
    }// for all track points
    
    std::sort(dir_v.begin(), dir_v.end());

    if (dir_v.size() == 0) return 0.;
      
    if (dir_v.size() % 2 == 1)
	medangle = dir_v[dir_v.size() / 2];
      else
	medangle = 0.5 * (dir_v[dir_v.size() / 2] + dir_v[dir_v.size() / 2 - 1]);

    //std::cout << "DAVIDC median angle : " << medangle << std::endl;
    
    return medangle;

  }
  
} // namespace searchingfornues

#endif
