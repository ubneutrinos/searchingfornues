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


    if (tkcalo->ResidualRange().size() == 0)
    {
      dedx  = std::numeric_limits<float>::lowest();
      nhits = std::numeric_limits<int>::lowest();
      return;
    }// if no points with which to calculate dE/dx

    // vector where to store dE/dx
    std::vector<float> dedxNcm;

    for (size_t ic = 0; ic < tkcalo->ResidualRange().size(); ++ic)
    {
      // check that at least cmskip cm away from vertex
      if ((tkcalo->ResidualRange().back() - tkcalo->ResidualRange()[ic]) >= cmskip)
      {
        // check that no more then cmskip + cmlen away from vertex
        if ((tkcalo->ResidualRange().back() - tkcalo->ResidualRange()[ic]) < (cmskip + cmlen) )
        {
          if (localdEdx) dedxNcm.push_back(tkcalo->dEdx()[ic]); // if we want to use local dEdx from calo
          else dedxNcm.push_back(tkcalo->dQdx()[ic]); // if we should use Q and convert to MeV
        }
      }
    }

    float dedxNcm_med = -1.;

    if (dedxNcm.size() > 0)
    {
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
   * @brief Return dE/dx in MeV/cm for shower track-fit
   * @input dedx_v     : dedx values associated to fitted track
   * @input rr_v     : rr values associated to fitted track
   * @input cmskip     : centimeters to skip from vertex for dE/dx calculation
   * @input cmlen      : centimeters over which dE/dx should be calculated (starting from cmskip)
   * @input localdEdx  : use local dE/dx from calorimetry? True -> yes. False: constant Q -> MeV conversion
   * @output dE/dx in MeV/cm
   * @output number of hits used for dE/dx calculation
   */
  void GetTrackFitdEdx(const std::vector<float> dedx_v,
           const std::vector<float> rr_v,
           const float& cmskip, const float& cmlen,
           float& dedx, int& nhits) {


    if (rr_v.size() == 0)
    {
      dedx  = std::numeric_limits<float>::lowest();
      nhits = std::numeric_limits<int>::lowest();
      return;
    }// if no points with which to calculate dE/dx

    // vector where to store dE/dx
    std::vector<float> dedxNcm;

    for (size_t ic = 0; ic < rr_v.size(); ++ic)
    {
      // check that at least cmskip cm away from vertex
      if ((rr_v.back() - rr_v[ic]) >= cmskip)
      {
        // check that no more then cmskip + cmlen away from vertex
        if ((rr_v.back() - rr_v[ic]) < (cmskip + cmlen) )
        {
          dedxNcm.push_back(dedx_v[ic]);
        }
      }
    }

    float dedxNcm_med = -1.;

    if (dedxNcm.size() > 0)
    {

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
   * @brief Return dE/dx in MeV/cm for shower track-fit
   * @input tkcalo     : calorimetry object associated to fitted track
   * @input cmskip     : centimeters to skip from vertex for dE/dx calculation
   * @input cmlen      : centimeters over which dE/dx should be calculated (starting from cmskip)
   * @input localdEdx  : use local dE/dx from calorimetry? True -> yes. False: constant Q -> MeV conversion
   * @input shrstartx  : start x coordinate [with SCE corrections!]
   * @input shrstarty  : start y coordinate [with SCE corrections!]
   * @input shrstartz  : start z coordinate [with SCE corrections!]
   * @output dE/dx in MeV/cm
   * @output number of hits used for dE/dx calculation
   */
  void GetTrackFitdEdx(const art::Ptr<anab::Calorimetry>& tkcalo,
           const float& cmskip, const float& cmlen, const bool& localdEdx,
           const float& shrstartx, const float& shrstarty, const float& shrstartz,
           float& dedx, int& nhits) {


    if (tkcalo->ResidualRange().size() == 0)
    {
      dedx  = std::numeric_limits<float>::lowest();
      nhits = std::numeric_limits<int>::lowest();
      return;
    }// if no points with which to calculate dE/dx

    // vector where to store dE/dx
    std::vector<float> dedxNcm;

    for (size_t ic = 0; ic < tkcalo->XYZ().size(); ++ic)
    {

      // calculate 3D distance to start point
      float d3d = sqrt( ( pow(tkcalo->XYZ()[ic].X() - shrstartx, 2) ) +
      ( pow(tkcalo->XYZ()[ic].Y() - shrstarty, 2) ) +
      ( pow(tkcalo->XYZ()[ic].Z() - shrstartz, 2) ) );

      // check that at least cmskip cm away from vertex
      if ( d3d >= cmskip)
      {
        // check that no more then cmskip + cmlen away from vertex
        if ( d3d < (cmskip + cmlen) )
        {
          if (localdEdx)
            dedxNcm.push_back(tkcalo->dEdx()[ic]); // if we want to use local dEdx from calo
          else
            dedxNcm.push_back(tkcalo->dQdx()[ic]); // if we should use Q and convert to MeV
        }
      }
    }

    float dedxNcm_med = -1.;

    if (dedxNcm.size() > 0)
    {
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
   * @brief Return rms angular deviation of a track in a given cm interval
   * @input trk  : track being provided as input
   * @input dmax : distance from track start point over which to calculate the angular deviations
   * @output rms deflection angle [radins]
   */
  float GetTrackRMSDeflection(const searchingfornues::ProxyCaloElem_t& trk,
            const float dmax) {

    float medangle = 0.;
    float rmsangle = 0.;
    std::vector<float> dir_v; // vector of all directions
    //int npoints = 0; // how many direction vectors have we looked at?

    // store vertex of track
    // needed to understand if we are within distance dmax

    auto vtx = trk->Vertex();

    TVector3 trkdir(0,0,0);

    for(size_t i=0; i < trk->NumberTrajectoryPoints(); i++)
    {

      if (trk->HasValidPoint(i))
      { // check this point is valid
        auto pt = trk->LocationAtPoint(i);

        auto dvtx = sqrt( ( (pt.X() - vtx.X()) * (pt.X() - vtx.X()) ) +
        ( (pt.Y() - vtx.Y()) * (pt.Y() - vtx.Y()) ) +
        ( (pt.Z() - vtx.Z()) * (pt.Z() - vtx.Z()) ) );

        // quit if we've passed maximum distance
        if (dvtx > dmax) break;

        auto dir = trk->DirectionAtPoint(i);

        TVector3 thisdir( dir.X(), dir.Y(), dir.Z() );

        // if we are at the first point, nothing to compare to
        if ( (trkdir.X() == 0) && (trkdir.X() == 0) && (trkdir.X() == 0) )
        {
          trkdir = thisdir;
        }
        // if we found a new point and can compare to the previous!
        else
        {
          float dot = trkdir.Dot(thisdir) / ( trkdir.Mag() * thisdir.Mag() );
          trkdir = thisdir;
          float angle = acos(dot);
          if (angle > 1e-5)
            dir_v.push_back(angle);
        }
      }// if point is valid
    }// for all track points

    if (dir_v.size() == 0) return 0.;

    // calculate average...
    for (size_t d=0; d < dir_v.size(); d++)
      medangle += dir_v[d];
    medangle /= dir_v.size();
    // ... and RMS
    for (size_t d=0; d < dir_v.size(); d++)
      rmsangle += (dir_v[d] - medangle) * (dir_v[d] - medangle);
    rmsangle /= sqrt( dir_v.size() );

    return rmsangle;
  }

  /**
   * @brief Return rms angular deviation of a shower's spacepoints w.r.t. direction
   * @input shower proxy  : track being provided as input
   */
  void GetMoliereRadius(const searchingfornues::ProxyPfpElem_t pfp_pxy,
      float& medangle, float& rmsangle) {

    auto shower_v  = pfp_pxy.get<recob::Shower>();
    auto spcpnts_v = pfp_pxy.get<recob::SpacePoint>();

    if (shower_v.size() != 1) return;

    auto shower = shower_v[0];

    auto dir3D = shower->Direction();
    TVector3 shrdir(dir3D.X(), dir3D.Y(), dir3D.Z() );

    auto vtx3D = shower->ShowerStart();
    TVector3 shrvtx(vtx3D.X(), vtx3D.Y(), vtx3D.Z() );

    std::vector<float> angle_v;

    medangle = 0;
    rmsangle = 0;

    for (auto &sp : spcpnts_v)
    {
      auto spxyz = sp->XYZ();

      TVector3 sppos(spxyz[0],spxyz[1],spxyz[2]);

      TVector3 sptovtx = sppos-shrvtx;

      if (sptovtx.Mag() == 0) continue;

      float angle = acos( sptovtx.Dot( shrdir ) / ( sptovtx.Mag() * shrdir.Mag() ) );
      angle *= (180./3.14);

      angle_v.push_back( angle );

    }// for all spacepoints

    if (angle_v.size()==0) {
      medangle = std::numeric_limits<float>::max();
      rmsangle = std::numeric_limits<float>::max();
      return;
    }

    // calculate average...
    for (size_t d=0; d < angle_v.size(); d++)
    {
      medangle += angle_v[d];
    }
    medangle /= angle_v.size();
    // ... and RMS
    for (size_t d=0; d < angle_v.size(); d++)
    {
      rmsangle += (angle_v[d] - medangle) * (angle_v[d] - medangle);
    }
    rmsangle /= sqrt( angle_v.size() );

    return;

  }// end of function

  /**
   * @brief Return rms angular deviation of a shower's spacepoints w.r.t. direction; this version is for a shower to be merged with a second pfp.
   * @input shower proxy  : track being provided as input
   */
  void GetMoliereRadiusMergedShowers(const searchingfornues::ProxyPfpElem_t pfp_pxy, const searchingfornues::ProxyPfpElem_t pfp2_pxy,
      float& medangle, float& rmsangle) {

    auto shower_v  = pfp_pxy.get<recob::Shower>();
    auto spcpnts_v = pfp_pxy.get<recob::SpacePoint>();

    auto shower2_v  = pfp2_pxy.get<recob::Shower>();
    auto spcpnts2_v = pfp2_pxy.get<recob::SpacePoint>();

    if (shower_v.size() != 1) return;
    if (shower2_v.size() != 1) return;

    auto shower = shower_v[0];
    auto shower2 = shower2_v[0];

    auto dir3D = shower->Direction();
    TVector3 shrdir(dir3D.X(), dir3D.Y(), dir3D.Z() );

    auto vtx3D = shower->ShowerStart();
    TVector3 shrvtx(vtx3D.X(), vtx3D.Y(), vtx3D.Z() );

    //if the 2nd pfp is upstream of the shower, use it for vtx and dir
    auto vtx3D_shr2 = shower2->ShowerStart();
    if ( (vtx3D_shr2-vtx3D).Dot(shower->Direction()) < 0. ) {
      shrvtx = TVector3(vtx3D_shr2.X(), vtx3D_shr2.Y(), vtx3D_shr2.Z() );
      shrdir = TVector3(shower2->Direction().X(),shower2->Direction().Y(),shower2->Direction().Z());
    }

    std::vector<float> angle_v;

    medangle = 0;
    rmsangle = 0;

    for (auto &sp : spcpnts_v)
    {
      auto spxyz = sp->XYZ();

      TVector3 sppos(spxyz[0],spxyz[1],spxyz[2]);

      TVector3 sptovtx = sppos-shrvtx;

      if (sptovtx.Mag() == 0) continue;

      float angle = acos( sptovtx.Dot( shrdir ) / ( sptovtx.Mag() * shrdir.Mag() ) );
      angle *= (180./3.14);

      angle_v.push_back( angle );

    }// for all spacepoints

    //add also the spacepoints from the 2nd shower
    for (auto &sp : spcpnts2_v)
    {
      auto spxyz = sp->XYZ();

      TVector3 sppos(spxyz[0],spxyz[1],spxyz[2]);

      TVector3 sptovtx = sppos-shrvtx;

      if (sptovtx.Mag() == 0) continue;

      float angle = acos( sptovtx.Dot( shrdir ) / ( sptovtx.Mag() * shrdir.Mag() ) );
      angle *= (180./3.14);

      angle_v.push_back( angle );

    }// for all spacepoints

    if (angle_v.size()==0) {
      medangle = std::numeric_limits<float>::max();
      rmsangle = std::numeric_limits<float>::max();
      return;
    }

    // calculate average...
    for (size_t d=0; d < angle_v.size(); d++)
    {
      medangle += angle_v[d];
    }
    medangle /= angle_v.size();
    // ... and RMS
    for (size_t d=0; d < angle_v.size(); d++)
    {
      rmsangle += (angle_v[d] - medangle) * (angle_v[d] - medangle);
    }
    rmsangle /= sqrt( angle_v.size() );

    return;

  }// end of function


} // namespace searchingfornues

#endif
