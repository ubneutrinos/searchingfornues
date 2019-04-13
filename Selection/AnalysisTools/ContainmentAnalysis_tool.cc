#ifndef SELECTION_CONTAINMENTSELECTION_CXX
#define SELECTION_CONTAINMENTSELECTION_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"

// backtracking tools
#include "ubana/ubana/searchingfornues/Selection/CommonDefs/BacktrackingFuncs.h"
#include "larcore/Geometry/Geometry.h"

namespace analysis
{
////////////////////////////////////////////////////////////////////////
//
// Class:       ContainmentAnalysis
// File:        ContainmentAnalysis.cc
//
//              A basic selection example
//
// Configuration parameters:
//
// TBD
//
// Created by David Caratelli (davidc@fnal.gov) on 01/30/2019
//
////////////////////////////////////////////////////////////////////////

class ContainmentAnalysis : public AnalysisToolBase
{

public:
  /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
  ContainmentAnalysis(const fhicl::ParameterSet &pset);

  /**
     *  @brief  Destructor
     */
  ~ContainmentAnalysis(){};

  // provide for initialization
  void configure(fhicl::ParameterSet const &pset);

  /**
     * @brief Analysis function
     */
  void analyzeEvent(art::Event const &e, bool fData) override;

  /**
     * @brief Selection function
     */
  void analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;

  /**
     * @brief set branches for TTree
     */
  void setBranches(TTree *_tree);

  /**
     * @brief reset ttree branches
     */
  void resetTTree(TTree *_tree);

private:
  float DistFiducial(float x, float y, float z);
  void DistFiducialBoundaries(float x, float y, float z, std::vector<std::vector<double>> &dboundaries);

  float _FV; // FV boundary to apply

  // TTree variables
  float _rc_vtx_x, _rc_vtx_y, _rc_vtx_z;
  float _dvtx; // smallest distance between vertex and any boundary
  float _dtrk; // smallest distance between any track start/end point and any boundary

  std::vector<std::vector<double>> _dtrk_boundary; // smallest distance between any track start/end point and each boundary
  std::vector<std::vector<double>> _dvtx_boundary; // smallest distance between the neutrino vertex and each boundary
  std::vector<std::vector<double>> _dshr_boundary; // smallest distance between any shower start point and each boundary
};

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
ContainmentAnalysis::ContainmentAnalysis(const fhicl::ParameterSet &pset)
{
  configure(pset);
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void ContainmentAnalysis::configure(fhicl::ParameterSet const &pset)
{

  _FV = pset.get<float>("FV");
}

void ContainmentAnalysis::analyzeEvent(art::Event const &e, bool fData)
{

  return;
}

//----------------------------------------------------------------------------
/// selectEvent
///
/// Arguments:
///
/// art::Event
/// slice track pointer vector
/// slice shower pointer vector
///
void ContainmentAnalysis::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
{

  TVector3 nuvtx;
  Double_t xyz[3] = {};

  _dvtx = 1e3;
  _dtrk = 1e3;
  _dtrk_boundary.clear();
  _dvtx_boundary.clear();
  _dshr_boundary.clear();

  _dtrk_boundary.resize(3, std::vector<double>(2, std::numeric_limits<double>::max()));
  _dvtx_boundary.resize(3, std::vector<double>(2, std::numeric_limits<double>::max()));
  _dshr_boundary.resize(3, std::vector<double>(2, std::numeric_limits<double>::max()));

  for (const auto &pfp_pxy : slice_pfp_v)
  {

    auto PDG = fabs(pfp_pxy->PdgCode());

    if ((PDG == 12) || (PDG == 14))
    {

      // grab vertex
      auto vtx = pfp_pxy.get<recob::Vertex>();
      if (vtx.size() != 1)
      {
        std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
        return;
      }

      // save vertex to array
      vtx.at(0)->XYZ(xyz);
      nuvtx = TVector3(xyz[0], xyz[1], xyz[2]);

      _dvtx = DistFiducial(nuvtx.X(), nuvtx.Y(), nuvtx.Z());
      DistFiducialBoundaries(xyz[0], xyz[1], xyz[2], _dvtx_boundary);

    } // if neutrino PFP

    else
    { // if not the neutrino PFP

      auto ntrk = pfp_pxy.get<recob::Track>().size();

      if (ntrk == 1)
      {

        auto trk = pfp_pxy.get<recob::Track>().at(0);

        auto trkstart = trk->Vertex();
        auto trkend = trk->End();

        float dstrt = DistFiducial(trkstart.X(), trkstart.Y(), trkstart.Z());
        if (dstrt < _dtrk)
          _dtrk = dstrt;

        float dend = DistFiducial(trkend.X(), trkend.Y(), trkend.Z());
        if (dend < _dtrk)
          _dtrk = dend;

        DistFiducialBoundaries(trkend.X(), trkend.Y(), trkend.Z(), _dtrk_boundary);
        DistFiducialBoundaries(trkstart.X(), trkstart.Y(), trkstart.Z(), _dtrk_boundary);

      }

      auto nshr = pfp_pxy.get<recob::Shower>().size();

      if (nshr == 1) {
        auto shr = pfp_pxy.get<recob::Shower>().at(0);
        auto shrstart = shr->ShowerStart();
        DistFiducialBoundaries(shrstart.X(), shrstart.Y(), shrstart.Z(), _dshr_boundary);
      }


    } // if not the neutrino PFP

  } // for all PFP

}

void ContainmentAnalysis::setBranches(TTree *_tree)
{

  _tree->Branch("rc_vtx_x", &_rc_vtx_x, "rc_vtx_x/F");
  _tree->Branch("rc_vtx_y", &_rc_vtx_y, "rc_vtx_y/F");
  _tree->Branch("rc_vtx_z", &_rc_vtx_z, "rc_vtx_z/F");
  _tree->Branch("dvtx", &_dvtx, "dvtx/F");
  _tree->Branch("dtrk", &_dtrk, "dtrk/F");
  _tree->Branch("dtrk_boundary", "std::vector < std::vector < double > >", &_dtrk_boundary);
  _tree->Branch("dvtx_boundary", "std::vector < std::vector < double > >", &_dvtx_boundary);
  _tree->Branch("dshr_boundary", "std::vector < std::vector < double > >", &_dshr_boundary);

  return;
}

void ContainmentAnalysis::resetTTree(TTree *_tree)
{

  _rc_vtx_x = std::numeric_limits<int>::min();
  _rc_vtx_y = std::numeric_limits<int>::min();
  _rc_vtx_z = std::numeric_limits<int>::min();

  _dvtx = std::numeric_limits<int>::max();
  _dtrk = std::numeric_limits<int>::max();

  return;
}

void ContainmentAnalysis::DistFiducialBoundaries(float x, float y, float z, std::vector< std::vector < double > > &dist_boundaries)
{

  art::ServiceHandle<geo::Geometry> geo;
  std::vector<double> bnd = {
      0., 2. * geo->DetHalfWidth(),
      -geo->DetHalfHeight(), geo->DetHalfHeight(),
      0., geo->DetLength()
    };

  auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

  if (SCE->EnableCalSpatialSCE() == true)
  {

    auto offset = SCE->GetPosOffsets(geo::Point_t(x, y, z));
    x += offset.X();
    y -= offset.Y();
    z -= offset.Z();

  } // if spatial offset calibrations are enabled

  std::cout << x << " " << bnd[0] << " " << x - bnd[0] <<std::endl;
  if (x - bnd[0] < dist_boundaries[0][0])
  {
    dist_boundaries[0][0] = x - bnd[0];
  }

  if (bnd[1] - x < dist_boundaries[0][1])
  {
    dist_boundaries[0][1] = bnd[1] - x;
  }

  if (y - bnd[2] < dist_boundaries[1][0])
  {
     dist_boundaries[1][0] = y - bnd[2];
  }

  if (bnd[3] - y < dist_boundaries[1][1])
  {

    dist_boundaries[1][1] = bnd[3] - y;
  }

  if (z - bnd[4] < dist_boundaries[2][0])
  {
    dist_boundaries[2][0] = z - bnd[4];
  }

  if (bnd[5] - z < dist_boundaries[2][1])
  {
    dist_boundaries[2][1] = bnd[5] - z;
  }
}

float ContainmentAnalysis::DistFiducial(float x, float y, float z)
{

  float dmin = 1e3;

  auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

  if (SCE->EnableCalSpatialSCE() == true)
  {

    auto offset = SCE->GetPosOffsets(geo::Point_t(x, y, z));
    x += offset.X();
    y -= offset.Y();
    z -= offset.Z();

  } // if spatial offset calibrations are enabled

  double dxl = x - 0.;
  if (dxl < dmin)
    dmin = dxl;

  double dxh = 256. - x;
  if (dxh < dmin)
    dmin = dxh;

  double dyl = y - -116.;
  if (dyl < dmin)
    dmin = dyl;

  double dyh = 116. - y;
  if (dyh < dmin)
    dmin = dyh;

  double dzl = z - 0.;
  if (dzl < dmin)
    dmin = dzl;

  double dzh = 1036. - z;
  if (dzh < dmin)
    dmin = dzh;

  return dmin;
}

DEFINE_ART_CLASS_TOOL(ContainmentAnalysis)
} // namespace analysis

#endif
