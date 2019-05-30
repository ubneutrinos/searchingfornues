#ifndef ANALYSIS_SHOWERSTARTPOINT_CXX
#define ANALYSIS_SHOWERSTARTPOINT_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "ubobj/CRT/CRTHit.hh"

#include "nusimdata/SimulationBase/MCTruth.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// backtracking tools
#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"
#include "../CommonDefs/Geometry.h"

namespace analysis
{
////////////////////////////////////////////////////////////////////////
//
// Class:       ShowerStartPoint
// File:        ShowerStartPoint.cc
//
//              A basic analysis example
//
// Configuration parameters:
//
// TBD
//
// Created by Giuseppe Cerati (cerati@fnal.gov) on 03/15/2019
//
////////////////////////////////////////////////////////////////////////

class ShowerStartPoint : public AnalysisToolBase
{

public:
  /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
  ShowerStartPoint(const fhicl::ParameterSet &pset);

  /**
     *  @brief  Destructor
     */
  ~ShowerStartPoint(){};

  // provide for initialization
  void configure(fhicl::ParameterSet const &pset);

  /**
     * @brief Analysis function
     */
  void analyzeEvent(art::Event const &e, bool fData) override;

  /**
     * @brief Analyze slice
     */
  void analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;

  /**
     * @brief fill default for TTree
     */
  void fillDefault();

  /**
     * @brief set branches for TTree
     */
  void setBranches(TTree *_tree) override;

  /**
     * @brief reset ttree branches
     */
  void resetTTree(TTree *_tree) override;

private:
  std::vector<float> _shr_spacepoint_start_x_v;
  std::vector<float> _shr_spacepoint_start_y_v;
  std::vector<float> _shr_spacepoint_start_z_v;

  std::vector<float> _shr_hits_start_U_wire_v;
  std::vector<float> _shr_hits_start_U_x_v;
  std::vector<float> _shr_hits_start_V_wire_v;
  std::vector<float> _shr_hits_start_V_x_v;
  std::vector<float> _shr_hits_start_Y_wire_v;
  std::vector<float> _shr_hits_start_Y_x_v;

  detinfo::DetectorProperties const *detprop;

  float wireSpacing;
};

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
ShowerStartPoint::ShowerStartPoint(const fhicl::ParameterSet &p)
{
  detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  wireSpacing = 0.3;
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void ShowerStartPoint::configure(fhicl::ParameterSet const &p)
{
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void ShowerStartPoint::analyzeEvent(art::Event const &e, bool fData)
{
  std::cout << "analyze event" << std::endl;
}

void ShowerStartPoint::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
{
  Double_t reco_vtx[3] = {};

  for (size_t i_pfp = 0; i_pfp < slice_pfp_v.size(); i_pfp++)
  {
    auto PDG = fabs(slice_pfp_v[i_pfp]->PdgCode());

    if (PDG == 12 || PDG == 14)
    {
      // grab vertex
      auto vtx = slice_pfp_v[i_pfp].get<recob::Vertex>();
      if (vtx.size() != 1)
      {
        std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
      }
      else
      {
        // save vertex to array
        vtx.at(0)->XYZ(reco_vtx);
      }
      break;
    }
  }

  for (size_t i_pfp = 0; i_pfp < slice_pfp_v.size(); i_pfp++)
  {
    auto PDG = fabs(slice_pfp_v[i_pfp]->PdgCode());

    if (PDG == 12 || PDG == 14)
      continue;

    fillDefault();
    // spacepoints start point
    auto spacepoints = slice_pfp_v[i_pfp].get<recob::SpacePoint>();
    if (spacepoints.size() == 0)
    {
      continue;
    }
    else
    {
      //loop on spacepoints
      float smallest_sp_distance = std::numeric_limits<float>::max();
      size_t index_smallest_distance;
      for (size_t i = 0; i < spacepoints.size(); i++)
      {
        const auto &spacepoint = spacepoints[i];

        float sp_x = spacepoint->XYZ()[0];
        float sp_y = spacepoint->XYZ()[1];
        float sp_z = spacepoint->XYZ()[2];

        float distance_wrt_vertex = searchingfornues::distance3d(sp_x, sp_y, sp_z,
                                               reco_vtx[0], reco_vtx[1], reco_vtx[2]);

        if (distance_wrt_vertex < smallest_sp_distance)
        {
          smallest_sp_distance = distance_wrt_vertex;
          index_smallest_distance = i;
        }
      }

      _shr_spacepoint_start_x_v.back() = spacepoints[index_smallest_distance]->XYZ()[0];
      _shr_spacepoint_start_y_v.back() = spacepoints[index_smallest_distance]->XYZ()[1];
      _shr_spacepoint_start_z_v.back() = spacepoints[index_smallest_distance]->XYZ()[2];
    }

    // cluster per pfparticle
    auto clusters = slice_pfp_v[i_pfp].get<recob::Cluster>();
    for (const auto cluster : clusters)
    {
      int i_plane = cluster->Plane().Plane;
      float wire_reco_vtx = searchingfornues::YZtoPlanecoordinate(reco_vtx[1], reco_vtx[2], i_plane);

      //loop on the hits
      float smallest_hit_distance = std::numeric_limits<float>::max();
      float wire_coord_hit_min = std::numeric_limits<float>::max();
      float x_hit_min = std::numeric_limits<float>::max();
      auto hits = cluster.get<recob::Hit>();
      for (size_t i = 0; i < hits.size(); i++)
      {
        const auto &hit = hits[i];

        float wire_coord_hit = hit->WireID().Wire * wireSpacing;
        float x_hit = detprop->ConvertTicksToX(hit->PeakTime(), cluster->Plane());

        float distance_wrt_vertex = searchingfornues::distance2d(x_hit, wire_coord_hit,
                                               reco_vtx[0], wire_reco_vtx);

        if (distance_wrt_vertex < smallest_hit_distance)
        {
          smallest_hit_distance = distance_wrt_vertex;
          wire_coord_hit_min = wire_coord_hit;
          x_hit_min = x_hit;
        }
      }

      if (i_plane == 0)
      {
        _shr_hits_start_U_wire_v.back() = wire_coord_hit_min;
        _shr_hits_start_U_x_v.back() = x_hit_min;
      }
      if (i_plane == 1)
      {
        _shr_hits_start_V_wire_v.back() = wire_coord_hit_min;
        _shr_hits_start_V_x_v.back() = x_hit_min;
      }
      if (i_plane == 2)
      {
        _shr_hits_start_Y_wire_v.back() = wire_coord_hit_min;
        _shr_hits_start_Y_x_v.back() = x_hit_min;
      }
    }
  }
}

void ShowerStartPoint::fillDefault()
{
  _shr_spacepoint_start_x_v.push_back(std::numeric_limits<float>::lowest());
  _shr_spacepoint_start_y_v.push_back(std::numeric_limits<float>::lowest());
  _shr_spacepoint_start_z_v.push_back(std::numeric_limits<float>::lowest());

  _shr_hits_start_U_wire_v.push_back(std::numeric_limits<float>::lowest());
  _shr_hits_start_U_x_v.push_back(std::numeric_limits<float>::lowest());
  _shr_hits_start_V_wire_v.push_back(std::numeric_limits<float>::lowest());
  _shr_hits_start_V_x_v.push_back(std::numeric_limits<float>::lowest());
  _shr_hits_start_Y_wire_v.push_back(std::numeric_limits<float>::lowest());
  _shr_hits_start_Y_x_v.push_back(std::numeric_limits<float>::lowest());
}

void ShowerStartPoint::setBranches(TTree *_tree)
{
  _tree->Branch("shr_spacepoint_start_x_v", "std::vector<float>", &_shr_spacepoint_start_x_v);
  _tree->Branch("shr_spacepoint_start_y_v", "std::vector<float>", &_shr_spacepoint_start_y_v);
  _tree->Branch("shr_spacepoint_start_z_v", "std::vector<float>", &_shr_spacepoint_start_z_v);

  _tree->Branch("shr_hits_start_U_wire_v", "std::vector<float>", &_shr_hits_start_U_wire_v);
  _tree->Branch("shr_hits_start_U_x_v", "std::vector<float>", &_shr_hits_start_U_x_v);
  _tree->Branch("shr_hits_start_V_wire_v", "std::vector<float>", &_shr_hits_start_V_wire_v);
  _tree->Branch("shr_hits_start_V_x_v", "std::vector<float>", &_shr_hits_start_V_x_v);
  _tree->Branch("shr_hits_start_Y_wire_v", "std::vector<float>", &_shr_hits_start_Y_wire_v);
  _tree->Branch("shr_hits_start_Y_x_v", "std::vector<float>", &_shr_hits_start_Y_x_v);
}

void ShowerStartPoint::resetTTree(TTree *_tree)
{
  _shr_spacepoint_start_x_v.clear();
  _shr_spacepoint_start_y_v.clear();
  _shr_spacepoint_start_z_v.clear();

  _shr_hits_start_U_wire_v.clear();
  _shr_hits_start_U_x_v.clear();
  _shr_hits_start_V_wire_v.clear();
  _shr_hits_start_V_x_v.clear();
  _shr_hits_start_Y_wire_v.clear();
  _shr_hits_start_Y_x_v.clear();
}

DEFINE_ART_CLASS_TOOL(ShowerStartPoint)
} // namespace analysis

#endif
