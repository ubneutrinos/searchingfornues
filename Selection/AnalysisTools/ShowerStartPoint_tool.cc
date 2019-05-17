#ifndef ANALYSIS_SHOWERSTARTPOINT_CXX
#define ANALYSIS_SHOWERSTARTPOINT_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "ubobj/CRT/CRTHit.hh"

#include "nusimdata/SimulationBase/MCTruth.h"

// backtracking tools
#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"
#include "../CommonDefs/Geometry.h"
#include "../CommonDefs/SCECorrections.h"

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
  art::InputTag fTRKproducer;
  art::InputTag fCALproducer;

  std::vector<float> _shr_true_start_x_v;
  std::vector<float> _shr_true_start_y_v;
  std::vector<float> _shr_true_start_z_v;

  std::vector<float> _shr_true_start_U_v;
  std::vector<float> _shr_true_start_V_v;
  std::vector<float> _shr_true_start_Y_v;

  std::vector<float> _shr_true_sce_start_x_v;
  std::vector<float> _shr_true_sce_start_y_v;
  std::vector<float> _shr_true_sce_start_z_v;

  std::vector<float> _shr_true_sce_start_U_v;
  std::vector<float> _shr_true_sce_start_V_v;
  std::vector<float> _shr_true_sce_start_Y_v;

  std::vector<float> _shr_spacepoint_start_x_v;
  std::vector<float> _shr_spacepoint_start_y_v;
  std::vector<float> _shr_spacepoint_start_z_v;

  std::vector<float> _shr_hits_start_U_v;
  std::vector<float> _shr_hits_start_V_v;
  std::vector<float> _shr_hits_start_Y_v;
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
  fTRKproducer = p.get<art::InputTag>("TRKproducer", "");
  fCALproducer = p.get<art::InputTag>("CALproducer", "");
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

    //fill by default
    fillDefault();
    // auto nshr = slice_pfp_v[i_pfp].get<recob::Shower>().size();
    // if (nshr == 1)
    // {
    //   //need association with true shower
    //   float true_start_x, true_start_y, true_start_z;
    //   float true_start_U = searchingfornues::YZtoUcoordinate(true_start_y, true_start_z);
    //   float true_start_V = searchingfornues::YZtoVcoordinate(true_start_y, true_start_z);
    //   float true_start_Y = searchingfornues::YZtoYcoordinate(true_start_y, true_start_z);
    //
    //   _shr_true_start_x_v.push_back(true_start_x);
    //   _shr_true_start_y_v.push_back(true_start_y);
    //   _shr_true_start_z_v.push_back(true_start_z);
    //   _shr_true_start_U_v.push_back(true_start_U);
    //   _shr_true_start_V_v.push_back(true_start_V);
    //   _shr_true_start_Y_v.push_back(true_start_Y);
    //
    //   searchingfornues::True2RecoMappingXYZ(true_start_x, true_start_y, true_start_z);
    //   true_start_U = searchingfornues::YZtoUcoordinate(true_start_y, true_start_z);
    //   true_start_V = searchingfornues::YZtoVcoordinate(true_start_y, true_start_z);
    //   true_start_Y = searchingfornues::YZtoYcoordinate(true_start_y, true_start_z);
    //
    //   _shr_true_sce_start_x_v.push_back(true_start_x);
    //   _shr_true_sce_start_y_v.push_back(true_start_y);
    //   _shr_true_sce_start_z_v.push_back(true_start_z);
    //   _shr_true_sce_start_U_v.push_back(true_start_U);
    //   _shr_true_sce_start_V_v.push_back(true_start_V);
    //   _shr_true_sce_start_Y_v.push_back(true_start_Y);
    // }

    auto spacepoints = slice_pfp_v[i_pfp].get<recob::SpacePoint>();

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

    _shr_spacepoint_start_x_v.push_back(spacepoints[index_smallest_distance]->XYZ()[0]);
    _shr_spacepoint_start_y_v.push_back(spacepoints[index_smallest_distance]->XYZ()[1]);
    _shr_spacepoint_start_z_v.push_back(spacepoints[index_smallest_distance]->XYZ()[2]);
  }
}

void ShowerStartPoint::fillDefault()
{
  _shr_true_start_x_v.push_back(std::numeric_limits<float>::lowest());
  _shr_true_start_y_v.push_back(std::numeric_limits<float>::lowest());
  _shr_true_start_z_v.push_back(std::numeric_limits<float>::lowest());
  _shr_true_start_U_v.push_back(std::numeric_limits<float>::lowest());
  _shr_true_start_V_v.push_back(std::numeric_limits<float>::lowest());
  _shr_true_start_Y_v.push_back(std::numeric_limits<float>::lowest());
  _shr_true_sce_start_x_v.push_back(std::numeric_limits<float>::lowest());
  _shr_true_sce_start_y_v.push_back(std::numeric_limits<float>::lowest());
  _shr_true_sce_start_z_v.push_back(std::numeric_limits<float>::lowest());
  _shr_true_sce_start_U_v.push_back(std::numeric_limits<float>::lowest());
  _shr_true_sce_start_V_v.push_back(std::numeric_limits<float>::lowest());
  _shr_true_sce_start_Y_v.push_back(std::numeric_limits<float>::lowest());
  _shr_spacepoint_start_x_v.push_back(std::numeric_limits<float>::lowest());
  _shr_spacepoint_start_y_v.push_back(std::numeric_limits<float>::lowest());
  _shr_spacepoint_start_z_v.push_back(std::numeric_limits<float>::lowest());
  _shr_hits_start_U_v.push_back(std::numeric_limits<float>::lowest());
  _shr_hits_start_V_v.push_back(std::numeric_limits<float>::lowest());
  _shr_hits_start_Y_v.push_back(std::numeric_limits<float>::lowest());
}

void ShowerStartPoint::setBranches(TTree *_tree)
{
  _tree->Branch("_shr_true_start_x_v", "std::vector<float>", &_shr_true_start_x_v);
  _tree->Branch("_shr_true_start_y_v", "std::vector<float>", &_shr_true_start_y_v);
  _tree->Branch("_shr_true_start_z_v", "std::vector<float>", &_shr_true_start_z_v);

  _tree->Branch("_shr_true_start_U_v", "std::vector<float>", &_shr_true_start_U_v);
  _tree->Branch("_shr_true_start_V_v", "std::vector<float>", &_shr_true_start_V_v);
  _tree->Branch("_shr_true_start_Y_v", "std::vector<float>", &_shr_true_start_Y_v);

  _tree->Branch("_shr_true_sce_start_x_v", "std::vector<float>", &_shr_true_sce_start_x_v);
  _tree->Branch("_shr_true_sce_start_y_v", "std::vector<float>", &_shr_true_sce_start_y_v);
  _tree->Branch("_shr_true_sce_start_z_v", "std::vector<float>", &_shr_true_sce_start_z_v);

  _tree->Branch("_shr_true_sce_start_U_v", "std::vector<float>", &_shr_true_sce_start_U_v);
  _tree->Branch("_shr_true_sce_start_V_v", "std::vector<float>", &_shr_true_sce_start_V_v);
  _tree->Branch("_shr_true_sce_start_Y_v", "std::vector<float>", &_shr_true_sce_start_Y_v);

  _tree->Branch("_shr_spacepoint_start_x_v", "std::vector<float>", &_shr_spacepoint_start_x_v);
  _tree->Branch("_shr_spacepoint_start_y_v", "std::vector<float>", &_shr_spacepoint_start_y_v);
  _tree->Branch("_shr_spacepoint_start_z_v", "std::vector<float>", &_shr_spacepoint_start_z_v);

  _tree->Branch("_shr_hits_start_U_v", "std::vector<float>", &_shr_hits_start_U_v);
  _tree->Branch("_shr_hits_start_V_v", "std::vector<float>", &_shr_hits_start_V_v);
  _tree->Branch("_shr_hits_start_Y_v", "std::vector<float>", &_shr_hits_start_Y_v);
}

void ShowerStartPoint::resetTTree(TTree *_tree)
{
  _shr_true_start_x_v.clear();
  _shr_true_start_y_v.clear();
  _shr_true_start_z_v.clear();
  _shr_true_start_U_v.clear();
  _shr_true_start_V_v.clear();
  _shr_true_start_Y_v.clear();
  _shr_true_sce_start_x_v.clear();
  _shr_true_sce_start_y_v.clear();
  _shr_true_sce_start_z_v.clear();
  _shr_true_sce_start_U_v.clear();
  _shr_true_sce_start_V_v.clear();
  _shr_true_sce_start_Y_v.clear();
  _shr_spacepoint_start_x_v.clear();
  _shr_spacepoint_start_y_v.clear();
  _shr_spacepoint_start_z_v.clear();
  _shr_hits_start_U_v.clear();
  _shr_hits_start_V_v.clear();
  _shr_hits_start_Y_v.clear();
}

DEFINE_ART_CLASS_TOOL(ShowerStartPoint)
} // namespace analysis

#endif
