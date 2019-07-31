#ifndef ANALYSIS_VERTEXANALYSIS_CXX
#define ANALYSIS_VERTEXANALYSIS_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "../CommonDefs/Typedefs.h"

#include "larreco/RecoAlg/Geometric3DVertexFitter.h"

namespace analysis
{
////////////////////////////////////////////////////////////////////////
//
// Class:       VertexAnalysis
// File:        VertexAnalysis.cc
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

class VertexAnalysis : public AnalysisToolBase
{

public:
  /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
  VertexAnalysis(const fhicl::ParameterSet &pset);

  /**
     *  @brief  Destructor
     */
  ~VertexAnalysis(){};

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
     * @brief set branches for TTree
     */
  void setBranches(TTree *_tree) override;

  /**
     * @brief reset ttree branches
     */
  void resetTTree(TTree *_tree) override;

private:
  art::InputTag fTRKFITproducer;
  trkf::Geometric3DVertexFitter *fitter;

  int _n_tracks_pandora;
  bool _vtx_fit_pandora_is_valid;
  float _vtx_fit_pandora_x;
  float _vtx_fit_pandora_y;
  float _vtx_fit_pandora_z;

  int _n_tracks_tkfit;
  bool _vtx_fit_tkfit_is_valid;
  float _vtx_fit_tkfit_x;
  float _vtx_fit_tkfit_y;
  float _vtx_fit_tkfit_z;
};

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///

VertexAnalysis::VertexAnalysis(const fhicl::ParameterSet &p)
{
  fTRKFITproducer = p.get<art::InputTag>("TRKFITproducer");
  fhicl::Table<trkf::Geometric3DVertexFitter::Config> options(p.get<fhicl::ParameterSet>("options"));
  fhicl::Table<trkf::TrackStatePropagator::Config> propagator(p.get<fhicl::ParameterSet>("propagator"));
  fitter = new trkf::Geometric3DVertexFitter(options, propagator);
}


void VertexAnalysis::configure(fhicl::ParameterSet const &p)
{
}

void VertexAnalysis::analyzeEvent(art::Event const &e, bool fData)
{
}

void VertexAnalysis::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
{
  auto trk_fit_h = e.getValidHandle<std::vector<recob::Track>>(fTRKFITproducer);
  std::vector<art::Ptr<recob::Track> > trk_fit_v;
  art::fill_ptr_vector(trk_fit_v, trk_fit_h);

  std::vector< art::Ptr<recob::Track> > pandora_tracks;
  std::vector< art::Ptr<recob::Track> > fitted_tracks;

  // look for the neutrino
  size_t primary_index = 100000;
  for (auto pfp : slice_pfp_v)
  {
    if (pfp->IsPrimary())
    {
      primary_index = pfp->Self();
    }
  }
  std::cout << "[VERTEX ANALYZER] primary index " << primary_index << std::endl;

  // loop on the daughters of the neutrino
  for (auto pfp : slice_pfp_v)
  {
    if (pfp->Parent() != primary_index)
      continue;

    // fill pandora
    auto trk_v = pfp.get<recob::Track>();
    if (trk_v.size() == 1)
    {
      auto trk = trk_v.at(0);
      pandora_tracks.push_back(trk);
    }

    // fill trk fit
    for (const auto tk: trk_fit_v)
    {
      if (tk->ID() == int(pfp.index()))
      {
        fitted_tracks.push_back(tk);
      }
    }
  }

  _n_tracks_pandora = pandora_tracks.size();
  _n_tracks_tkfit = fitted_tracks.size();

  if (_n_tracks_pandora > 1)
  {
    trkf::VertexWrapper track_vtx = fitter->fitTracks(pandora_tracks);
    _vtx_fit_pandora_is_valid = track_vtx.isValid();
    if (_vtx_fit_pandora_is_valid)
    {
      _vtx_fit_pandora_x = track_vtx.position().X();
      _vtx_fit_pandora_y = track_vtx.position().Y();
      _vtx_fit_pandora_z = track_vtx.position().Z();
    }
  }
  if (_n_tracks_tkfit > 1)
  {
    trkf::VertexWrapper trackfit_vtx = fitter->fitTracks(fitted_tracks);
    _vtx_fit_tkfit_is_valid = trackfit_vtx.isValid();
    if (_vtx_fit_tkfit_is_valid)
    {
      _vtx_fit_tkfit_x = trackfit_vtx.position().X();
      _vtx_fit_tkfit_y = trackfit_vtx.position().Y();
      _vtx_fit_tkfit_z = trackfit_vtx.position().Z();
    }
  }
}

void VertexAnalysis::setBranches(TTree *_tree)
{
  _tree->Branch("n_tracks_pandora", &_n_tracks_pandora, "n_tracks_pandora/i");
  _tree->Branch("vtx_fit_pandora_is_valid", &_vtx_fit_pandora_is_valid, "vtx_fit_pandora_is_valid/O");
  _tree->Branch("vtx_fit_pandora_x", &_vtx_fit_pandora_x, "vtx_fit_pandora_x/f");
  _tree->Branch("vtx_fit_pandora_y", &_vtx_fit_pandora_y, "vtx_fit_pandora_y/f");
  _tree->Branch("vtx_fit_pandora_z", &_vtx_fit_pandora_z, "vtx_fit_pandora_z/f");

  _tree->Branch("n_tracks_tkfit", &_n_tracks_tkfit, "n_tracks_tkfit/i");
  _tree->Branch("vtx_fit_tkfit_is_valid", &_vtx_fit_tkfit_is_valid, "vtx_fit_tkfit_is_valid/O");
  _tree->Branch("vtx_fit_tkfit_x", &_vtx_fit_tkfit_x, "vtx_fit_tkfit_x/f");
  _tree->Branch("vtx_fit_tkfit_y", &_vtx_fit_tkfit_y, "vtx_fit_tkfit_y/f");
  _tree->Branch("vtx_fit_tkfit_z", &_vtx_fit_tkfit_z, "vtx_fit_tkfit_z/f");
}

void VertexAnalysis::resetTTree(TTree *_tree)
{
  _n_tracks_pandora = -1;
  _vtx_fit_pandora_is_valid = false;
  _vtx_fit_pandora_x = -10000;
  _vtx_fit_pandora_y = -10000;
  _vtx_fit_pandora_z = -10000;

  _n_tracks_tkfit = -1;
  _vtx_fit_tkfit_is_valid = false;
  _vtx_fit_tkfit_x = -10000;
  _vtx_fit_tkfit_y = -10000;
  _vtx_fit_tkfit_z = -10000;
}

DEFINE_ART_CLASS_TOOL(VertexAnalysis)
} // namespace analysis

#endif
