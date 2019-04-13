#ifndef ANALYSIS_SHOWERANALYSIS_CXX
#define ANALYSIS_SHOWERANALYSIS_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "ubobj/CRT/CRTHit.hh"
#include "lardataobj/RecoBase/OpFlash.h"

#include "nusimdata/SimulationBase/MCTruth.h"

// backtracking tools
#include "ubana/ubana/searchingfornues/Selection/CommonDefs/BacktrackingFuncs.h"
#include "ubana/ubana/searchingfornues/Selection/CommonDefs/TrackShowerScoreFuncs.h"

namespace analysis
{
////////////////////////////////////////////////////////////////////////
//
// Class:       ShowerAnalysis
// File:        ShowerAnalysis.cc
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

class ShowerAnalysis : public AnalysisToolBase
{

public:
  /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
  ShowerAnalysis(const fhicl::ParameterSet &pset);

  /**
     *  @brief  Destructor
     */
  ~ShowerAnalysis(){};

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
     * @brief Save truth info for event associated to neutrino
     */
  void SaveTruth(art::Event const &e);

  /**
     * @brief set branches for TTree
     */
  void setBranches(TTree *_tree) override;

  /**
     * @brief reset ttree branches
     */
  void resetTTree(TTree *_tree) override;

private:
  unsigned int _n_showers;

  std::vector<std::vector<double>> _shr_energy;
  std::vector<std::vector<double>> _shr_dedx;
  std::vector<size_t> _shr_pfp_id;

  std::vector<double> _shr_start_x;
  std::vector<double> _shr_start_y;
  std::vector<double> _shr_start_z;

  std::vector<double> _shr_theta;
  std::vector<double> _shr_phi;
  std::vector<double> _trkshr_score;
};

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
ShowerAnalysis::ShowerAnalysis(const fhicl::ParameterSet &p)
{

}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void ShowerAnalysis::configure(fhicl::ParameterSet const &p)
{
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void ShowerAnalysis::analyzeEvent(art::Event const &e, bool fData)
{
    std::cout << "analyze event" << std::endl;
}

void ShowerAnalysis::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
{

  for (size_t i_pfp = 0; i_pfp < slice_pfp_v.size(); i_pfp++)
  {
    auto PDG = fabs(slice_pfp_v[i_pfp]->PdgCode());

    if (PDG == 12 || PDG == 14)
      continue;

    for (const auto &shr : slice_pfp_v[i_pfp].get<recob::Shower>())
    {
      _n_showers++;
      _shr_dedx.push_back(shr->dEdx());
      _shr_energy.push_back(shr->Energy());
      _shr_pfp_id.push_back(i_pfp);
      _shr_phi.push_back(shr->Direction().Phi());
      _shr_theta.push_back(shr->Direction().Theta());
      _shr_start_x.push_back(shr->ShowerStart().X());
      _shr_start_y.push_back(shr->ShowerStart().Y());
      _shr_start_z.push_back(shr->ShowerStart().Z());
      _trkshr_score.push_back(searchingfornues::GetTrackShowerScore(slice_pfp_v[i_pfp]));
    }

  }

}

void ShowerAnalysis::setBranches(TTree *_tree)
{
  _tree->Branch("shr_dedx", "std::vector< std::vector< double > >", &_shr_dedx);
  _tree->Branch("shr_energy", "std::vector< std::vector< double > >", &_shr_energy);
  _tree->Branch("shr_pfp_id", "std::vector< size_t >", &_shr_pfp_id);

  _tree->Branch("shr_start_x", "std::vector< double >", &_shr_start_x);
  _tree->Branch("shr_start_y", "std::vector< double >", &_shr_start_y);
  _tree->Branch("shr_start_z", "std::vector< double >", &_shr_start_z);

  _tree->Branch("shr_theta", "std::vector< double >", &_shr_theta);
  _tree->Branch("shr_phi", "std::vector< double >", &_shr_phi);

  _tree->Branch("trkshr_score", "std::vector< double >", &_trkshr_score);

  _tree->Branch("n_showers", &_n_showers, "n_showers/i");
}

void ShowerAnalysis::resetTTree(TTree *_tree)
{
  _shr_energy.clear();
  _shr_dedx.clear();
  _shr_pfp_id.clear();

  _shr_start_x.clear();
  _shr_start_y.clear();
  _shr_start_z.clear();

  _shr_theta.clear();
  _shr_phi.clear();

  _n_showers = 0;
  _trkshr_score.clear();
}

DEFINE_ART_CLASS_TOOL(ShowerAnalysis)
} // namespace analysis

#endif