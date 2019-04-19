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
  art::InputTag fTRKproducer;
  art::InputTag fCALproducer;

  unsigned int _n_showers;

  std::vector<std::vector<double>> _shr_energy;
  std::vector<std::vector<double>> _shr_dedx;
  std::vector<size_t> _shr_pfp_id;

  std::vector<double> _shr_start_x;
  std::vector<double> _shr_start_y;
  std::vector<double> _shr_start_z;

  std::vector<double> _shr_px;
  std::vector<double> _shr_py;
  std::vector<double> _shr_pz;

  std::vector<double> _shr_theta;
  std::vector<double> _shr_phi;
  std::vector<double> _trkshr_score;

  std::vector<int> _shr_tkfit_nhits;
  std::vector<double> _shr_tkfit_start_x;
  std::vector<double> _shr_tkfit_start_y;
  std::vector<double> _shr_tkfit_start_z;
  std::vector<double> _shr_tkfit_theta;
  std::vector<double> _shr_tkfit_phi;
  std::vector<std::vector<double>> _shr_tkfit_dedx;
  std::vector<std::vector<int>> _shr_tkfit_dedx_nhits;
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

  searchingfornues::ProxyCaloColl_t const *tkcalo_proxy = NULL;
  if (fTRKproducer != "")
  {
    tkcalo_proxy = new searchingfornues::ProxyCaloColl_t(proxy::getCollection<std::vector<recob::Track>>(e, fTRKproducer, proxy::withAssociated<anab::Calorimetry>(fCALproducer)));
  }

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

      //fill dummy track fit values, overwrite them later
      _shr_tkfit_nhits.push_back(std::numeric_limits<double>::min());
      _shr_tkfit_start_x.push_back(std::numeric_limits<double>::min());
      _shr_tkfit_start_y.push_back(std::numeric_limits<double>::min());
      _shr_tkfit_start_z.push_back(std::numeric_limits<double>::min());
      _shr_tkfit_phi.push_back(std::numeric_limits<double>::min());
      _shr_tkfit_theta.push_back(std::numeric_limits<double>::min());
      _shr_tkfit_dedx.push_back(std::vector<double>(3, std::numeric_limits<double>::min()));
      _shr_tkfit_dedx_nhits.push_back(std::vector<int>(3, -1));

      if (tkcalo_proxy == NULL)
        continue;

      for (const searchingfornues::ProxyCaloElem_t &tk : *tkcalo_proxy)
      {

        // find track with ID matching the pfp index (this convention apparently works only for shower fits...)
        if (tk->ID() != int(slice_pfp_v[i_pfp].index()))
          continue;

        _shr_tkfit_nhits.back() = tk->CountValidPoints();
        _shr_tkfit_start_x.back() = tk->Start().X();
        _shr_tkfit_start_y.back() = tk->Start().Y();
        _shr_tkfit_start_z.back() = tk->Start().Z();
        _shr_tkfit_phi.back() = tk->StartDirection().Phi();
        _shr_tkfit_theta.back() = tk->StartDirection().Theta();

        auto const tkcalos = tk.get<anab::Calorimetry>();
        for (const auto &tkcalo : tkcalos)
        {
          if (tkcalo->ResidualRange().size() == 0)
            continue;
          std::vector<float> dedx4cm;
          for (size_t ic = 0; ic < tkcalo->ResidualRange().size(); ++ic)
          {
            if ((tkcalo->ResidualRange().back() - tkcalo->ResidualRange()[ic]) < 4.)
            {
              dedx4cm.push_back(tkcalo->dEdx()[ic]);
            }
          }
          float dedx4cm_med = -1.;
          if (dedx4cm.size() > 0)
          {
            std::sort(dedx4cm.begin(), dedx4cm.end());
            if (dedx4cm.size() % 2 == 1)
              dedx4cm_med = dedx4cm[dedx4cm.size() / 2];
            else
              dedx4cm_med = 0.5 * (dedx4cm[dedx4cm.size() / 2] + dedx4cm[dedx4cm.size() / 2 - 1]);
          }
          _shr_tkfit_dedx.back()[tkcalo->PlaneID().Plane] = dedx4cm_med;
          _shr_tkfit_dedx_nhits.back()[tkcalo->PlaneID().Plane] = dedx4cm.size();
        }
      }
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

  _tree->Branch("shr_px", &_shr_px, "shr_px/F");
  _tree->Branch("shr_py", &_shr_py, "shr_py/F");
  _tree->Branch("shr_pz", &_shr_pz, "shr_pz/F");

  _tree->Branch("shr_theta", "std::vector< double >", &_shr_theta);
  _tree->Branch("shr_phi", "std::vector< double >", &_shr_phi);

  _tree->Branch("trkshr_score", "std::vector< double >", &_trkshr_score);

  _tree->Branch("n_showers", &_n_showers, "n_showers/i");

  _tree->Branch("shr_tkfit_nhits", "std::vector< int >", &_shr_tkfit_nhits);
  _tree->Branch("shr_tkfit_start_x", "std::vector< double >", &_shr_tkfit_start_x);
  _tree->Branch("shr_tkfit_start_y", "std::vector< double >", &_shr_tkfit_start_y);
  _tree->Branch("shr_tkfit_start_z", "std::vector< double >", &_shr_tkfit_start_z);
  _tree->Branch("shr_tkfit_theta", "std::vector< double >", &_shr_tkfit_theta);
  _tree->Branch("shr_tkfit_phi", "std::vector< double >", &_shr_tkfit_phi);
  _tree->Branch("shr_tkfit_dedx", "std::vector< std::vector< double > >", &_shr_tkfit_dedx);
  _tree->Branch("shr_tkfit_dedx_nhits", "std::vector< std::vector< int > >", &_shr_tkfit_dedx_nhits);
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

  _shr_px.clear();
  _shr_py.clear();
  _shr_pz.clear();

  _n_showers = 0;
  _trkshr_score.clear();

  _shr_tkfit_nhits.clear();
  _shr_tkfit_start_x.clear();
  _shr_tkfit_start_y.clear();
  _shr_tkfit_start_z.clear();
  _shr_tkfit_theta.clear();
  _shr_tkfit_phi.clear();
  _shr_tkfit_dedx.clear();
  _shr_tkfit_dedx_nhits.clear();
}

DEFINE_ART_CLASS_TOOL(ShowerAnalysis)
} // namespace analysis

#endif
