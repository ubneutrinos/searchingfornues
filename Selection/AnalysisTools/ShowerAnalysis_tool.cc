#ifndef ANALYSIS_SHOWERANALYSIS_CXX
#define ANALYSIS_SHOWERANALYSIS_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "ubobj/CRT/CRTHit.hh"

#include "nusimdata/SimulationBase/MCTruth.h"

// backtracking tools
#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"

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

  std::vector<double> _shr_energy_u_v;
  std::vector<double> _shr_energy_v_v;
  std::vector<double> _shr_energy_y_v;

  std::vector<double> _shr_dedx_u_v;
  std::vector<double> _shr_dedx_v_v;
  std::vector<double> _shr_dedx_y_v;

  std::vector<size_t> _shr_pfp_id_v;

  std::vector<double> _shr_start_x_v;
  std::vector<double> _shr_start_y_v;
  std::vector<double> _shr_start_z_v;

  std::vector<double> _shr_dist_v;

  std::vector<double> _shr_px_v;
  std::vector<double> _shr_py_v;
  std::vector<double> _shr_pz_v;

  std::vector<double> _shr_theta_v;
  std::vector<double> _shr_phi_v;
  std::vector<double> _trkshr_score_v;

  std::vector<int> _shr_tkfit_nhits_v;
  std::vector<double> _shr_tkfit_start_x_v;
  std::vector<double> _shr_tkfit_start_y_v;
  std::vector<double> _shr_tkfit_start_z_v;
  std::vector<double> _shr_tkfit_theta_v;
  std::vector<double> _shr_tkfit_phi_v;
  std::vector<std::vector<double>> _shr_tkfit_dedx_v;
  std::vector<std::vector<int>> _shr_tkfit_dedx_nhits_v;
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

  TVector3 nuvtx;

  for (size_t i_pfp = 0; i_pfp < slice_pfp_v.size(); i_pfp++)
  {
    auto PDG = fabs(slice_pfp_v[i_pfp]->PdgCode());

    if (PDG == 12 || PDG == 14)
    {
      // grab vertex
      Double_t xyz[3] = {};

      auto vtx = slice_pfp_v[i_pfp].get<recob::Vertex>();
      if (vtx.size() != 1)
      {
        std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
      }
      else
      {
        // save vertex to array
        vtx.at(0)->XYZ(xyz);
        nuvtx.SetXYZ(xyz[0], xyz[1], xyz[2]);
      }

      break;
    }
  }

  for (size_t i_pfp = 0; i_pfp < slice_pfp_v.size(); i_pfp++)
  {
    auto PDG = fabs(slice_pfp_v[i_pfp]->PdgCode());

    if (PDG == 12 || PDG == 14)
      continue;

    for (const auto &shr : slice_pfp_v[i_pfp].get<recob::Shower>())
    {
      _n_showers++;
      _shr_dedx_u_v.push_back(shr->dEdx()[0]);
      _shr_dedx_v_v.push_back(shr->dEdx()[1]);
      _shr_dedx_y_v.push_back(shr->dEdx()[2]);

      _shr_energy_u_v.push_back(shr->Energy()[0]);
      _shr_energy_v_v.push_back(shr->Energy()[1]);
      _shr_energy_y_v.push_back(shr->Energy()[2]);

      _shr_pfp_id_v.push_back(i_pfp);
      _shr_phi_v.push_back(shr->Direction().Phi());
      _shr_theta_v.push_back(shr->Direction().Theta());
      _shr_start_x_v.push_back(shr->ShowerStart().X());
      _shr_start_y_v.push_back(shr->ShowerStart().Y());
      _shr_start_z_v.push_back(shr->ShowerStart().Z());
      _trkshr_score_v.push_back(searchingfornues::GetTrackShowerScore(slice_pfp_v[i_pfp]));

      //fill dummy track fit values, overwrite them later
      _shr_tkfit_nhits_v.push_back(std::numeric_limits<int>::lowest());
      _shr_tkfit_start_x_v.push_back(std::numeric_limits<double>::lowest());
      _shr_tkfit_start_y_v.push_back(std::numeric_limits<double>::lowest());
      _shr_tkfit_start_z_v.push_back(std::numeric_limits<double>::lowest());
      _shr_tkfit_phi_v.push_back(std::numeric_limits<double>::lowest());
      _shr_tkfit_theta_v.push_back(std::numeric_limits<double>::lowest());
      _shr_tkfit_dedx_v.push_back(std::vector<double>(3, std::numeric_limits<double>::lowest()));
      _shr_tkfit_dedx_nhits_v.push_back(std::vector<int>(3, -1));
      _shr_px_v.push_back(shr->Direction().X());
      _shr_py_v.push_back(shr->Direction().Y());
      _shr_pz_v.push_back(shr->Direction().Z());

      _shr_dist_v.push_back((shr->ShowerStart() - nuvtx).Mag());

      if (tkcalo_proxy == NULL)
        continue;

      for (const searchingfornues::ProxyCaloElem_t &tk : *tkcalo_proxy)
      {

        // find track with ID matching the pfp index (this convention apparently works only for shower fits...)
        if (tk->ID() != int(slice_pfp_v[i_pfp].index()))
          continue;

        _shr_tkfit_nhits_v.back() = tk->CountValidPoints();
        _shr_tkfit_start_x_v.back() = tk->Start().X();
        _shr_tkfit_start_y_v.back() = tk->Start().Y();
        _shr_tkfit_start_z_v.back() = tk->Start().Z();
        _shr_tkfit_phi_v.back() = tk->StartDirection().Phi();
        _shr_tkfit_theta_v.back() = tk->StartDirection().Theta();

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
          _shr_tkfit_dedx_v.back()[tkcalo->PlaneID().Plane] = dedx4cm_med;
          _shr_tkfit_dedx_nhits_v.back()[tkcalo->PlaneID().Plane] = dedx4cm.size();
        }
      }
    }
  }
}

void ShowerAnalysis::setBranches(TTree *_tree)
{
  _tree->Branch("shr_dedx_u_v", "std::vector< double >", &_shr_dedx_u_v);
  _tree->Branch("shr_dedx_v_v", "std::vector< double >", &_shr_dedx_v_v);
  _tree->Branch("shr_dedx_y_v", "std::vector< double >", &_shr_dedx_y_v);

  _tree->Branch("shr_energy_u_v", "std::vector< double >", &_shr_energy_u_v);
  _tree->Branch("shr_energy_v_v", "std::vector< double >", &_shr_energy_v_v);
  _tree->Branch("shr_energy_y_v", "std::vector< double >", &_shr_energy_y_v);

  _tree->Branch("shr_pfp_id_v", "std::vector< size_t >", &_shr_pfp_id_v);

  _tree->Branch("shr_start_x_v", "std::vector< double >", &_shr_start_x_v);
  _tree->Branch("shr_start_y_v", "std::vector< double >", &_shr_start_y_v);
  _tree->Branch("shr_start_z_v", "std::vector< double >", &_shr_start_z_v);

  _tree->Branch("shr_dist_v", "std::vector< double >", &_shr_dist_v);


  _tree->Branch("shr_px_v", "std::vector< double >", &_shr_px_v);
  _tree->Branch("shr_py_v", "std::vector< double >", &_shr_py_v);
  _tree->Branch("shr_pz_v", "std::vector< double >", &_shr_pz_v);

  _tree->Branch("shr_theta_v", "std::vector< double >", &_shr_theta_v);
  _tree->Branch("shr_phi_v", "std::vector< double >", &_shr_phi_v);

  _tree->Branch("trkshr_score_v", "std::vector< double >", &_trkshr_score_v);

  _tree->Branch("n_showers", &_n_showers, "n_showers/i");

  _tree->Branch("shr_tkfit_nhits_v", "std::vector< int >", &_shr_tkfit_nhits_v);
  _tree->Branch("shr_tkfit_start_x_v", "std::vector< double >", &_shr_tkfit_start_x_v);
  _tree->Branch("shr_tkfit_start_y_v", "std::vector< double >", &_shr_tkfit_start_y_v);
  _tree->Branch("shr_tkfit_start_z_v", "std::vector< double >", &_shr_tkfit_start_z_v);
  _tree->Branch("shr_tkfit_theta_v", "std::vector< double >", &_shr_tkfit_theta_v);
  _tree->Branch("shr_tkfit_phi_v", "std::vector< double >", &_shr_tkfit_phi_v);
  _tree->Branch("shr_tkfit_dedx_v", "std::vector< std::vector< double > >", &_shr_tkfit_dedx_v);
  _tree->Branch("shr_tkfit_dedx_nhits_v", "std::vector< std::vector< int > >", &_shr_tkfit_dedx_nhits_v);
}

void ShowerAnalysis::resetTTree(TTree *_tree)
{
  _shr_energy_u_v.clear();
  _shr_energy_v_v.clear();
  _shr_energy_y_v.clear();

  _shr_dedx_u_v.clear();
  _shr_dedx_v_v.clear();
  _shr_dedx_y_v.clear();

  _shr_pfp_id_v.clear();

  _shr_start_x_v.clear();
  _shr_start_y_v.clear();
  _shr_start_z_v.clear();

  _shr_theta_v.clear();
  _shr_phi_v.clear();
  _shr_dist_v.clear();

  _shr_px_v.clear();
  _shr_py_v.clear();
  _shr_pz_v.clear();

  _n_showers = 0;
  _trkshr_score_v.clear();

  _shr_tkfit_nhits_v.clear();
  _shr_tkfit_start_x_v.clear();
  _shr_tkfit_start_y_v.clear();
  _shr_tkfit_start_z_v.clear();
  _shr_tkfit_theta_v.clear();
  _shr_tkfit_phi_v.clear();
  _shr_tkfit_dedx_v.clear();
  _shr_tkfit_dedx_nhits_v.clear();
}

DEFINE_ART_CLASS_TOOL(ShowerAnalysis)
} // namespace analysis

#endif
