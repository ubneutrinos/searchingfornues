#ifndef ANALYSIS_SHOWERANALYSIS_CXX
#define ANALYSIS_SHOWERANALYSIS_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "ubobj/CRT/CRTHit.hh"

#include "nusimdata/SimulationBase/MCTruth.h"

// backtracking tools
#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/Geometry.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"
#include "../CommonDefs/ProximityClustering.h"
#include "../CommonDefs/TrackFitterFunctions.h"

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
     * @brief Fill Default info for event associated to neutrino
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
  int _run, _sub, _evt;

  // input variables for tool
  art::InputTag fTRKproducer;
  art::InputTag fCALproducer;
  float fdEdxcmSkip, fdEdxcmLen;
  bool fLocaldEdx;

  std::vector<float> _shr_energy_u_v;
  std::vector<float> _shr_energy_v_v;
  std::vector<float> _shr_energy_y_v;

  std::vector<float> _shr_dedx_u_v;
  std::vector<float> _shr_dedx_v_v;
  std::vector<float> _shr_dedx_y_v;

  std::vector<size_t> _shr_pfp_id_v;

  std::vector<float> _shr_start_x_v;
  std::vector<float> _shr_start_y_v;
  std::vector<float> _shr_start_z_v;

  std::vector<float> _shr_start_U_v;
  std::vector<float> _shr_start_V_v;
  std::vector<float> _shr_dist_v;

  std::vector<float> _shr_px_v;
  std::vector<float> _shr_py_v;
  std::vector<float> _shr_pz_v;

  std::vector<float> _shr_theta_v;
  std::vector<float> _shr_phi_v;

  std::vector<float> _shr_pitch_u_v;
  std::vector<float> _shr_pitch_v_v;
  std::vector<float> _shr_pitch_y_v;

  std::vector<float> _shr_openangle_v;

  std::vector<int> _shr_tkfit_nhits_v;
  std::vector<float> _shr_tkfit_start_x_v;
  std::vector<float> _shr_tkfit_start_y_v;
  std::vector<float> _shr_tkfit_start_z_v;

  std::vector<float> _shr_tkfit_start_U_v;
  std::vector<float> _shr_tkfit_start_V_v;

  std::vector<float> _shr_tkfit_theta_v;
  std::vector<float> _shr_tkfit_phi_v;

  std::vector<float> _shr_tkfit_pitch_u_v;
  std::vector<float> _shr_tkfit_pitch_v_v;
  std::vector<float> _shr_tkfit_pitch_y_v;

  std::vector<float> _shr_tkfit_dedx_u_v;
  std::vector<float> _shr_tkfit_dedx_v_v;
  std::vector<float> _shr_tkfit_dedx_y_v;

  std::vector<float> _shr_tkfit_gap10_dedx_u_v;
  std::vector<float> _shr_tkfit_gap10_dedx_v_v;
  std::vector<float> _shr_tkfit_gap10_dedx_y_v;

  std::vector<int> _shr_tkfit_dedx_nhits_u_v;
  std::vector<int> _shr_tkfit_dedx_nhits_v_v;
  std::vector<int> _shr_tkfit_dedx_nhits_y_v;

  std::vector<float> _shr_moliere_avg_v;
  std::vector<float> _shr_moliere_rms_v;
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
  fdEdxcmSkip = p.get<float>("dEdxcmSkip", 0.0); // how many cm to skip @ vtx for dE/dx calculation
  fdEdxcmLen = p.get<float>("dEdxcmLen", 4.0);   // how long the dE/dx segment should be
  fLocaldEdx = p.get<bool>("LocaldEdx", true);   // use dE/dx from calo?

  // load proximity clustering algorithm
  //PrxyCluster = new searchingfornues::ProximityClustering();
  //PrxyCluster->initialize();
  //PrxyCluster->setRadius(2.0);
  //PrxyCluster->setCellSize(2.0);
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
  _evt = e.event();
  _sub = e.subRun();
  _run = e.run();
  std::cout << "[ShowerAnalysis::analyzeEvent] Run: " << _run << ", SubRun: " << _sub << ", Event: " << _evt << std::endl;
}

void ShowerAnalysis::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
{

  art::InputTag clusproducer("pandora");
  ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, clusproducer, proxy::withAssociated<recob::Hit>(clusproducer));

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
      double xyz[3] = {};

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

    size_t n_shw = slice_pfp_v[i_pfp].get<recob::Shower>().size();
    if (n_shw == 0)
    {
      fillDefault();
    }
    else if (n_shw == 1)
    {
      auto shr = slice_pfp_v[i_pfp].get<recob::Shower>()[0];

      _shr_dedx_u_v.push_back(shr->dEdx()[0]);
      _shr_dedx_v_v.push_back(shr->dEdx()[1]);
      _shr_dedx_y_v.push_back(shr->dEdx()[2]);

      _shr_energy_u_v.push_back(shr->Energy()[0]);
      _shr_energy_v_v.push_back(shr->Energy()[1]);
      _shr_energy_y_v.push_back(shr->Energy()[2]);

      _shr_pfp_id_v.push_back(i_pfp);
      _shr_openangle_v.push_back(shr->OpenAngle());
      _shr_phi_v.push_back(shr->Direction().Phi());
      _shr_theta_v.push_back(shr->Direction().Theta());

      _shr_pitch_u_v.push_back(searchingfornues::getPitch(shr->Direction().Y(), shr->Direction().Z(), 0));
      _shr_pitch_v_v.push_back(searchingfornues::getPitch(shr->Direction().Y(), shr->Direction().Z(), 1));
      _shr_pitch_y_v.push_back(searchingfornues::getPitch(shr->Direction().Y(), shr->Direction().Z(), 2));

      _shr_start_x_v.push_back(shr->ShowerStart().X());
      _shr_start_y_v.push_back(shr->ShowerStart().Y());
      _shr_start_z_v.push_back(shr->ShowerStart().Z());

      _shr_start_U_v.push_back(searchingfornues::YZtoPlanecoordinate(shr->ShowerStart().Y(), shr->ShowerStart().Z(), 0));
      _shr_start_V_v.push_back(searchingfornues::YZtoPlanecoordinate(shr->ShowerStart().Y(), shr->ShowerStart().Z(), 1));

      //fill dummy track fit values, overwrite them later
      _shr_tkfit_nhits_v.push_back(std::numeric_limits<int>::lowest());
      _shr_tkfit_start_x_v.push_back(std::numeric_limits<float>::lowest());
      _shr_tkfit_start_y_v.push_back(std::numeric_limits<float>::lowest());
      _shr_tkfit_start_z_v.push_back(std::numeric_limits<float>::lowest());

      _shr_tkfit_start_U_v.push_back(std::numeric_limits<float>::lowest());
      _shr_tkfit_start_V_v.push_back(std::numeric_limits<float>::lowest());

      _shr_tkfit_phi_v.push_back(std::numeric_limits<float>::lowest());
      _shr_tkfit_theta_v.push_back(std::numeric_limits<float>::lowest());

      _shr_tkfit_pitch_u_v.push_back(std::numeric_limits<float>::lowest());
      _shr_tkfit_pitch_v_v.push_back(std::numeric_limits<float>::lowest());
      _shr_tkfit_pitch_y_v.push_back(std::numeric_limits<float>::lowest());

      _shr_tkfit_dedx_u_v.push_back(std::numeric_limits<float>::lowest());
      _shr_tkfit_dedx_v_v.push_back(std::numeric_limits<float>::lowest());
      _shr_tkfit_dedx_y_v.push_back(std::numeric_limits<float>::lowest());

      _shr_tkfit_gap10_dedx_u_v.push_back(std::numeric_limits<float>::lowest());
      _shr_tkfit_gap10_dedx_v_v.push_back(std::numeric_limits<float>::lowest());
      _shr_tkfit_gap10_dedx_y_v.push_back(std::numeric_limits<float>::lowest());

      _shr_tkfit_dedx_nhits_u_v.push_back(std::numeric_limits<int>::lowest());
      _shr_tkfit_dedx_nhits_v_v.push_back(std::numeric_limits<int>::lowest());
      _shr_tkfit_dedx_nhits_y_v.push_back(std::numeric_limits<int>::lowest());

      _shr_px_v.push_back(shr->Direction().X());
      _shr_py_v.push_back(shr->Direction().Y());
      _shr_pz_v.push_back(shr->Direction().Z());

      _shr_dist_v.push_back((shr->ShowerStart() - nuvtx).Mag());

      float _shrmoliereavg; /**< avg of moliere angle */
      float _shrmoliererms; /**< rms of moliere angle */
      searchingfornues::GetMoliereRadius(slice_pfp_v[i_pfp], _shrmoliereavg, _shrmoliererms);
      _shr_moliere_rms_v.push_back(_shrmoliererms);
      _shr_moliere_avg_v.push_back(_shrmoliereavg);

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

        _shr_tkfit_start_U_v.back() = searchingfornues::YZtoPlanecoordinate(_shr_tkfit_start_y_v.back(), _shr_tkfit_start_z_v.back(), 0);
        _shr_tkfit_start_V_v.back() = searchingfornues::YZtoPlanecoordinate(_shr_tkfit_start_y_v.back(), _shr_tkfit_start_z_v.back(), 1);

        _shr_tkfit_phi_v.back() = tk->StartDirection().Phi();
        _shr_tkfit_theta_v.back() = tk->StartDirection().Theta();

        _shr_tkfit_pitch_u_v.push_back(searchingfornues::getPitch(tk->StartDirection().Y(), tk->StartDirection().Z(), 0));
        _shr_tkfit_pitch_v_v.push_back(searchingfornues::getPitch(tk->StartDirection().Y(), tk->StartDirection().Z(), 1));
        _shr_tkfit_pitch_y_v.push_back(searchingfornues::getPitch(tk->StartDirection().Y(), tk->StartDirection().Z(), 2));

        auto const tkcalos = tk.get<anab::Calorimetry>();

        float calodEdx; // dEdx computed for track-fitter
        int caloNpts;   // number of track-fitter dE/dx hits

        for (const auto &tkcalo : tkcalos)
        {
          // using function from CommonDefs/TrackFitterFunctions.h
          searchingfornues::GetTrackFitdEdx(tkcalo, fdEdxcmSkip, fdEdxcmLen, fLocaldEdx, calodEdx, caloNpts);

          if (tkcalo->PlaneID().Plane == 0)
          {
            _shr_tkfit_dedx_u_v.back() = calodEdx;
            _shr_tkfit_dedx_nhits_u_v.back() = caloNpts;
          }

          if (tkcalo->PlaneID().Plane == 1)
          {
            _shr_tkfit_dedx_v_v.back() = calodEdx;
            _shr_tkfit_dedx_nhits_v_v.back() = caloNpts;
          }

          if (tkcalo->PlaneID().Plane == 2)
          {
            _shr_tkfit_dedx_y_v.back() = calodEdx;
            _shr_tkfit_dedx_nhits_y_v.back() = caloNpts;
          }

          // Gap 1.0 cm
          searchingfornues::GetTrackFitdEdx(tkcalo, 1.0, fdEdxcmLen, fLocaldEdx, calodEdx, caloNpts);
          if (tkcalo->PlaneID().Plane == 2)
          {
            _shr_tkfit_gap10_dedx_y_v.back() = calodEdx;
          }
          else if (tkcalo->PlaneID().Plane == 1)
          {
            _shr_tkfit_gap10_dedx_v_v.back() = calodEdx;
          }
          else if (tkcalo->PlaneID().Plane == 0)
          {
            _shr_tkfit_gap10_dedx_u_v.back() = calodEdx;
          }
        } // for all calorimetry objects
      }
    }
  }
} // analyzeSlice

void ShowerAnalysis::setBranches(TTree *_tree)
{
  _tree->Branch("shr_dedx_u_v", "std::vector< float >", &_shr_dedx_u_v);
  _tree->Branch("shr_dedx_v_v", "std::vector< float >", &_shr_dedx_v_v);
  _tree->Branch("shr_dedx_y_v", "std::vector< float >", &_shr_dedx_y_v);

  _tree->Branch("shr_energy_u_v", "std::vector< float >", &_shr_energy_u_v);
  _tree->Branch("shr_energy_v_v", "std::vector< float >", &_shr_energy_v_v);
  _tree->Branch("shr_energy_y_v", "std::vector< float >", &_shr_energy_y_v);

  _tree->Branch("shr_pfp_id_v", "std::vector< size_t >", &_shr_pfp_id_v);

  _tree->Branch("shr_start_x_v", "std::vector< float >", &_shr_start_x_v);
  _tree->Branch("shr_start_y_v", "std::vector< float >", &_shr_start_y_v);
  _tree->Branch("shr_start_z_v", "std::vector< float >", &_shr_start_z_v);
  _tree->Branch("shr_dist_v", "std::vector< float >", &_shr_dist_v);

  _tree->Branch("shr_start_U_v", "std::vector< float >", &_shr_start_U_v);
  _tree->Branch("shr_start_V_v", "std::vector< float >", &_shr_start_V_v);

  _tree->Branch("shr_px_v", "std::vector< float >", &_shr_px_v);
  _tree->Branch("shr_py_v", "std::vector< float >", &_shr_py_v);
  _tree->Branch("shr_pz_v", "std::vector< float >", &_shr_pz_v);

  _tree->Branch("shr_openangle_v", "std::vector< float >", &_shr_openangle_v);
  _tree->Branch("shr_theta_v", "std::vector< float >", &_shr_theta_v);
  _tree->Branch("shr_phi_v", "std::vector< float >", &_shr_phi_v);

  _tree->Branch("shr_pitch_u_v", "std::vector<float>", &_shr_pitch_u_v);
  _tree->Branch("shr_pitch_v_v", "std::vector<float>", &_shr_pitch_v_v);
  _tree->Branch("shr_pitch_y_v", "std::vector<float>", &_shr_pitch_y_v);

  _tree->Branch("shr_tkfit_nhits_v", "std::vector< int >", &_shr_tkfit_nhits_v);
  _tree->Branch("shr_tkfit_start_x_v", "std::vector< float >", &_shr_tkfit_start_x_v);
  _tree->Branch("shr_tkfit_start_y_v", "std::vector< float >", &_shr_tkfit_start_y_v);
  _tree->Branch("shr_tkfit_start_z_v", "std::vector< float >", &_shr_tkfit_start_z_v);

  _tree->Branch("shr_tkfit_start_U_v", "std::vector< float >", &_shr_tkfit_start_U_v);
  _tree->Branch("shr_tkfit_start_V_v", "std::vector< float >", &_shr_tkfit_start_V_v);

  _tree->Branch("shr_tkfit_theta_v", "std::vector< float >", &_shr_tkfit_theta_v);
  _tree->Branch("shr_tkfit_phi_v", "std::vector< float >", &_shr_tkfit_phi_v);

  _tree->Branch("shr_tkfit_pitch_u_v", "std::vector<float>", &_shr_tkfit_pitch_u_v);
  _tree->Branch("shr_tkfit_pitch_v_v", "std::vector<float>", &_shr_tkfit_pitch_v_v);
  _tree->Branch("shr_tkfit_pitch_y_v", "std::vector<float>", &_shr_tkfit_pitch_y_v);

  _tree->Branch("shr_tkfit_dedx_u_v", "std::vector< float >", &_shr_tkfit_dedx_u_v);
  _tree->Branch("shr_tkfit_dedx_v_v", "std::vector< float >", &_shr_tkfit_dedx_v_v);
  _tree->Branch("shr_tkfit_dedx_y_v", "std::vector< float >", &_shr_tkfit_dedx_y_v);

  _tree->Branch("shr_tkfit_gap10_dedx_u_v", "std::vector< float >", &_shr_tkfit_gap10_dedx_u_v);
  _tree->Branch("shr_tkfit_gap10_dedx_v_v", "std::vector< float >", &_shr_tkfit_gap10_dedx_v_v);
  _tree->Branch("shr_tkfit_gap10_dedx_y_v", "std::vector< float >", &_shr_tkfit_gap10_dedx_y_v);

  _tree->Branch("shr_tkfit_dedx_nhits_u_v", "std::vector< int >", &_shr_tkfit_dedx_nhits_u_v);
  _tree->Branch("shr_tkfit_dedx_nhits_v_v", "std::vector< int >", &_shr_tkfit_dedx_nhits_v_v);
  _tree->Branch("shr_tkfit_dedx_nhits_y_v", "std::vector< int >", &_shr_tkfit_dedx_nhits_y_v);

  _tree->Branch("shr_moliere_avg_v", "std::vector< float >", &_shr_moliere_avg_v);
  _tree->Branch("shr_moliere_rms_v", "std::vector< float >", &_shr_moliere_rms_v);
}

void ShowerAnalysis::fillDefault()
{
  _shr_energy_u_v.push_back(std::numeric_limits<float>::lowest());
  _shr_energy_v_v.push_back(std::numeric_limits<float>::lowest());
  _shr_energy_y_v.push_back(std::numeric_limits<float>::lowest());

  _shr_dedx_u_v.push_back(std::numeric_limits<float>::lowest());
  _shr_dedx_v_v.push_back(std::numeric_limits<float>::lowest());
  _shr_dedx_y_v.push_back(std::numeric_limits<float>::lowest());

  _shr_pfp_id_v.push_back(std::numeric_limits<int>::lowest());

  _shr_start_x_v.push_back(std::numeric_limits<float>::lowest());
  _shr_start_y_v.push_back(std::numeric_limits<float>::lowest());
  _shr_start_z_v.push_back(std::numeric_limits<float>::lowest());

  _shr_start_U_v.push_back(std::numeric_limits<float>::lowest());
  _shr_start_V_v.push_back(std::numeric_limits<float>::lowest());

  _shr_openangle_v.push_back(std::numeric_limits<float>::lowest());
  _shr_theta_v.push_back(std::numeric_limits<float>::lowest());
  _shr_phi_v.push_back(std::numeric_limits<float>::lowest());

  _shr_pitch_u_v.push_back(std::numeric_limits<float>::lowest());
  _shr_pitch_v_v.push_back(std::numeric_limits<float>::lowest());
  _shr_pitch_y_v.push_back(std::numeric_limits<float>::lowest());

  _shr_dist_v.push_back(std::numeric_limits<float>::lowest());

  _shr_px_v.push_back(std::numeric_limits<float>::lowest());
  _shr_py_v.push_back(std::numeric_limits<float>::lowest());
  _shr_pz_v.push_back(std::numeric_limits<float>::lowest());

  _shr_tkfit_nhits_v.push_back(std::numeric_limits<int>::lowest());
  _shr_tkfit_start_x_v.push_back(std::numeric_limits<float>::lowest());
  _shr_tkfit_start_y_v.push_back(std::numeric_limits<float>::lowest());
  _shr_tkfit_start_z_v.push_back(std::numeric_limits<float>::lowest());

  _shr_tkfit_start_U_v.push_back(std::numeric_limits<float>::lowest());
  _shr_tkfit_start_V_v.push_back(std::numeric_limits<float>::lowest());

  _shr_tkfit_theta_v.push_back(std::numeric_limits<float>::lowest());
  _shr_tkfit_phi_v.push_back(std::numeric_limits<float>::lowest());

  _shr_tkfit_pitch_u_v.push_back(std::numeric_limits<float>::lowest());
  _shr_tkfit_pitch_v_v.push_back(std::numeric_limits<float>::lowest());
  _shr_tkfit_pitch_y_v.push_back(std::numeric_limits<float>::lowest());

  _shr_tkfit_dedx_u_v.push_back(std::numeric_limits<float>::lowest());
  _shr_tkfit_dedx_v_v.push_back(std::numeric_limits<float>::lowest());
  _shr_tkfit_dedx_y_v.push_back(std::numeric_limits<float>::lowest());

  _shr_tkfit_gap10_dedx_u_v.push_back(std::numeric_limits<float>::lowest());
  _shr_tkfit_gap10_dedx_v_v.push_back(std::numeric_limits<float>::lowest());
  _shr_tkfit_gap10_dedx_y_v.push_back(std::numeric_limits<float>::lowest());

  _shr_tkfit_dedx_nhits_u_v.push_back(std::numeric_limits<int>::lowest());
  _shr_tkfit_dedx_nhits_v_v.push_back(std::numeric_limits<int>::lowest());
  _shr_tkfit_dedx_nhits_y_v.push_back(std::numeric_limits<int>::lowest());

  _shr_moliere_rms_v.push_back(std::numeric_limits<float>::lowest());
  _shr_moliere_avg_v.push_back(std::numeric_limits<float>::lowest());
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

  _shr_start_U_v.clear();
  _shr_start_V_v.clear();

  _shr_openangle_v.clear();
  _shr_theta_v.clear();
  _shr_phi_v.clear();

  _shr_pitch_u_v.clear();
  _shr_pitch_v_v.clear();
  _shr_pitch_y_v.clear();

  _shr_dist_v.clear();

  _shr_px_v.clear();
  _shr_py_v.clear();
  _shr_pz_v.clear();

  _shr_tkfit_nhits_v.clear();
  _shr_tkfit_start_x_v.clear();
  _shr_tkfit_start_y_v.clear();
  _shr_tkfit_start_z_v.clear();

  _shr_tkfit_start_U_v.clear();
  _shr_tkfit_start_V_v.clear();

  _shr_tkfit_theta_v.clear();
  _shr_tkfit_phi_v.clear();

  _shr_tkfit_pitch_u_v.clear();
  _shr_tkfit_pitch_v_v.clear();
  _shr_tkfit_pitch_y_v.clear();

  _shr_tkfit_dedx_u_v.clear();
  _shr_tkfit_dedx_v_v.clear();
  _shr_tkfit_dedx_y_v.clear();

  _shr_tkfit_gap10_dedx_u_v.clear();
  _shr_tkfit_gap10_dedx_v_v.clear();
  _shr_tkfit_gap10_dedx_y_v.clear();

  _shr_tkfit_dedx_nhits_u_v.clear();
  _shr_tkfit_dedx_nhits_v_v.clear();
  _shr_tkfit_dedx_nhits_y_v.clear();

  _shr_moliere_rms_v.clear();
  _shr_moliere_avg_v.clear();
}

DEFINE_ART_CLASS_TOOL(ShowerAnalysis)
} // namespace analysis

#endif
