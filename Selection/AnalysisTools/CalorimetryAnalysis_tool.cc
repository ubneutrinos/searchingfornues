#ifndef ANALYSIS_CALORIMETRYANALYSIS_CXX
#define ANALYSIS_CALORIMETRYANALYSIS_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include "../CommonDefs/Typedefs.h"

// backtracking tools
#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/Geometry.h"
#include "../CommonDefs/SCECorrections.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"
#include "../CommonDefs/PIDFuncs.h"

#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

namespace analysis
{
////////////////////////////////////////////////////////////////////////
//
// Class:       CalorimetryAnalysis
// File:        CalorimetryAnalysis.cc
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

class CalorimetryAnalysis : public AnalysisToolBase
{

public:
  /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
  CalorimetryAnalysis(const fhicl::ParameterSet &pset);

  /**
     *  @brief  Destructor
     */
  ~CalorimetryAnalysis(){};

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
  trkf::TrackMomentumCalculator _trkmom;

  TParticlePDG *proton = TDatabasePDG::Instance()->GetParticle(2212);
  TParticlePDG *muon = TDatabasePDG::Instance()->GetParticle(13);

  art::InputTag fCALOproducer;
  art::InputTag fPIDproducer;
  art::InputTag fTRKproducer;
  art::InputTag fCLSproducer;
  art::InputTag fBacktrackTag;
  art::InputTag fHproducer;
  art::InputTag fMCRproducer;

  TTree* _calo_tree;

  int _run, _sub, _evt;
  int _isprimary;
  // backtracking information
  int _backtracked_pdg;            // PDG code of backtracked particle
  float _backtracked_e;            // energy of backtracked particle
  float _backtracked_purity;       // purity of backtracking
  float _backtracked_completeness; // completeness of backtracking
  float _backtracked_overlay_purity; // purity of overlay

  float _backtracked_px;
  float _backtracked_py;
  float _backtracked_pz;

  float _backtracked_start_x;
  float _backtracked_start_y;
  float _backtracked_start_z;
  float _backtracked_start_t;
  float _backtracked_start_U;
  float _backtracked_start_V;
  float _backtracked_start_Y;
  float _backtracked_sce_start_x;
  float _backtracked_sce_start_y;
  float _backtracked_sce_start_z;
  float _backtracked_sce_start_U;
  float _backtracked_sce_start_V;
  float _backtracked_sce_start_Y;

  // track information
  int _nplanehits_U;
  int _nplanehits_V;
  int _nplanehits_Y;
  float _trk_score;

  float _trk_start_x;
  float _trk_start_y;
  float _trk_start_z;

  float _trk_theta;
  float _trk_phi;

  float _trk_dir_x;
  float _trk_dir_y;
  float _trk_dir_z;

  float _trk_end_x;
  float _trk_end_y;
  float _trk_end_z;

  float _trk_len;

  float _trk_bragg_p;
  float _trk_bragg_mu;
  float _trk_bragg_mip;
  float _trk_pid_chipr;
  float _trk_pid_chika;
  float _trk_pid_chipi;
  float _trk_pid_chimu;
  float _trk_pida;

  float _trk_bragg_p_u;
  float _trk_bragg_mu_u;
  float _trk_bragg_mip_u;
  float _trk_pid_chipr_u;
  float _trk_pid_chika_u;
  float _trk_pid_chipi_u;
  float _trk_pid_chimu_u;
  float _trk_pida_u;

  float _trk_bragg_p_v;
  float _trk_bragg_mu_v;
  float _trk_bragg_mip_v;
  float _trk_pid_chipr_v;
  float _trk_pid_chika_v;
  float _trk_pid_chipi_v;
  float _trk_pid_chimu_v;
  float _trk_pida_v;

  float _trk_mcs_muon_mom;
  float _trk_energy_proton;
  float _trk_energy_muon;

  // dedx vector
  std::vector<float> _dqdx_u;
  std::vector<float> _dqdx_v;
  std::vector<float> _dqdx_y;

  std::vector<float> _dedx_u;
  std::vector<float> _dedx_v;
  std::vector<float> _dedx_y;

  std::vector<float> _rr_u;
  std::vector<float> _rr_v;
  std::vector<float> _rr_y;

  std::vector<float> _pitch_u;
  std::vector<float> _pitch_v;
  std::vector<float> _pitch_y;
};

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
CalorimetryAnalysis::CalorimetryAnalysis(const fhicl::ParameterSet &p)
{
  fCALOproducer = p.get< art::InputTag > ("CALOproducer");
  fPIDproducer  = p.get< art::InputTag > ("PIDproducer" );
  fTRKproducer  = p.get< art::InputTag > ("TRKproducer" );

  fCLSproducer = p.get<art::InputTag>("CLSproducer");
  fBacktrackTag = p.get<art::InputTag>("BacktrackTag");
  fHproducer = p.get<art::InputTag>("Hproducer");
  fMCRproducer = p.get<art::InputTag>("MCRproducer");

  art::ServiceHandle<art::TFileService> tfs;

  _calo_tree = tfs->make<TTree>("CalorimetryAnalyzer", "Calo Tree");
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void CalorimetryAnalysis::configure(fhicl::ParameterSet const &p)
{
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void CalorimetryAnalysis::analyzeEvent(art::Event const &e, bool fData)
{
  _evt = e.event();
  _sub = e.subRun();
  _run = e.run();
}

void CalorimetryAnalysis::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
{
  // load clusters
  ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, fCLSproducer,
                                                                                        proxy::withAssociated<recob::Hit>(fCLSproducer));

  // load backtrack information
  std::vector<searchingfornues::BtPart> btparts_v;
  std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> assocMCPart;
  if (!fData)
  {
    const std::vector<sim::MCShower> &inputMCShower = *(e.getValidHandle<std::vector<sim::MCShower>>(fMCRproducer));
    const std::vector<sim::MCTrack> &inputMCTrack = *(e.getValidHandle<std::vector<sim::MCTrack>>(fMCRproducer));
    art::ValidHandle<std::vector<recob::Hit>> inputHits = e.getValidHandle<std::vector<recob::Hit>>(fHproducer);
    assocMCPart = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(inputHits, e, fBacktrackTag));
    btparts_v = searchingfornues::initBacktrackingParticleVec(inputMCShower, inputMCTrack, *inputHits, assocMCPart);
  }

  // load calorimetry
  searchingfornues::ProxyCaloColl_t const& calo_proxy = proxy::getCollection<std::vector<recob::Track> >(e, fTRKproducer,
                           proxy::withAssociated<anab::Calorimetry>(fCALOproducer));

  searchingfornues::ProxyPIDColl_t const& pid_proxy = proxy::getCollection<std::vector<recob::Track> >(e, fTRKproducer,
                           proxy::withAssociated<anab::ParticleID>(fPIDproducer));

  for (size_t i_pfp = 0; i_pfp < slice_pfp_v.size(); i_pfp++)
  {
    fillDefault();

    auto pfp = slice_pfp_v[i_pfp];
    if (pfp->IsPrimary())
      continue;

    auto trk_v = pfp.get<recob::Track>();
    if (trk_v.size() != 1)
      continue;
    auto trk = trk_v.at(0);

    //add is primary daughter

    // store pfp information
    _trk_score = searchingfornues::GetTrackShowerScore(pfp);
    // get hits associated to this PFParticle through the clusters
    std::vector<art::Ptr<recob::Hit>> hit_v;
    auto clus_pxy_v = pfp.get<recob::Cluster>();
    for (auto ass_clus : clus_pxy_v)
    {
      // get cluster proxy
      const auto &clus = clus_proxy[ass_clus.key()];
      auto clus_hit_v = clus.get<recob::Hit>();
      auto nhits = clus_hit_v.size();

      if (clus->Plane().Plane == 0)
      {
        _nplanehits_U = nhits;
      }
      else if (clus->Plane().Plane == 1)
      {
        _nplanehits_V = nhits;
      }
      else if (clus->Plane().Plane == 2)
      {
        _nplanehits_Y = nhits;
      }
      for (const auto &hit : clus_hit_v)
      {
        hit_v.push_back(hit);
      }
    } // for all clusters associated to PFP

    // store Backtracking
    if (!fData)
    {
      if (clus_pxy_v.size() != 0)
      {
        float purity = 0., completeness = 0., overlay_purity = 0.;
        int ibt = searchingfornues::getAssocBtPart(hit_v, assocMCPart, btparts_v, purity, completeness, overlay_purity);
        if (ibt >= 0)
        {
          auto &mcp = btparts_v[ibt];
          _backtracked_e = mcp.e;
          _backtracked_pdg = mcp.pdg;
          _backtracked_purity = purity;
          _backtracked_completeness = completeness;
          _backtracked_overlay_purity = overlay_purity;

          _backtracked_px = mcp.px;
          _backtracked_py = mcp.py;
          _backtracked_pz = mcp.pz;
          _backtracked_start_x = mcp.start_x;
          _backtracked_start_y = mcp.start_y;
          _backtracked_start_z = mcp.start_z;
          _backtracked_start_t = mcp.start_t;

          _backtracked_start_U = searchingfornues::YZtoPlanecoordinate(mcp.start_y, mcp.start_z, 0);
          _backtracked_start_V = searchingfornues::YZtoPlanecoordinate(mcp.start_y, mcp.start_z, 1);
          _backtracked_start_Y = searchingfornues::YZtoPlanecoordinate(mcp.start_y, mcp.start_z, 2);

          float reco_st[3] = {mcp.start_x, mcp.start_y, mcp.start_z};

          if (mcp.pdg == 11 || mcp.pdg == 22)
          {
            reco_st[0] += searchingfornues::x_offset(mcp.start_t);
          }
          else
          {
            searchingfornues::True2RecoMappingXYZ(mcp.start_t, mcp.start_x, mcp.start_y, mcp.start_z, reco_st);
          }
          _backtracked_sce_start_x = reco_st[0];
          _backtracked_sce_start_y = reco_st[1];
          _backtracked_sce_start_z = reco_st[2];

          _backtracked_sce_start_U = searchingfornues::YZtoPlanecoordinate(reco_st[1], reco_st[2], 0);
          _backtracked_sce_start_V = searchingfornues::YZtoPlanecoordinate(reco_st[1], reco_st[2], 1);
          _backtracked_sce_start_Y = searchingfornues::YZtoPlanecoordinate(reco_st[1], reco_st[2], 2);
        }
      } // if there are associated clusters
    }

    // get trk proxy in order to fetch PID
    auto trkpxy2 = pid_proxy[trk.key()];
    auto pidpxy_v = trkpxy2.get<anab::ParticleID>();

    //collection plane
    _trk_bragg_p = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 2212, 2),
                              searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 2212, 2));
    _trk_bragg_mu = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 13, 2),
                                searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 13, 2));
    _trk_bragg_mip = searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 0, 2);
    _trk_pid_chipr = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 2212, 2);
    _trk_pid_chimu = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 13, 2);
    _trk_pid_chipi = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 211, 2);
    _trk_pid_chika = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 321, 2);
    _trk_pida_v = searchingfornues::PID(pidpxy_v[0], "PIDA_mean", anab::kPIDA, anab::kForward, 0, 2);

    //u plane
    _trk_bragg_p_u = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 2212, 0),
                              searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 2212, 0));
    _trk_bragg_mu_u = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 13, 0),
                                searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 13, 0));
    _trk_bragg_mip_u = searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 0, 0);
    _trk_pid_chipr_u = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 2212, 0);
    _trk_pid_chimu_u = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 13, 0);
    _trk_pid_chipi_u = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 211, 0);
    _trk_pid_chika_u = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 321, 0);
    _trk_pida_u = searchingfornues::PID(pidpxy_v[0], "PIDA_mean", anab::kPIDA, anab::kForward, 0, 0);

    //v plane
    _trk_bragg_p_v = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 2212, 1),
                              searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 2212, 1));
    _trk_bragg_mu_v = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 13, 1),
                                searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 13, 1));
    _trk_bragg_mip_v = searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 0, 1);
    _trk_pid_chipr_v = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 2212, 1);
    _trk_pid_chimu_v = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 13, 1);
    _trk_pid_chipi_v = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 211, 1);
    _trk_pid_chika_v = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 321, 1);
    _trk_pida_v = searchingfornues::PID(pidpxy_v[0], "PIDA_mean", anab::kPIDA, anab::kForward, 0, 1);


    // Kinetic energy using tabulated stopping power (GeV)
    float mcs_momentum_muon = _trkmom.GetTrackMomentum(trk->Length(), 13);
    float energy_proton = std::sqrt(std::pow(_trkmom.GetTrackMomentum(trk->Length(), 2212), 2) + std::pow(proton->Mass(), 2)) - proton->Mass();
    float energy_muon = std::sqrt(std::pow(mcs_momentum_muon, 2) + std::pow(muon->Mass(), 2)) - muon->Mass();

    _trk_mcs_muon_mom = mcs_momentum_muon;
    _trk_energy_proton = energy_proton;
    _trk_energy_muon = energy_muon;

    _trk_dir_x = trk->StartDirection().X();
    _trk_dir_y = trk->StartDirection().Y();
    _trk_dir_z = trk->StartDirection().Z();

    _trk_start_x = trk->Start().X();
    _trk_start_y = trk->Start().Y();
    _trk_start_z = trk->Start().Z();

    _trk_end_x = trk->End().X();
    _trk_end_y = trk->End().Y();
    _trk_end_z = trk->End().Z();

    _trk_theta = trk->Theta();
    _trk_phi = trk->Phi();

    _trk_len = trk->Length();

    // fill Calorimetry
    auto calo_v = calo_proxy[trk.key()].get<anab::Calorimetry>();
    for (auto const& calo : calo_v)
    {
      auto const& plane = calo->PlaneID().Plane;
      if (plane == 0)
      {
        _dqdx_u = calo->dQdx();
        _dedx_u = calo->dEdx();
	_rr_u = calo->ResidualRange();
	_pitch_u = calo->TrkPitchVec();
      }
      else if (plane == 1)
      {
	_dqdx_v = calo->dQdx();
        _dedx_v = calo->dEdx();
	_rr_v = calo->ResidualRange();
	_pitch_v = calo->TrkPitchVec();
      }
      else if (plane == 2) //collection
      {
	_dqdx_y = calo->dQdx();
        _dedx_y = calo->dEdx();
	_rr_y = calo->ResidualRange();
	_pitch_y = calo->TrkPitchVec();
      }
    }
    _calo_tree->Fill();
  } // for all PFParticles
}

void CalorimetryAnalysis::fillDefault()
{
  _isprimary = std::numeric_limits<int>::lowest();

  // backtracking information
  _backtracked_pdg = std::numeric_limits<int>::lowest();            // PDG code of backtracked particle
  _backtracked_e = std::numeric_limits<float>::lowest();            // energy of backtracked particle
  _backtracked_purity = std::numeric_limits<float>::lowest();       // purity of backtracking
  _backtracked_completeness = std::numeric_limits<float>::lowest(); // completeness of backtracking
  _backtracked_overlay_purity = std::numeric_limits<float>::lowest(); // purity of overlay

  _backtracked_px = std::numeric_limits<float>::lowest();
  _backtracked_py = std::numeric_limits<float>::lowest();
  _backtracked_pz = std::numeric_limits<float>::lowest();

  _backtracked_start_x = std::numeric_limits<float>::lowest();
  _backtracked_start_y = std::numeric_limits<float>::lowest();
  _backtracked_start_z = std::numeric_limits<float>::lowest();
  _backtracked_start_t = std::numeric_limits<float>::lowest();
  _backtracked_start_U = std::numeric_limits<float>::lowest();
  _backtracked_start_V = std::numeric_limits<float>::lowest();
  _backtracked_start_Y = std::numeric_limits<float>::lowest();
  _backtracked_sce_start_x = std::numeric_limits<float>::lowest();
  _backtracked_sce_start_y = std::numeric_limits<float>::lowest();
  _backtracked_sce_start_z = std::numeric_limits<float>::lowest();
  _backtracked_sce_start_U = std::numeric_limits<float>::lowest();
  _backtracked_sce_start_V = std::numeric_limits<float>::lowest();
  _backtracked_sce_start_Y = std::numeric_limits<float>::lowest();

  // track information
  _nplanehits_U = std::numeric_limits<int>::lowest();
  _nplanehits_V = std::numeric_limits<int>::lowest();
  _nplanehits_Y = std::numeric_limits<int>::lowest();
  _trk_score = std::numeric_limits<float>::lowest();

  _trk_start_x = std::numeric_limits<float>::lowest();
  _trk_start_y = std::numeric_limits<float>::lowest();
  _trk_start_z = std::numeric_limits<float>::lowest();

  _trk_theta = std::numeric_limits<float>::lowest();
  _trk_phi = std::numeric_limits<float>::lowest();

  _trk_dir_x = std::numeric_limits<float>::lowest();
  _trk_dir_y = std::numeric_limits<float>::lowest();
  _trk_dir_z = std::numeric_limits<float>::lowest();

  _trk_end_x = std::numeric_limits<float>::lowest();
  _trk_end_y = std::numeric_limits<float>::lowest();
  _trk_end_z = std::numeric_limits<float>::lowest();

  _trk_len = std::numeric_limits<float>::lowest();

  _trk_bragg_p = std::numeric_limits<float>::lowest();
  _trk_bragg_mu = std::numeric_limits<float>::lowest();
  _trk_bragg_mip = std::numeric_limits<float>::lowest();
  _trk_pid_chipr = std::numeric_limits<float>::lowest();
  _trk_pid_chika = std::numeric_limits<float>::lowest();
  _trk_pid_chipi = std::numeric_limits<float>::lowest();
  _trk_pid_chimu = std::numeric_limits<float>::lowest();
  _trk_pida = std::numeric_limits<float>::lowest();

  _trk_bragg_p_u = std::numeric_limits<float>::lowest();
  _trk_bragg_mu_u = std::numeric_limits<float>::lowest();
  _trk_bragg_mip_u = std::numeric_limits<float>::lowest();
  _trk_pid_chipr_u = std::numeric_limits<float>::lowest();
  _trk_pid_chika_u = std::numeric_limits<float>::lowest();
  _trk_pid_chipi_u = std::numeric_limits<float>::lowest();
  _trk_pid_chimu_u = std::numeric_limits<float>::lowest();
  _trk_pida_u = std::numeric_limits<float>::lowest();

  _trk_bragg_p_v = std::numeric_limits<float>::lowest();
  _trk_bragg_mu_v = std::numeric_limits<float>::lowest();
  _trk_bragg_mip_v = std::numeric_limits<float>::lowest();
  _trk_pid_chipr_v = std::numeric_limits<float>::lowest();
  _trk_pid_chika_v = std::numeric_limits<float>::lowest();
  _trk_pid_chipi_v = std::numeric_limits<float>::lowest();
  _trk_pid_chimu_v = std::numeric_limits<float>::lowest();
  _trk_pida_v = std::numeric_limits<float>::lowest();

  _trk_mcs_muon_mom = std::numeric_limits<float>::lowest();
  _trk_energy_proton = std::numeric_limits<float>::lowest();
  _trk_energy_muon = std::numeric_limits<float>::lowest();

  // dedx vector
  _dqdx_u.clear();
  _dqdx_v.clear();
  _dqdx_y.clear();

  _dedx_u.clear();
  _dedx_v.clear();
  _dedx_y.clear();

  _rr_u.clear();
  _rr_v.clear();
  _rr_y.clear();

  _pitch_u.clear();
  _pitch_v.clear();
  _pitch_y.clear();
}

void CalorimetryAnalysis::setBranches(TTree *_tree)
{
  _calo_tree->Branch("run", &_run, "run/i");
  _calo_tree->Branch("sub", &_sub, "sub/i");
  _calo_tree->Branch("evt", &_evt, "evt/i");
  _calo_tree->Branch("isprimary", &_isprimary, "isprimary/i");
  // backtracking information
  _calo_tree->Branch("backtracked_pdg", &_backtracked_pdg, "backtracked_pdg/i");            // PDG code of backtracked particle
  _calo_tree->Branch("backtracked_e", &_backtracked_e, "backtracked_e/f");            // energy of backtracked particle
  _calo_tree->Branch("backtracked_purity", &_backtracked_purity, "backtracked_purity/f");       // purity of backtracking
  _calo_tree->Branch("backtracked_completeness", &_backtracked_completeness, "backtracked_completeness/f"); // completeness of backtracking
  _calo_tree->Branch("backtracked_overlay_purity", &_backtracked_overlay_purity, "backtracked_overlay_purity/f"); // purity of overlay

  _calo_tree->Branch("backtracked_px", &_backtracked_px, "backtracked_px/f");
  _calo_tree->Branch("backtracked_py", &_backtracked_py, "backtracked_py/f");
  _calo_tree->Branch("backtracked_pz", &_backtracked_pz, "backtracked_pz/f");

  _calo_tree->Branch("backtracked_start_x", &_backtracked_start_x, "backtracked_start_x/f");
  _calo_tree->Branch("backtracked_start_y", &_backtracked_start_y, "backtracked_start_y/f");
  _calo_tree->Branch("backtracked_start_z", &_backtracked_start_z, "backtracked_start_z/f");
  _calo_tree->Branch("backtracked_start_t", &_backtracked_start_t, "backtracked_start_t/f");
  _calo_tree->Branch("backtracked_start_U", &_backtracked_start_U, "backtracked_start_U/f");
  _calo_tree->Branch("backtracked_start_V", &_backtracked_start_V, "backtracked_start_V/f");
  _calo_tree->Branch("backtracked_start_Y", &_backtracked_start_Y, "backtracked_start_Y/f");
  _calo_tree->Branch("backtracked_sce_start_x", &_backtracked_sce_start_x, "backtracked_sce_start_x/f");
  _calo_tree->Branch("backtracked_sce_start_y", &_backtracked_sce_start_y, "backtracked_sce_start_y/f");
  _calo_tree->Branch("backtracked_sce_start_z", &_backtracked_sce_start_z, "backtracked_sce_start_z/f");
  _calo_tree->Branch("backtracked_sce_start_U", &_backtracked_sce_start_U, "backtracked_sce_start_U/f");
  _calo_tree->Branch("backtracked_sce_start_V", &_backtracked_sce_start_V, "backtracked_sce_start_V/f");
  _calo_tree->Branch("backtracked_sce_start_Y", &_backtracked_sce_start_Y, "backtracked_sce_start_Y/f");

  // track information
  _calo_tree->Branch("nplanehits_U", &_nplanehits_U, "nplanehits_U/i");
  _calo_tree->Branch("nplanehits_V", &_nplanehits_V, "nplanehits_V/i");
  _calo_tree->Branch("nplanehits_Y", &_nplanehits_Y, "nplanehits_Y/i");
  _calo_tree->Branch("trk_score", &_trk_score, "trk_score/f");

  _calo_tree->Branch("trk_start_x", &_trk_start_x, "trk_start_x/f");
  _calo_tree->Branch("trk_start_y", &_trk_start_y, "trk_start_y/f");
  _calo_tree->Branch("trk_start_z", &_trk_start_z, "trk_start_z/f");

  _calo_tree->Branch("trk_theta", &_trk_theta, "trk_theta/f");
  _calo_tree->Branch("trk_phi", &_trk_phi, "trk_phi/f");

  _calo_tree->Branch("trk_dir_x", &_trk_dir_x, "trk_dir_x/f");
  _calo_tree->Branch("trk_dir_y", &_trk_dir_y, "trk_dir_y/f");
  _calo_tree->Branch("trk_dir_z", &_trk_dir_z, "trk_dir_z/f");

  _calo_tree->Branch("trk_end_x", &_trk_end_x, "trk_end_x/f");
  _calo_tree->Branch("trk_end_y", &_trk_end_y, "trk_end_y/f");
  _calo_tree->Branch("trk_end_z", &_trk_end_z, "trk_end_z/f");

  _calo_tree->Branch("trk_len", &_trk_len, "trk_len/f");

  _calo_tree->Branch("trk_bragg_p", &_trk_bragg_p, "trk_bragg_p/f");
  _calo_tree->Branch("trk_bragg_mu", &_trk_bragg_mu, "trk_bragg_mu/f");
  _calo_tree->Branch("trk_bragg_mip", &_trk_bragg_mip, "trk_bragg_mip/f");
  _calo_tree->Branch("trk_pid_chipr", &_trk_pid_chipr, "trk_pid_chipr/f");
  _calo_tree->Branch("trk_pid_chika", &_trk_pid_chika, "trk_pid_chika/f");
  _calo_tree->Branch("trk_pid_chipi", &_trk_pid_chipi, "trk_pid_chipi/f");
  _calo_tree->Branch("trk_pid_chimu", &_trk_pid_chimu, "trk_pid_chimu/f");
  _calo_tree->Branch("trk_pida", &_trk_pida, "trk_pida/f");

  _calo_tree->Branch("trk_bragg_p_u", &_trk_bragg_p_u, "trk_bragg_p_u/f");
  _calo_tree->Branch("trk_bragg_mu_u", &_trk_bragg_mu_u, "trk_bragg_mu_u/f");
  _calo_tree->Branch("trk_bragg_mip_u", &_trk_bragg_mip_u, "trk_bragg_mip_u/f");
  _calo_tree->Branch("trk_pid_chipr_u", &_trk_pid_chipr_u, "trk_pid_chipr_u/f");
  _calo_tree->Branch("trk_pid_chika_u", &_trk_pid_chika_u, "trk_pid_chika_u/f");
  _calo_tree->Branch("trk_pid_chipi_u", &_trk_pid_chipi_u, "trk_pid_chipi_u/f");
  _calo_tree->Branch("trk_pid_chimu_u", &_trk_pid_chimu_u, "trk_pid_chimu_u/f");
  _calo_tree->Branch("trk_pida_u", &_trk_pida_u, "trk_pida_u/f");

  _calo_tree->Branch("trk_bragg_p_v", &_trk_bragg_p_v, "trk_bragg_p_v/f");
  _calo_tree->Branch("trk_bragg_mu_v", &_trk_bragg_mu_v, "trk_bragg_mu_v/f");
  _calo_tree->Branch("trk_bragg_mip_v", &_trk_bragg_mip_v, "trk_bragg_mip_v/f");
  _calo_tree->Branch("trk_pid_chipr_v", &_trk_pid_chipr_v, "trk_pid_chipr_v/f");
  _calo_tree->Branch("trk_pid_chika_v", &_trk_pid_chika_v, "trk_pid_chika_v/f");
  _calo_tree->Branch("trk_pid_chipi_v", &_trk_pid_chipi_v, "trk_pid_chipi_v/f");
  _calo_tree->Branch("trk_pid_chimu_v", &_trk_pid_chimu_v, "trk_pid_chimu_v/f");
  _calo_tree->Branch("trk_pida_v", &_trk_pida_v, "trk_pida_v/f");

  _calo_tree->Branch("trk_mcs_muon_mom", &_trk_mcs_muon_mom, "trk_mcs_muon_mom/f");
  _calo_tree->Branch("trk_energy_proton", &_trk_energy_proton, "trk_energy_proton/f");
  _calo_tree->Branch("trk_energy_muon", &_trk_energy_muon, "trk_energy_muon/f");

  // dedx vector
  _calo_tree->Branch("dqdx_u", "std::vector<float>", &_dqdx_u);
  _calo_tree->Branch("dqdx_v", "std::vector<float>", &_dqdx_v);
  _calo_tree->Branch("dqdx_y", "std::vector<float>", &_dqdx_y);

  _calo_tree->Branch("dedx_u", "std::vector<float>", &_dedx_u);
  _calo_tree->Branch("dedx_v", "std::vector<float>", &_dedx_v);
  _calo_tree->Branch("dedx_y", "std::vector<float>", &_dedx_y);

  _calo_tree->Branch("rr_u", "std::vector<float>", &_rr_u);
  _calo_tree->Branch("rr_v", "std::vector<float>", &_rr_v);
  _calo_tree->Branch("rr_y", "std::vector<float>", &_rr_y);

  _calo_tree->Branch("pitch_u", "std::vector<float>", &_pitch_u);
  _calo_tree->Branch("pitch_v", "std::vector<float>", &_pitch_v);
  _calo_tree->Branch("pitch_y", "std::vector<float>", &_pitch_y);
}

void CalorimetryAnalysis::resetTTree(TTree *_tree)
{
}

DEFINE_ART_CLASS_TOOL(CalorimetryAnalysis)
} // namespace analysis

#endif
