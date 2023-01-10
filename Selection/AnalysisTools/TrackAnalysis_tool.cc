#ifndef ANALYSIS_TRACKANALYSIS_CXX
#define ANALYSIS_TRACKANALYSIS_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include "../CommonDefs/Typedefs.h"

// backtracking tools
#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"
#include "../CommonDefs/PIDFuncs.h"
#include "../CommonDefs/SCECorrections.h"
#include "../CommonDefs/Geometry.h"

#include "../CommonDefs/LLR_PID.h"
#include "../CommonDefs/LLRPID_proton_muon_lookup.h"
#include "../CommonDefs/LLRPID_correction_lookup.h"
#include "../CommonDefs/CalibrationFuncs.h"

#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

namespace analysis
{
////////////////////////////////////////////////////////////////////////
//
// Class:       TrackAnalysis
// File:        TrackAnalysis.cc
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

class TrackAnalysis : public AnalysisToolBase
{

public:
  /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
  TrackAnalysis(const fhicl::ParameterSet &pset);

  /**
     *  @brief  Destructor
     */
  ~TrackAnalysis(){};

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

  float CalculateTrackTrunkdEdx(const std::vector<float> &dEdx_values);

  const trkf::TrackMomentumCalculator _trkmom;
  const trkf::TrajectoryMCSFitter _mcsfitter;

  TParticlePDG *proton = TDatabasePDG::Instance()->GetParticle(2212);
  TParticlePDG *muon = TDatabasePDG::Instance()->GetParticle(13);

  searchingfornues::LLRPID llr_pid_calculator;
  searchingfornues::ProtonMuonLookUpParameters protonmuon_parameters;
  searchingfornues::CorrectionLookUpParameters correction_parameters;

  art::InputTag fCALOproducer;
  art::InputTag fPIDproducer;
  art::InputTag fTRKproducer;
  art::InputTag fBacktrackTag;
  art::InputTag fHproducer;
  bool fRecalibrateHits;
  float fEnergyThresholdForMCHits;
  std::vector<float> fADCtoE; // vector of ADC to # of e- conversion [to be taken from production reco2 fhicl files]

  int _run, _sub, _evt;

  std::vector<size_t> _trk_pfp_id_v;

  std::vector<float> _trk_start_x_v;
  std::vector<float> _trk_start_y_v;
  std::vector<float> _trk_start_z_v;

  std::vector<float> _trk_sce_start_x_v;
  std::vector<float> _trk_sce_start_y_v;
  std::vector<float> _trk_sce_start_z_v;

  std::vector<float> _trk_distance_v;

  std::vector<float> _trk_theta_v;
  std::vector<float> _trk_phi_v;

  std::vector<float> _trk_dir_x_v;
  std::vector<float> _trk_dir_y_v;
  std::vector<float> _trk_dir_z_v;

  std::vector<float> _trk_end_x_v;
  std::vector<float> _trk_end_y_v;
  std::vector<float> _trk_end_z_v;

  std::vector<float> _trk_sce_end_x_v;
  std::vector<float> _trk_sce_end_y_v;
  std::vector<float> _trk_sce_end_z_v;

  std::vector<float> _trk_len_v;

  std::vector<float> _trk_bragg_p_v;
  std::vector<float> _trk_bragg_mu_v;
  std::vector<float> _trk_bragg_pion_v;
  std::vector<float> _trk_bragg_mip_v;
  std::vector<float> _trk_pid_chipr_v;
  std::vector<float> _trk_pid_chika_v;
  std::vector<float> _trk_pid_chipi_v;
  std::vector<float> _trk_pid_chimu_v;
  std::vector<float> _trk_pida_v;

  std::vector<float> _trk_bragg_p_u_v;
  std::vector<float> _trk_bragg_mu_u_v;
  std::vector<float> _trk_bragg_pion_u_v;
  std::vector<float> _trk_bragg_mip_u_v;
  std::vector<float> _trk_pid_chipr_u_v;
  std::vector<float> _trk_pid_chika_u_v;
  std::vector<float> _trk_pid_chipi_u_v;
  std::vector<float> _trk_pid_chimu_u_v;
  std::vector<float> _trk_pida_u_v;

  std::vector<float> _trk_bragg_p_v_v;
  std::vector<float> _trk_bragg_mu_v_v;
  std::vector<float> _trk_bragg_pion_v_v;
  std::vector<float> _trk_bragg_mip_v_v;
  std::vector<float> _trk_pid_chipr_v_v;
  std::vector<float> _trk_pid_chika_v_v;
  std::vector<float> _trk_pid_chipi_v_v;
  std::vector<float> _trk_pid_chimu_v_v;
  std::vector<float> _trk_pida_v_v;

  std::vector<float> _trk_llr_pid_u_v;
  std::vector<float> _trk_llr_pid_v_v;
  std::vector<float> _trk_llr_pid_y_v;
  std::vector<float> _trk_llr_pid_v;

  std::vector<float> _trk_llr_pid_score_v;

  std::vector<float> _trk_mcs_muon_mom_v;
  std::vector<float> _trk_range_muon_mom_v;
  std::vector<float> _trk_energy_proton_v;
  std::vector<float> _trk_energy_muon_v;
  std::vector<float> _trk_calo_energy_u_v;
  std::vector<float> _trk_calo_energy_v_v;
  std::vector<float> _trk_calo_energy_y_v;

  std::vector<float> _trk_trunk_dEdx_u_v;
  std::vector<float> _trk_trunk_dEdx_v_v;
  std::vector<float> _trk_trunk_dEdx_y_v;
};

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
TrackAnalysis::TrackAnalysis(const fhicl::ParameterSet &p) : _mcsfitter(fhicl::Table<trkf::TrajectoryMCSFitter::Config>(p.get<fhicl::ParameterSet>("mcsfitmu")))
{
  fCALOproducer = p.get<art::InputTag>("CALOproducer");
  fPIDproducer = p.get<art::InputTag>("PIDproducer");
  fTRKproducer = p.get<art::InputTag>("TRKproducer");
  fBacktrackTag = p.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch");
  fHproducer = p.get<art::InputTag>("Hproducer", "gaushit");
  fEnergyThresholdForMCHits = p.get<float>("EnergyThresholdForMCHits", 0.1);
  fRecalibrateHits = p.get<bool>("RecalibrateHits", false);
  fADCtoE = p.get<std::vector<float>>("ADCtoE");

  // set dedx pdf parameters
  llr_pid_calculator.set_dedx_binning(0, protonmuon_parameters.dedx_edges_pl_0);
  llr_pid_calculator.set_par_binning(0, protonmuon_parameters.parameters_edges_pl_0);
  llr_pid_calculator.set_lookup_tables(0, protonmuon_parameters.dedx_pdf_pl_0);

  llr_pid_calculator.set_dedx_binning(1, protonmuon_parameters.dedx_edges_pl_1);
  llr_pid_calculator.set_par_binning(1, protonmuon_parameters.parameters_edges_pl_1);
  llr_pid_calculator.set_lookup_tables(1, protonmuon_parameters.dedx_pdf_pl_1);

  llr_pid_calculator.set_dedx_binning(2, protonmuon_parameters.dedx_edges_pl_2);
  llr_pid_calculator.set_par_binning(2, protonmuon_parameters.parameters_edges_pl_2);
  llr_pid_calculator.set_lookup_tables(2, protonmuon_parameters.dedx_pdf_pl_2);

  // set correction parameters
  if (fRecalibrateHits)
  {
    llr_pid_calculator.set_corr_par_binning(0, correction_parameters.parameter_correction_edges_pl_0);
    llr_pid_calculator.set_correction_tables(0, correction_parameters.correction_table_pl_0);

    llr_pid_calculator.set_corr_par_binning(1, correction_parameters.parameter_correction_edges_pl_1);
    llr_pid_calculator.set_correction_tables(1, correction_parameters.correction_table_pl_1);

    llr_pid_calculator.set_corr_par_binning(2, correction_parameters.parameter_correction_edges_pl_2);
    llr_pid_calculator.set_correction_tables(2, correction_parameters.correction_table_pl_2);
  }
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void TrackAnalysis::configure(fhicl::ParameterSet const &p)
{
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void TrackAnalysis::analyzeEvent(art::Event const &e, bool fData)
{
  _evt = e.event();
  _sub = e.subRun();
  _run = e.run();
  std::cout << "[TrackAnalysis::analyzeEvent] Run: " << _run << ", SubRun: " << _sub << ", Event: " << _evt << std::endl;
}

void TrackAnalysis::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
{
  searchingfornues::ProxyCaloColl_t const &calo_proxy = proxy::getCollection<std::vector<recob::Track>>(e, fTRKproducer,
                                                                                                        proxy::withAssociated<anab::Calorimetry>(fCALOproducer));

  searchingfornues::ProxyPIDColl_t const &pid_proxy = proxy::getCollection<std::vector<recob::Track>>(e, fTRKproducer,
                                                                                                      proxy::withAssociated<anab::ParticleID>(fPIDproducer));

  TVector3 nuvtx;
  for (auto pfp : slice_pfp_v)
  {
    if (pfp->IsPrimary())
    {
      double xyz[3] = {};
      auto vtx = pfp.get<recob::Vertex>();
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

  // grab hit backtracked information to be able to apply hit-by-hit re-calibrations
  // only to MC charge
  std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> assocMCPart;

  if (!fData)
  {
    art::ValidHandle<std::vector<recob::Hit>> inputHits = e.getValidHandle<std::vector<recob::Hit>>(fHproducer);
    assocMCPart = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(inputHits, e, fBacktrackTag));
  }

  for (size_t i_pfp = 0; i_pfp < slice_pfp_v.size(); i_pfp++)
  {
    auto pfp = slice_pfp_v[i_pfp];
    if (pfp->IsPrimary())
      continue;

    auto trk_v = pfp.get<recob::Track>();

    if (trk_v.size() == 1)
    {
      auto trk = trk_v.at(0);

      // get trk proxy in order to fetch PID
      auto trkpxy2 = pid_proxy[trk.key()];
      auto pidpxy_v = trkpxy2.get<anab::ParticleID>();

      //collection plane
      float bragg_p = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 2212, 2),
                               searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 2212, 2));

      float bragg_mu = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 13, 2),
                                searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 13, 2));

      float bragg_pion = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 211, 2),
                                  searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 211, 2));

      float bragg_mip = searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 0, 2);

      float pidchipr = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 2212, 2);
      float pidchimu = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 13, 2);
      float pidchipi = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 211, 2);
      float pidchika = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 321, 2);

      float pida_mean = searchingfornues::PID(pidpxy_v[0], "PIDA_mean", anab::kPIDA, anab::kForward, 0, 2);

      _trk_bragg_p_v.push_back(bragg_p);
      _trk_bragg_mu_v.push_back(bragg_mu);
      _trk_bragg_pion_v.push_back(bragg_pion);
      _trk_bragg_mip_v.push_back(bragg_mip);
      _trk_pid_chipr_v.push_back(pidchipr);
      _trk_pid_chimu_v.push_back(pidchimu);
      _trk_pid_chipi_v.push_back(pidchipi);
      _trk_pid_chika_v.push_back(pidchika);
      _trk_pida_v.push_back(pida_mean);

      //u plane
      float bragg_p_u = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 2212, 0),
                                 searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 2212, 0));

      float bragg_mu_u = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 13, 0),
                                  searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 13, 0));

      float bragg_pion_u = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 211, 0),
                                    searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 211, 0));

      float bragg_mip_u = searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 0, 0);

      float pidchipr_u = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 2212, 0);
      float pidchimu_u = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 13, 0);
      float pidchipi_u = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 211, 0);
      float pidchika_u = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 321, 0);

      float pida_mean_u = searchingfornues::PID(pidpxy_v[0], "PIDA_mean", anab::kPIDA, anab::kForward, 0, 0);

      _trk_bragg_p_u_v.push_back(bragg_p_u);
      _trk_bragg_mu_u_v.push_back(bragg_mu_u);
      _trk_bragg_pion_u_v.push_back(bragg_pion_u);
      _trk_bragg_mip_u_v.push_back(bragg_mip_u);
      _trk_pid_chipr_u_v.push_back(pidchipr_u);
      _trk_pid_chimu_u_v.push_back(pidchimu_u);
      _trk_pid_chipi_u_v.push_back(pidchipi_u);
      _trk_pid_chika_u_v.push_back(pidchika_u);
      _trk_pida_u_v.push_back(pida_mean_u);

      //v plane
      float bragg_p_v = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 2212, 1),
                                 searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 2212, 1));

      float bragg_mu_v = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 13, 1),
                                  searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 13, 1));

      float bragg_pion_v = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 211, 1),
                                    searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 211, 1));

      float bragg_mip_v = searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 0, 1);

      float pidchipr_v = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 2212, 1);
      float pidchimu_v = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 13, 1);
      float pidchipi_v = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 211, 1);
      float pidchika_v = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 321, 1);

      float pida_mean_v = searchingfornues::PID(pidpxy_v[0], "PIDA_mean", anab::kPIDA, anab::kForward, 0, 1);

      _trk_bragg_p_v_v.push_back(bragg_p_v);
      _trk_bragg_mu_v_v.push_back(bragg_mu_v);
      _trk_bragg_pion_v_v.push_back(bragg_pion_v);
      _trk_bragg_mip_v_v.push_back(bragg_mip_v);
      _trk_pid_chipr_v_v.push_back(pidchipr_v);
      _trk_pid_chimu_v_v.push_back(pidchimu_v);
      _trk_pid_chipi_v_v.push_back(pidchipi_v);
      _trk_pid_chika_v_v.push_back(pidchika_v);
      _trk_pida_v_v.push_back(pida_mean_v);

      // Kinetic energy using tabulated stopping power (GeV)
      float mcs_momentum_muon = _mcsfitter.fitMcs(trk->Trajectory(), 13).bestMomentum();
      float range_momentum_muon = _trkmom.GetTrackMomentum(searchingfornues::GetSCECorrTrackLength(trk), 13);
      float energy_proton = std::sqrt(std::pow(_trkmom.GetTrackMomentum(searchingfornues::GetSCECorrTrackLength(trk), 2212), 2) + std::pow(proton->Mass(), 2)) - proton->Mass();
      float energy_muon = std::sqrt(std::pow(mcs_momentum_muon, 2) + std::pow(muon->Mass(), 2)) - muon->Mass();

      _trk_mcs_muon_mom_v.push_back(mcs_momentum_muon);
      _trk_range_muon_mom_v.push_back(range_momentum_muon);
      _trk_energy_proton_v.push_back(energy_proton);
      _trk_energy_muon_v.push_back(energy_muon);
      _trk_calo_energy_u_v.push_back(-1);
      _trk_calo_energy_v_v.push_back(-1);
      _trk_calo_energy_y_v.push_back(-1);

      _trk_dir_x_v.push_back(trk->StartDirection().X());
      _trk_dir_y_v.push_back(trk->StartDirection().Y());
      _trk_dir_z_v.push_back(trk->StartDirection().Z());

      _trk_start_x_v.push_back(trk->Start().X());
      _trk_start_y_v.push_back(trk->Start().Y());
      _trk_start_z_v.push_back(trk->Start().Z());

      float _trk_start_sce[3];
      searchingfornues::ApplySCECorrectionXYZ(trk->Start().X(), trk->Start().Y(), trk->Start().Z(), _trk_start_sce);
      _trk_sce_start_x_v.push_back(_trk_start_sce[0]);
      _trk_sce_start_y_v.push_back(_trk_start_sce[1]);
      _trk_sce_start_z_v.push_back(_trk_start_sce[2]);

      _trk_end_x_v.push_back(trk->End().X());
      _trk_end_y_v.push_back(trk->End().Y());
      _trk_end_z_v.push_back(trk->End().Z());

      float _trk_end_sce[3];
      searchingfornues::ApplySCECorrectionXYZ(trk->End().X(), trk->End().Y(), trk->End().Z(), _trk_end_sce);
      _trk_sce_end_x_v.push_back(_trk_end_sce[0]);
      _trk_sce_end_y_v.push_back(_trk_end_sce[1]);
      _trk_sce_end_z_v.push_back(_trk_end_sce[2]);

      _trk_theta_v.push_back(trk->Theta());
      _trk_phi_v.push_back(trk->Phi());

      _trk_len_v.push_back(searchingfornues::GetSCECorrTrackLength(trk));

      TVector3 trk_vtx_v;
      trk_vtx_v.SetXYZ(trk->Start().X(), trk->Start().Y(), trk->Start().Z());
      trk_vtx_v -= nuvtx;
      _trk_distance_v.push_back(trk_vtx_v.Mag());

      _trk_pfp_id_v.push_back(i_pfp);

      //PID LLR calculator
      _trk_llr_pid_u_v.push_back(0);
      _trk_llr_pid_v_v.push_back(0);
      _trk_llr_pid_y_v.push_back(0);
      _trk_llr_pid_v.push_back(0);
      _trk_llr_pid_score_v.push_back(0);

      // track trunk dEdx
      _trk_trunk_dEdx_u_v.push_back(std::numeric_limits<float>::lowest());
      _trk_trunk_dEdx_v_v.push_back(std::numeric_limits<float>::lowest());
      _trk_trunk_dEdx_y_v.push_back(std::numeric_limits<float>::lowest());

      auto calo_v = calo_proxy[trk.key()].get<anab::Calorimetry>();
      for (auto const &calo : calo_v)
      {
        auto const &plane = calo->PlaneID().Plane;
        auto const &dqdx_values = calo->dQdx();
        auto const &dedx_values = calo->dEdx();
        auto const &rr = calo->ResidualRange();
        auto const &pitch = calo->TrkPitchVec();
        auto const& xyz_v = calo->XYZ();
        std::vector<std::vector<float>> par_values;
        par_values.push_back(rr);
        par_values.push_back(pitch);

        float calo_energy = 0;
        std::vector<float> dqdx_values_corrected, dedx_values_corrected;


        if (fData || !fRecalibrateHits)
        {
          dqdx_values_corrected = dqdx_values;
        }
        else
        {
          dqdx_values_corrected = llr_pid_calculator.correct_many_hits_one_plane(calo, trk.value(), assocMCPart, fRecalibrateHits, fEnergyThresholdForMCHits, false);
        }

        for (size_t i = 0; i < dqdx_values_corrected.size(); i++)
        {
          float aux_dedx;
          aux_dedx = searchingfornues::ModBoxCorrection(dqdx_values_corrected[i]*fADCtoE[plane], xyz_v[i].X(), xyz_v[i].Y(), xyz_v[i].Z());
          dedx_values_corrected.push_back(aux_dedx);
          calo_energy += aux_dedx * pitch[i];
        }

        float trk_trunk_dEdx = CalculateTrackTrunkdEdx(dedx_values);

        float llr_pid = llr_pid_calculator.LLR_many_hits_one_plane(dedx_values_corrected, par_values, plane);

        if (plane == 0)
        {
          _trk_llr_pid_u_v.back() = llr_pid;
          _trk_calo_energy_u_v.back() = calo_energy;
          _trk_trunk_dEdx_u_v.back() = trk_trunk_dEdx;
        }
        else if (plane == 1)
        {
          _trk_llr_pid_v_v.back() = llr_pid;
          _trk_calo_energy_v_v.back() = calo_energy;
          _trk_trunk_dEdx_v_v.back() = trk_trunk_dEdx;
        }
        else if (plane == 2)
        {
          _trk_llr_pid_y_v.back() = llr_pid;
          _trk_calo_energy_y_v.back() = calo_energy;
          _trk_trunk_dEdx_y_v.back() = trk_trunk_dEdx;
        }
        _trk_llr_pid_v.back() += llr_pid;
      }
      _trk_llr_pid_score_v.back() = atan(_trk_llr_pid_v.back() / 100.) * 2 / 3.14159266;
    }
    else
    {
      fillDefault();
    }
  } // for all PFParticles
}

void TrackAnalysis::fillDefault()
{
  _trk_pfp_id_v.push_back(std::numeric_limits<int>::lowest());

  _trk_distance_v.push_back(std::numeric_limits<float>::lowest());

  _trk_theta_v.push_back(std::numeric_limits<float>::lowest());
  _trk_phi_v.push_back(std::numeric_limits<float>::lowest());

  _trk_dir_x_v.push_back(std::numeric_limits<float>::lowest());
  _trk_dir_y_v.push_back(std::numeric_limits<float>::lowest());
  _trk_dir_z_v.push_back(std::numeric_limits<float>::lowest());

  _trk_start_x_v.push_back(std::numeric_limits<float>::lowest());
  _trk_start_y_v.push_back(std::numeric_limits<float>::lowest());
  _trk_start_z_v.push_back(std::numeric_limits<float>::lowest());

  _trk_sce_start_x_v.push_back(std::numeric_limits<float>::lowest());
  _trk_sce_start_y_v.push_back(std::numeric_limits<float>::lowest());
  _trk_sce_start_z_v.push_back(std::numeric_limits<float>::lowest());

  _trk_end_x_v.push_back(std::numeric_limits<float>::lowest());
  _trk_end_y_v.push_back(std::numeric_limits<float>::lowest());
  _trk_end_z_v.push_back(std::numeric_limits<float>::lowest());

  _trk_sce_end_x_v.push_back(std::numeric_limits<float>::lowest());
  _trk_sce_end_y_v.push_back(std::numeric_limits<float>::lowest());
  _trk_sce_end_z_v.push_back(std::numeric_limits<float>::lowest());

  _trk_len_v.push_back(std::numeric_limits<float>::lowest());

  _trk_bragg_p_v.push_back(std::numeric_limits<float>::lowest());
  _trk_bragg_mu_v.push_back(std::numeric_limits<float>::lowest());
  _trk_bragg_pion_v.push_back(std::numeric_limits<float>::lowest());
  _trk_bragg_mip_v.push_back(std::numeric_limits<float>::lowest());
  _trk_pid_chipr_v.push_back(std::numeric_limits<float>::lowest());
  _trk_pid_chika_v.push_back(std::numeric_limits<float>::lowest());
  _trk_pid_chipi_v.push_back(std::numeric_limits<float>::lowest());
  _trk_pid_chimu_v.push_back(std::numeric_limits<float>::lowest());
  _trk_pida_v.push_back(std::numeric_limits<float>::lowest());

  _trk_bragg_p_u_v.push_back(std::numeric_limits<float>::lowest());
  _trk_bragg_mu_u_v.push_back(std::numeric_limits<float>::lowest());
  _trk_bragg_pion_u_v.push_back(std::numeric_limits<float>::lowest());
  _trk_bragg_mip_u_v.push_back(std::numeric_limits<float>::lowest());
  _trk_pid_chipr_u_v.push_back(std::numeric_limits<float>::lowest());
  _trk_pid_chika_u_v.push_back(std::numeric_limits<float>::lowest());
  _trk_pid_chipi_u_v.push_back(std::numeric_limits<float>::lowest());
  _trk_pid_chimu_u_v.push_back(std::numeric_limits<float>::lowest());
  _trk_pida_u_v.push_back(std::numeric_limits<float>::lowest());

  _trk_bragg_p_v_v.push_back(std::numeric_limits<float>::lowest());
  _trk_bragg_mu_v_v.push_back(std::numeric_limits<float>::lowest());
  _trk_bragg_pion_v_v.push_back(std::numeric_limits<float>::lowest());
  _trk_bragg_mip_v_v.push_back(std::numeric_limits<float>::lowest());
  _trk_pid_chipr_v_v.push_back(std::numeric_limits<float>::lowest());
  _trk_pid_chika_v_v.push_back(std::numeric_limits<float>::lowest());
  _trk_pid_chipi_v_v.push_back(std::numeric_limits<float>::lowest());
  _trk_pid_chimu_v_v.push_back(std::numeric_limits<float>::lowest());
  _trk_pida_v_v.push_back(std::numeric_limits<float>::lowest());

  _trk_mcs_muon_mom_v.push_back(std::numeric_limits<float>::lowest());
  _trk_range_muon_mom_v.push_back(std::numeric_limits<float>::lowest());
  _trk_energy_proton_v.push_back(std::numeric_limits<float>::lowest());
  _trk_energy_muon_v.push_back(std::numeric_limits<float>::lowest());
  _trk_calo_energy_u_v.push_back(std::numeric_limits<float>::lowest());
  _trk_calo_energy_v_v.push_back(std::numeric_limits<float>::lowest());
  _trk_calo_energy_y_v.push_back(std::numeric_limits<float>::lowest());

  _trk_llr_pid_u_v.push_back(std::numeric_limits<float>::lowest());
  _trk_llr_pid_v_v.push_back(std::numeric_limits<float>::lowest());
  _trk_llr_pid_y_v.push_back(std::numeric_limits<float>::lowest());
  _trk_llr_pid_v.push_back(std::numeric_limits<float>::lowest());
  _trk_llr_pid_score_v.push_back(std::numeric_limits<float>::lowest());

  _trk_trunk_dEdx_u_v.push_back(std::numeric_limits<float>::lowest());
  _trk_trunk_dEdx_v_v.push_back(std::numeric_limits<float>::lowest());
  _trk_trunk_dEdx_y_v.push_back(std::numeric_limits<float>::lowest());
}

void TrackAnalysis::setBranches(TTree *_tree)
{
  _tree->Branch("trk_bragg_p_v", "std::vector< float >", &_trk_bragg_p_v);
  _tree->Branch("trk_bragg_mu_v", "std::vector< float >", &_trk_bragg_mu_v);
  _tree->Branch("trk_bragg_pion_v", "std::vector< float >", &_trk_bragg_pion_v);
  _tree->Branch("trk_bragg_mip_v", "std::vector< float >", &_trk_bragg_mip_v);
  _tree->Branch("trk_pida_v", "std::vector< float >", &_trk_pida_v);
  _tree->Branch("trk_pid_chipr_v", "std::vector< float >", &_trk_pid_chipr_v);
  _tree->Branch("trk_pid_chipi_v", "std::vector< float >", &_trk_pid_chipi_v);
  _tree->Branch("trk_pid_chika_v", "std::vector< float >", &_trk_pid_chika_v);
  _tree->Branch("trk_pid_chimu_v", "std::vector< float >", &_trk_pid_chimu_v);

  _tree->Branch("trk_bragg_p_u_v", "std::vector< float >", &_trk_bragg_p_u_v);
  _tree->Branch("trk_bragg_mu_u_v", "std::vector< float >", &_trk_bragg_mu_u_v);
  _tree->Branch("trk_bragg_pion_u_v", "std::vector< float >", &_trk_bragg_pion_u_v);
  _tree->Branch("trk_bragg_mip_u_v", "std::vector< float >", &_trk_bragg_mip_u_v);
  _tree->Branch("trk_pida_u_v", "std::vector< float >", &_trk_pida_u_v);
  _tree->Branch("trk_pid_chipr_u_v", "std::vector< float >", &_trk_pid_chipr_u_v);
  _tree->Branch("trk_pid_chipi_u_v", "std::vector< float >", &_trk_pid_chipi_u_v);
  _tree->Branch("trk_pid_chika_u_v", "std::vector< float >", &_trk_pid_chika_u_v);
  _tree->Branch("trk_pid_chimu_u_v", "std::vector< float >", &_trk_pid_chimu_u_v);

  _tree->Branch("trk_bragg_p_v_v", "std::vector< float >", &_trk_bragg_p_v_v);
  _tree->Branch("trk_bragg_mu_v_v", "std::vector< float >", &_trk_bragg_mu_v_v);
  _tree->Branch("trk_bragg_pion_v_v", "std::vector< float >", &_trk_bragg_pion_v_v);
  _tree->Branch("trk_bragg_mip_v_v", "std::vector< float >", &_trk_bragg_mip_v_v);
  _tree->Branch("trk_pida_v_v", "std::vector< float >", &_trk_pida_v_v);
  _tree->Branch("trk_pid_chipr_v_v", "std::vector< float >", &_trk_pid_chipr_v_v);
  _tree->Branch("trk_pid_chipi_v_v", "std::vector< float >", &_trk_pid_chipi_v_v);
  _tree->Branch("trk_pid_chika_v_v", "std::vector< float >", &_trk_pid_chika_v_v);
  _tree->Branch("trk_pid_chimu_v_v", "std::vector< float >", &_trk_pid_chimu_v_v);

  _tree->Branch("trk_pfp_id_v", "std::vector< size_t >", &_trk_pfp_id_v);
  _tree->Branch("trk_dir_x_v", "std::vector< float >", &_trk_dir_x_v);
  _tree->Branch("trk_dir_y_v", "std::vector< float >", &_trk_dir_y_v);
  _tree->Branch("trk_dir_z_v", "std::vector< float >", &_trk_dir_z_v);

  _tree->Branch("trk_start_x_v", "std::vector< float >", &_trk_start_x_v);
  _tree->Branch("trk_start_y_v", "std::vector< float >", &_trk_start_y_v);
  _tree->Branch("trk_start_z_v", "std::vector< float >", &_trk_start_z_v);

  _tree->Branch("trk_sce_start_x_v", "std::vector< float >", &_trk_sce_start_x_v);
  _tree->Branch("trk_sce_start_y_v", "std::vector< float >", &_trk_sce_start_y_v);
  _tree->Branch("trk_sce_start_z_v", "std::vector< float >", &_trk_sce_start_z_v);

  _tree->Branch("trk_end_x_v", "std::vector< float >", &_trk_end_x_v);
  _tree->Branch("trk_end_y_v", "std::vector< float >", &_trk_end_y_v);
  _tree->Branch("trk_end_z_v", "std::vector< float >", &_trk_end_z_v);

  _tree->Branch("trk_sce_end_x_v", "std::vector< float >", &_trk_sce_end_x_v);
  _tree->Branch("trk_sce_end_y_v", "std::vector< float >", &_trk_sce_end_y_v);
  _tree->Branch("trk_sce_end_z_v", "std::vector< float >", &_trk_sce_end_z_v);

  _tree->Branch("trk_distance_v", "std::vector< float >", &_trk_distance_v);
  _tree->Branch("trk_theta_v", "std::vector< float >", &_trk_theta_v);
  _tree->Branch("trk_phi_v", "std::vector< float >", &_trk_phi_v);

  _tree->Branch("trk_len_v", "std::vector< float >", &_trk_len_v);
  _tree->Branch("trk_mcs_muon_mom_v", "std::vector< float >", &_trk_mcs_muon_mom_v);
  _tree->Branch("trk_range_muon_mom_v", "std::vector< float >", &_trk_range_muon_mom_v);
  _tree->Branch("trk_energy_proton_v", "std::vector< float >", &_trk_energy_proton_v);
  _tree->Branch("trk_energy_muon_v", "std::vector< float >", &_trk_energy_muon_v);
  _tree->Branch("trk_calo_energy_u_v", "std::vector< float >", &_trk_calo_energy_u_v);
  _tree->Branch("trk_calo_energy_v_v", "std::vector< float >", &_trk_calo_energy_v_v);
  _tree->Branch("trk_calo_energy_y_v", "std::vector< float >", &_trk_calo_energy_y_v);

  _tree->Branch("trk_llr_pid_u_v", "std::vector<float>", &_trk_llr_pid_u_v);
  _tree->Branch("trk_llr_pid_v_v", "std::vector<float>", &_trk_llr_pid_v_v);
  _tree->Branch("trk_llr_pid_y_v", "std::vector<float>", &_trk_llr_pid_y_v);
  _tree->Branch("trk_llr_pid_v", "std::vector<float>", &_trk_llr_pid_v);
  _tree->Branch("trk_llr_pid_score_v", "std::vector<float>", &_trk_llr_pid_score_v);

  _tree->Branch("trk_trunk_dEdx_u_v", "std::vector<float>", &_trk_trunk_dEdx_u_v);
  _tree->Branch("trk_trunk_dEdx_v_v", "std::vector<float>", &_trk_trunk_dEdx_v_v);
  _tree->Branch("trk_trunk_dEdx_y_v", "std::vector<float>", &_trk_trunk_dEdx_y_v);
}

void TrackAnalysis::resetTTree(TTree *_tree)
{
  _trk_bragg_p_v.clear();
  _trk_bragg_mu_v.clear();
  _trk_bragg_pion_v.clear();
  _trk_bragg_mip_v.clear();
  _trk_pida_v.clear();
  _trk_pid_chipr_v.clear();
  _trk_pid_chika_v.clear();
  _trk_pid_chipi_v.clear();
  _trk_pid_chimu_v.clear();

  _trk_bragg_p_u_v.clear();
  _trk_bragg_mu_u_v.clear();
  _trk_bragg_pion_u_v.clear();
  _trk_bragg_mip_u_v.clear();
  _trk_pida_u_v.clear();
  _trk_pid_chipr_u_v.clear();
  _trk_pid_chika_u_v.clear();
  _trk_pid_chipi_u_v.clear();
  _trk_pid_chimu_u_v.clear();

  _trk_bragg_p_v_v.clear();
  _trk_bragg_mu_v_v.clear();
  _trk_bragg_pion_v_v.clear();
  _trk_bragg_mip_v_v.clear();
  _trk_pida_v_v.clear();
  _trk_pid_chipr_v_v.clear();
  _trk_pid_chika_v_v.clear();
  _trk_pid_chipi_v_v.clear();
  _trk_pid_chimu_v_v.clear();

  _trk_pfp_id_v.clear();

  _trk_start_x_v.clear();
  _trk_start_y_v.clear();
  _trk_start_z_v.clear();

  _trk_sce_start_x_v.clear();
  _trk_sce_start_y_v.clear();
  _trk_sce_start_z_v.clear();

  _trk_end_x_v.clear();
  _trk_end_y_v.clear();
  _trk_end_z_v.clear();

  _trk_sce_end_x_v.clear();
  _trk_sce_end_y_v.clear();
  _trk_sce_end_z_v.clear();

  _trk_dir_x_v.clear();
  _trk_dir_y_v.clear();
  _trk_dir_z_v.clear();
  _trk_distance_v.clear();

  _trk_theta_v.clear();
  _trk_phi_v.clear();

  _trk_len_v.clear();

  _trk_mcs_muon_mom_v.clear();
  _trk_range_muon_mom_v.clear();
  _trk_energy_muon_v.clear();
  _trk_energy_proton_v.clear();
  _trk_calo_energy_u_v.clear();
  _trk_calo_energy_v_v.clear();
  _trk_calo_energy_y_v.clear();

  _trk_llr_pid_u_v.clear();
  _trk_llr_pid_v_v.clear();
  _trk_llr_pid_y_v.clear();
  _trk_llr_pid_v.clear();
  _trk_llr_pid_score_v.clear();

  _trk_trunk_dEdx_u_v.clear();
  _trk_trunk_dEdx_v_v.clear();
  _trk_trunk_dEdx_y_v.clear();
}

float TrackAnalysis::CalculateTrackTrunkdEdx(const std::vector<float> &dEdx_values) {
  
  unsigned int trk_nhits = dEdx_values.size();

  // initial offset of 3 hits
  int firstHitIdx = trk_nhits - 3 - 1;

  // find max hit corresponding to first third of track hits
  int lastHitIdx = trk_nhits - (int)(trk_nhits/3) - 1; 

  // check at least 5 hits remain, otherwise set as invalid for this track (too short)
  if (firstHitIdx - lastHitIdx < 5) {
    return std::numeric_limits<float>::lowest();
  }
  else {
    // loop through hits extracting relevant dE/dx values
    std::vector<float> trk_trunk_dEdx_values;
    trk_trunk_dEdx_values.reserve(firstHitIdx - lastHitIdx); 

    for (int i = trk_nhits - 1; i >= 0; i--) {

      // skip first part of track
      if (i > firstHitIdx) continue;

      trk_trunk_dEdx_values.push_back(dEdx_values[i]);

      // skip last part of track
      if (i < lastHitIdx) break;  
    }

    // calculate mean, median and standard deviation
    // median
    float median;
    std::sort(trk_trunk_dEdx_values.begin(), trk_trunk_dEdx_values.end());
    if (trk_trunk_dEdx_values.size() % 2 == 0) median = 0.5 * (trk_trunk_dEdx_values[trk_trunk_dEdx_values.size()/2 - 1] + trk_trunk_dEdx_values[trk_trunk_dEdx_values.size()/2]);
    else median = trk_trunk_dEdx_values[trk_trunk_dEdx_values.size()/2];

    // mean
    double sum = std::accumulate(std::begin(trk_trunk_dEdx_values), std::end(trk_trunk_dEdx_values), 0.0);
    double m =  sum / trk_trunk_dEdx_values.size();
    
    // standard deviation
    double accum = 0.0;
    std::for_each(std::begin(trk_trunk_dEdx_values), std::end(trk_trunk_dEdx_values), [&](const double d) {accum += (d - m) * (d - m);});
    double stdev = sqrt(accum / (trk_trunk_dEdx_values.size()-1));

    // create trimmed dE/dx vector,  removing any dE/dx greater than 1 standard deviation above the median 
    std::vector<float> trk_trunk_dEdx_values_trimmed;
    trk_trunk_dEdx_values_trimmed.reserve(firstHitIdx - lastHitIdx);
    for (unsigned int i = 0; i < trk_trunk_dEdx_values.size(); i++) {
      if (trk_trunk_dEdx_values[i] <= median + stdev) trk_trunk_dEdx_values_trimmed.push_back(trk_trunk_dEdx_values[i]);
    }

    // calculate mean of trimmed dE/dx vector
    double sum_trimmed = std::accumulate(std::begin(trk_trunk_dEdx_values_trimmed), std::end(trk_trunk_dEdx_values_trimmed), 0.0);
    double trk_dEdx_trunk =  sum_trimmed / trk_trunk_dEdx_values_trimmed.size();

    return trk_dEdx_trunk;
  }
}

DEFINE_ART_CLASS_TOOL(TrackAnalysis)
} // namespace analysis

#endif
