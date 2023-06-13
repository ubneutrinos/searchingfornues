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
#include "lardataobj/RecoBase/SpacePoint.h"

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

  float CalculateTrackTrunkdEdx(const std::vector<float> &dedxPerHit, const std::vector<float> &residualRangePerHit);
  float GetMultiplaneBraggLikelihood(const art::Ptr<anab::ParticleID> pid, const int pdg, const anab::kTrackDir dir, const float yzAngle, const float sin2AngleThreshold, const int nHitsU, const int nHitsV);
  float GetBraggLikelihood(const art::Ptr<anab::ParticleID>& pid, const int& pdg, const geo::View_t& view, const anab::kTrackDir& dir);
  art::Ptr<anab::ParticleID> GetSingleAssociatedFromTrack(const art::Ptr<recob::Track> &track, const art::Event &event);

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
  art::InputTag fCLSproducer;

  bool fRecalibrateHits;
  float fEnergyThresholdForMCHits;
  std::vector<float> fADCtoE; // vector of ADC to # of e- conversion [to be taken from production reco2 fhicl files]
  float fEndSpacepointDistance;

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

  // The Bragg likelihoods using the collection plane when possible and induction planes as fallback
  // Forward (fwd) direction 
  std::vector<float> _trk_bragg_fwd_p_uvy_v; // Proton
  std::vector<float> _trk_bragg_fwd_mu_uvy_v; // Muon
  std::vector<float> _trk_bragg_fwd_pion_uvy_v; // Pion
  std::vector<float> _trk_bragg_fwd_mip_uvy_v; // MIP
  // Backward (bwd) direction
  std::vector<float> _trk_bragg_bwd_p_uvy_v; // Proton
  std::vector<float> _trk_bragg_bwd_mu_uvy_v; // Muon
  std::vector<float> _trk_bragg_bwd_pion_uvy_v; // Pion
  std::vector<float> _trk_bragg_bwd_mip_uvy_v; // MIP

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

  std::vector<int> _trk_nhits_u_v;
  std::vector<int> _trk_nhits_v_v;
  std::vector<int> _trk_nhits_y_v;

  std::vector<float> _trk_avg_deflection_mean_v;
  std::vector<float> _trk_avg_deflection_stdev_v;
  std::vector<float> _trk_avg_deflection_separation_mean_v;

  std::vector<int> _trk_end_spacepoints_v;
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
  fCLSproducer = p.get<art::InputTag>("CLSproducer", "pandora");
  fEnergyThresholdForMCHits = p.get<float>("EnergyThresholdForMCHits", 0.1);
  fRecalibrateHits = p.get<bool>("RecalibrateHits", false);
  fADCtoE = p.get<std::vector<float>>("ADCtoE");
  fEndSpacepointDistance = p.get<float>("EndSpacepointDistance", 5.0);

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

  // get spacepoint information to look for activity at end of tracks
  auto spacePointHandle = e.getValidHandle<std::vector<recob::SpacePoint>>(fCLSproducer);
  std::vector< art::Ptr<recob::SpacePoint> > spacePointCollection;
  for (size_t i_sp = 0; i_sp < spacePointHandle->size(); i_sp++) {
    spacePointCollection.emplace_back(spacePointHandle, i_sp);
  }
  art::ValidHandle<std::vector<recob::PFParticle>> inputPfParticle = e.getValidHandle<std::vector<recob::PFParticle>>(fCLSproducer);

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
      _trk_nhits_u_v.push_back(0);
      _trk_nhits_v_v.push_back(0);
      _trk_nhits_y_v.push_back(0);
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

        float trk_nhits = dedx_values.size();
        float trk_trunk_dEdx = CalculateTrackTrunkdEdx(dedx_values, rr);

        float llr_pid = llr_pid_calculator.LLR_many_hits_one_plane(dedx_values_corrected, par_values, plane);

        if (plane == 0)
        {
          _trk_llr_pid_u_v.back() = llr_pid;
          _trk_calo_energy_u_v.back() = calo_energy;
          _trk_nhits_u_v.back() = trk_nhits;
          _trk_trunk_dEdx_u_v.back() = trk_trunk_dEdx;
        }
        else if (plane == 1)
        {
          _trk_llr_pid_v_v.back() = llr_pid;
          _trk_calo_energy_v_v.back() = calo_energy;
          _trk_nhits_v_v.back() = trk_nhits;
          _trk_trunk_dEdx_v_v.back() = trk_trunk_dEdx;
        }
        else if (plane == 2)
        {
          _trk_llr_pid_y_v.back() = llr_pid;
          _trk_calo_energy_y_v.back() = calo_energy;
          _trk_nhits_y_v.back() = trk_nhits;
          _trk_trunk_dEdx_y_v.back() = trk_trunk_dEdx;
        }
        _trk_llr_pid_v.back() += llr_pid;
      }
      _trk_llr_pid_score_v.back() = atan(_trk_llr_pid_v.back() / 100.) * 2 / 3.14159266;  

      // track end spacepoints
      int nPoints = 0;
      float distSquared = fEndSpacepointDistance*fEndSpacepointDistance;
      TVector3 trkEnd(_trk_end_sce[0], _trk_end_sce[1], _trk_end_sce[2]);
      for (auto &sp : spacePointCollection) {
        float _sp_sce[3];
        searchingfornues::ApplySCECorrectionXYZ(sp->XYZ()[0], sp->XYZ()[1], sp->XYZ()[2], _sp_sce);
        TVector3 spacePoint(_sp_sce[0], _sp_sce[1], _sp_sce[2]);
        if ((trkEnd - spacePoint).Mag2() < distSquared) nPoints++;
      }
      _trk_end_spacepoints_v.push_back(nPoints);

      // Combined planes bragg lieklihoods
      const float sin2AngleThreshold = 0.175;
      const auto yzAngle = std::atan2(_trk_dir_z_v.back(), _trk_dir_y_v.back());

      const float bragg_fwd_p_uvy_v = GetMultiplaneBraggLikelihood(pidpxy_v[0], 2212, anab::kForward, yzAngle, sin2AngleThreshold, _trk_nhits_u_v.back(), _trk_nhits_v_v.back());
      const float bragg_fwd_mu_uvy_v = GetMultiplaneBraggLikelihood(pidpxy_v[0], 13, anab::kForward, yzAngle, sin2AngleThreshold, _trk_nhits_u_v.back(), _trk_nhits_v_v.back());
      const float bragg_fwd_pion_uvy_v = GetMultiplaneBraggLikelihood(pidpxy_v[0], 211, anab::kForward, yzAngle, sin2AngleThreshold, _trk_nhits_u_v.back(), _trk_nhits_v_v.back());
      const float bragg_fwd_mip_uvy_v = GetMultiplaneBraggLikelihood(pidpxy_v[0], 0, anab::kForward, yzAngle, sin2AngleThreshold, _trk_nhits_u_v.back(), _trk_nhits_v_v.back());

      const float bragg_bwd_p_uvy_v = GetMultiplaneBraggLikelihood(pidpxy_v[0], 2212, anab::kBackward, yzAngle, sin2AngleThreshold, _trk_nhits_u_v.back(), _trk_nhits_v_v.back());
      const float bragg_bwd_mu_uvy_v = GetMultiplaneBraggLikelihood(pidpxy_v[0], 13, anab::kBackward, yzAngle, sin2AngleThreshold, _trk_nhits_u_v.back(), _trk_nhits_v_v.back());
      const float bragg_bwd_pion_uvy_v = GetMultiplaneBraggLikelihood(pidpxy_v[0], 211, anab::kBackward, yzAngle, sin2AngleThreshold, _trk_nhits_u_v.back(), _trk_nhits_v_v.back());
      const float bragg_bwd_mip_uvy_v = GetMultiplaneBraggLikelihood(pidpxy_v[0], 0, anab::kBackward, yzAngle, sin2AngleThreshold, _trk_nhits_u_v.back(), _trk_nhits_v_v.back());


      _trk_bragg_fwd_p_uvy_v.push_back(bragg_fwd_p_uvy_v);
      _trk_bragg_fwd_mu_uvy_v.push_back(bragg_fwd_mu_uvy_v);
      _trk_bragg_fwd_pion_uvy_v.push_back(bragg_fwd_pion_uvy_v);
      _trk_bragg_fwd_mip_uvy_v.push_back(bragg_fwd_mip_uvy_v);
      _trk_bragg_bwd_p_uvy_v.push_back(bragg_bwd_p_uvy_v);
      _trk_bragg_bwd_mu_uvy_v.push_back(bragg_bwd_mu_uvy_v);
      _trk_bragg_bwd_pion_uvy_v.push_back(bragg_bwd_pion_uvy_v);
      _trk_bragg_bwd_mip_uvy_v.push_back(bragg_bwd_mip_uvy_v);

      // Average of deflections along track ("wiggliness")
      // Get track valid points
      std::vector<size_t> validPoints;
      auto firstValidPoint = trk->FirstValidPoint();
      validPoints.push_back(firstValidPoint);
      auto nextValidPoint = trk->NextValidPoint(firstValidPoint + 1);
      while (nextValidPoint != recob::TrackTrajectory::InvalidIndex)
      {
          validPoints.push_back(nextValidPoint);
          nextValidPoint = trk->NextValidPoint(nextValidPoint + 1);
      }
      // determine average deflection between sequential valid points
      if (validPoints.size() < 3) 
      {
        _trk_avg_deflection_mean_v.push_back(0);
        _trk_avg_deflection_stdev_v.push_back(0);
        _trk_avg_deflection_separation_mean_v.push_back(0);
      }
      else 
      {
        std::vector<float> thetaVector;
        float thetaSum = 0.f;
        float separationSum = 0.f;
        for (unsigned int i = 1; i < validPoints.size(); ++i) {
          auto dir = trk->DirectionAtPoint(validPoints.at(i));
          auto dirPrev = trk->DirectionAtPoint(validPoints.at(i - 1));

          // Bind between -1 and 1 at floating precision to avoid issues with cast from double
          auto cosTheta = std::min(1.f, std::max(-1.f, static_cast<float>(dir.Dot(dirPrev))));
          auto theta = std::acos(cosTheta);

          thetaSum += theta;
          thetaVector.push_back(theta);

          // separation between points
          TVector3 point(trk->LocationAtPoint(validPoints.at(i)).X(), trk->LocationAtPoint(validPoints.at(i)).Y(), trk->LocationAtPoint(validPoints.at(i)).Z());
          TVector3 pointPrev(trk->LocationAtPoint(validPoints.at(i - 1)).X(), trk->LocationAtPoint(validPoints.at(i - 1)).Y(), trk->LocationAtPoint(validPoints.at(i - 1)).Z());
          TVector3 separation = point - pointPrev; 
          separationSum += separation.Mag();
        }

        float thetaMean = thetaSum / static_cast<float>(thetaVector.size());
        float separationMean = separationSum / static_cast<float>(thetaVector.size());

        float thetaDiffSum = 0.f;
        for (const auto &theta : thetaVector) 
        {
          thetaDiffSum += std::pow(theta - thetaMean, 2);
        }

        auto variance = thetaDiffSum / static_cast<float>(thetaVector.size() - 1);

        _trk_avg_deflection_mean_v.push_back(thetaMean);
        _trk_avg_deflection_stdev_v.push_back(std::sqrt(variance));
        _trk_avg_deflection_separation_mean_v.push_back(separationMean);
      }
      
      // count number of spacepoints at end of track
      int nPoints = 0;
      const auto distSquared = fEndSpacepointDistance*fEndSpacepointDistance;
      TVector3 trkEnd(_trk_end_sce[0], _trk_end_sce[1], _trk_end_sce[2]);
      for (auto &sp : spacePointCollection) {
        float _sp_sce[3];
        searchingfornues::ApplySCECorrectionXYZ(sp.XYZ()[0], sp.XYZ()[1], sp.XYZ()[2], _sp_sce);
        const TVector3 spacePoint(_sp_sce[0], _sp_sce[1], _sp_sce[2]);
        if ((trkEnd - spacePoint).Mag2() < distSquared) nPoints++;
      }
      _trk_end_spacepoints_v.push_back(nPoints);
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

  _trk_bragg_fwd_p_uvy_v.push_back(std::numeric_limits<float>::lowest());
  _trk_bragg_fwd_mu_uvy_v.push_back(std::numeric_limits<float>::lowest());
  _trk_bragg_fwd_pion_uvy_v.push_back(std::numeric_limits<float>::lowest());
  _trk_bragg_fwd_mip_uvy_v.push_back(std::numeric_limits<float>::lowest());
  _trk_bragg_bwd_p_uvy_v.push_back(std::numeric_limits<float>::lowest());
  _trk_bragg_bwd_mu_uvy_v.push_back(std::numeric_limits<float>::lowest());
  _trk_bragg_bwd_pion_uvy_v.push_back(std::numeric_limits<float>::lowest());
  _trk_bragg_bwd_mip_uvy_v.push_back(std::numeric_limits<float>::lowest());

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

  _trk_nhits_u_v.push_back(std::numeric_limits<int>::lowest());
  _trk_nhits_v_v.push_back(std::numeric_limits<int>::lowest());
  _trk_nhits_y_v.push_back(std::numeric_limits<int>::lowest());

  _trk_avg_deflection_mean_v.push_back(std::numeric_limits<float>::lowest());
  _trk_avg_deflection_stdev_v.push_back(std::numeric_limits<float>::lowest());
  _trk_avg_deflection_separation_mean_v.push_back(std::numeric_limits<float>::lowest());

  _trk_end_spacepoints_v.push_back(std::numeric_limits<int>::lowest());


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

  _tree->Branch("trk_bragg_fwd_p_uvy_v", "std::vector< float >", &_trk_bragg_fwd_p_uvy_v);
  _tree->Branch("trk_bragg_fwd_mu_uvy_v", "std::vector< float >", &_trk_bragg_fwd_mu_uvy_v);
  _tree->Branch("trk_bragg_fwd_pion_uvy_v", "std::vector< float >", &_trk_bragg_fwd_pion_uvy_v);
  _tree->Branch("trk_bragg_fwd_mip_uvy_v", "std::vector< float >", &_trk_bragg_fwd_mip_uvy_v);
  _tree->Branch("trk_bragg_bwd_p_uvy_v", "std::vector< float >", &_trk_bragg_bwd_p_uvy_v);
  _tree->Branch("trk_bragg_bwd_mu_uvy_v", "std::vector< float >", &_trk_bragg_bwd_mu_uvy_v);
  _tree->Branch("trk_bragg_bwd_pion_uvy_v", "std::vector< float >", &_trk_bragg_bwd_pion_uvy_v);
  _tree->Branch("trk_bragg_bwd_mip_uvy_v", "std::vector< float >", &_trk_bragg_bwd_mip_uvy_v);

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

  _tree->Branch("trk_nhits_u_v", "std::vector<int>", &_trk_nhits_u_v);
  _tree->Branch("trk_nhits_v_v", "std::vector<int>", &_trk_nhits_v_v);
  _tree->Branch("trk_nhits_y_v", "std::vector<int>", &_trk_nhits_y_v);

  _tree->Branch("trk_avg_deflection_mean_v", "std::vector<float>", &_trk_avg_deflection_mean_v);
  _tree->Branch("trk_avg_deflection_stdev_v", "std::vector<float>", &_trk_avg_deflection_stdev_v);
  _tree->Branch("trk_avg_deflection_separation_mean_v", "std::vector<float>", &_trk_avg_deflection_separation_mean_v);

  _tree->Branch("trk_end_spacepoints_v", "std::vector<int>", &_trk_end_spacepoints_v);
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

  _trk_bragg_fwd_p_uvy_v.clear();
  _trk_bragg_fwd_mu_uvy_v.clear();
  _trk_bragg_fwd_pion_uvy_v.clear();
  _trk_bragg_fwd_mip_uvy_v.clear();
  _trk_bragg_bwd_p_uvy_v.clear();
  _trk_bragg_bwd_mu_uvy_v.clear();
  _trk_bragg_bwd_pion_uvy_v.clear();
  _trk_bragg_bwd_mip_uvy_v.clear();

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

  _trk_nhits_u_v.clear();
  _trk_nhits_v_v.clear();
  _trk_nhits_y_v.clear();

  _trk_avg_deflection_mean_v.clear();
  _trk_avg_deflection_stdev_v.clear();
  _trk_avg_deflection_separation_mean_v.clear();

  _trk_end_spacepoints_v.clear();
}

float TrackAnalysis::CalculateTrackTrunkdEdx(const std::vector<float> &dedxPerHit, const std::vector<float> &residualRangePerHit)
{
    const auto nHitsToSkip = 3u;
    const auto lengthFraction = 1.f/3;
    if (dedxPerHit.size() != residualRangePerHit.size())
        throw cet::exception("RecoHelper::GetTruncatedMeandEdxAtTrackStart") << " - dEdx per hit and residual range vectors have different sizes" << std::endl;

    // Check if the variable is calculable
    if (residualRangePerHit.size() <= nHitsToSkip)
        return -std::numeric_limits<float>::max();

    // Make the vector of pairs to keep track of the incides, and find the maximum residual range
    std::vector<std::pair<float, unsigned int> > residualRangeIndices;
    float maxResidualRange = -std::numeric_limits<float>::max();
    for (unsigned int i = 0; i < residualRangePerHit.size(); ++i)
    {
        const auto residualRange = residualRangePerHit.at(i);
        maxResidualRange = std::max(maxResidualRange, residualRange);

        residualRangeIndices.emplace_back(residualRange, i);
    }

    const auto residualRangeCutoff = maxResidualRange * lengthFraction;

    // Sort the residual ranges such that the largest residual range (closest to the start of the track) is first
    std::sort(residualRangeIndices.begin(), residualRangeIndices.end(), [](auto &a, auto &b) {
        return a.first > b.first;
    });

    // Get the dEdx of the hits at the start of the track
    std::vector<float> dedxPerHitAtStart;
    for (unsigned int i = nHitsToSkip; i < residualRangeIndices.size(); ++i)
    {
        const auto entry = residualRangeIndices.at(i);
        const auto residualRange = entry.first;
        const auto hitIndex = entry.second;

        // ATTN small residual ranges are at the start of the track
        if (residualRange < residualRangeCutoff)
            continue;

        dedxPerHitAtStart.push_back(dedxPerHit.at(hitIndex));
    }

    const auto nHits = dedxPerHitAtStart.size();
    if (nHits == 0)
        return -std::numeric_limits<float>::max();

    // Sort the dEdx so we can find the median
    std::sort(dedxPerHitAtStart.begin(), dedxPerHitAtStart.end());
    const auto median = dedxPerHitAtStart.at(nHits / 2);

    //  #### Find the mean #### 
    float total = 0.f;
    for (const auto &dEdx : dedxPerHitAtStart)
        total += dEdx;
    const auto mean = total / static_cast<float>(nHits);

    // #### Find the variance ####
    float squareSum = 0.f;
    for (const auto &dEdx : dedxPerHitAtStart)
        squareSum += std::pow(dEdx - mean, 2);
    const auto variance = squareSum / static_cast<float>(nHits);

    // Get the mean dEdx of the hits within one standard deviation of the median
    float truncatedTotal = 0.f;
    unsigned int nTruncatedHits = 0;
    for (const auto &dEdx : dedxPerHitAtStart)
    {
        if (std::pow(dEdx - median, 2) > variance)
            continue;

        truncatedTotal += dEdx;
        nTruncatedHits++;
    }

    if (nTruncatedHits == 0)
        return -std::numeric_limits<float>::max();

    return truncatedTotal / static_cast<float>(nTruncatedHits);
}



float TrackAnalysis::CalculateTrackTrunkdEdx2(const std::vector<float> &dedxPerHit, const std::vector<float> &residualRangePerHit)
{
    const auto nHitsToSkip = 3u;
    const auto lengthFraction = 0.333333333f;
    if (dedxPerHit.size() != residualRangePerHit.size())
        throw cet::exception("RecoHelper::GetTruncatedMeandEdxAtTrackStart") << " - dEdx per hit and residual range vectors have different sizes" << std::endl;

    // Check if the variable is calculable
    if (residualRangePerHit.size() <= nHitsToSkip)
        return -std::numeric_limits<float>::max();

    // Make the vector of pairs to keep track of the incides, and find the maximum residual range
    std::vector<std::pair<float, unsigned int> > residualRangeIndices;
    float maxResidualRange = -std::numeric_limits<float>::max();
    for (unsigned int i = 0; i < residualRangePerHit.size(); ++i)
    {
        const auto residualRange = residualRangePerHit.at(i);
        maxResidualRange = std::max(maxResidualRange, residualRange);

        residualRangeIndices.emplace_back(residualRange, i);
    }

    const auto residualRangeCutoff = maxResidualRange * lengthFraction;

    // Sort the residual ranges such that the largest residual range (closest to the start of the track) is first
    std::sort(residualRangeIndices.begin(), residualRangeIndices.end(), [](auto &a, auto &b) {
        return a.first > b.first;
    });

    // Get the dEdx of the hits at the start of the track
    std::vector<float> dedxPerHitAtStart;
    for (unsigned int i = nHitsToSkip; i < residualRangeIndices.size(); ++i)
    {
        const auto entry = residualRangeIndices.at(i);
        const auto residualRange = entry.first;
        const auto hitIndex = entry.second;

        // ATTN small residual ranges are at the start of the track
        if (residualRange < residualRangeCutoff)
            continue;

        dedxPerHitAtStart.push_back(dedxPerHit.at(hitIndex));
    }

    const auto nHits = dedxPerHitAtStart.size();
    if (nHits == 0)
        return -std::numeric_limits<float>::max();

    // Sort the dEdx so we can find the median
    std::sort(dedxPerHitAtStart.begin(), dedxPerHitAtStart.end());
    const auto median = dedxPerHitAtStart.at(nHits / 2);

    //  #### Find the mean #### 
    float total = 0.f;
    for (const auto &dEdx : dedxPerHitAtStart)
        total += dEdx;
    const auto mean = total / static_cast<float>(nHits);

    // #### Find the variance ####
    float squareSum = 0.f;
    for (const auto &dEdx : dedxPerHitAtStart)
        squareSum += std::pow(dEdx - mean, 2);
    const auto variance = squareSum / static_cast<float>(nHits);

    // Get the mean dEdx of the hits within one standard deviation of the median
    float truncatedTotal = 0.f;
    unsigned int nTruncatedHits = 0;
    for (const auto &dEdx : dedxPerHitAtStart)
    {
        if (std::pow(dEdx - median, 2) > variance)
            continue;

        truncatedTotal += dEdx;
        nTruncatedHits++;
    }

    if (nTruncatedHits == 0)
        return -std::numeric_limits<float>::max();

    return truncatedTotal / static_cast<float>(nTruncatedHits);
}





float TrackAnalysis::GetBraggLikelihood(const art::Ptr<anab::ParticleID>& pid, const int& pdg, const geo::View_t& view, const anab::kTrackDir& dir) {
    bool found = false;
    float score = -std::numeric_limits<float>::max();

    for (const auto& algo : pid->ParticleIDAlgScores()) {
        const std::bitset<8> planeMask(algo.fPlaneMask);

        const bool usesW = planeMask.test(2);
        const bool usesU = planeMask.test(0);
        const bool usesV = planeMask.test(1);

        if (algo.fAlgName == "BraggPeakLLH" && algo.fTrackDir == dir && algo.fAssumedPdg == pdg && 
          ((usesW && !usesU && !usesV && view == geo::kW) || (usesU && !usesW && !usesV && view == geo::kU) || (usesV && !usesU && !usesW && view == geo::kV))) 
        {
            if (found) {
                throw cet::exception("RecoHelper::GetBraggLikelihood") << " - Ambiguous criteria supplied." << std::endl;
            }

            found = true;
            score = algo.fValue;
        }
    }

    return score;
}


  // double PID(art::Ptr<anab::ParticleID> selected_pid, std::string AlgName, anab::kVariableType VariableType, anab::kTrackDir TrackDirection, int pdgCode, int selectedPlane)
  // {
  //   std::vector<anab::sParticleIDAlgScores> AlgScoresVec = selected_pid->ParticleIDAlgScores();
  //   for (size_t i_algscore = 0; i_algscore < AlgScoresVec.size(); i_algscore++)
  //   {
  //     anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
  //     int planeid = UBPID::uB_getSinglePlane(AlgScore.fPlaneMask);

  //     if (selectedPlane == 0 || selectedPlane == 1 || selectedPlane == 2)
  //     {
  //       if (selectedPlane != planeid) continue;
  //     }

  //     if (AlgScore.fAlgName == AlgName)
  //     {
  //       if (anab::kVariableType(AlgScore.fVariableType) == VariableType && anab::kTrackDir(AlgScore.fTrackDir) == TrackDirection) // diff: VariableType!!!, kTrackdir ?
  //       {
  //         if (AlgScore.fAssumedPdg == pdgCode)
  //         {
  //           double alg_value = AlgScore.fValue;
  //           return alg_value; // diff: returns first matching value
  //         }
  //       }
  //     }
  //   }
  //   return std::numeric_limits<double>::lowest();
  // }

// Definition of a function that returns an art pointer to an anab::ParticleID object that is associated to a single input track
art::Ptr<anab::ParticleID> TrackAnalysis::GetSingleAssociatedFromTrack(const art::Ptr<recob::Track> &track, const art::Event &event) {
    // Definition of labels for the tracks and the particle identification object (PID) labels
    std::cout<<"DEBUG GSAFT Point 0"<<std::endl;
    art::InputTag TrackLabel("pandoraTrack");
    art::InputTag PIDLabel("pandoraTrackcalipid");
    // const std::string TrackLabel = "passsndora";
    // const std::string PIDLabel = "pandoracalipidSCE";
    std::cout<<"DEBUG GSAFT Point 1"<<std::endl;

    // Definition of an unordered map that associates each track with a vector of PIDs
    using Association = std::unordered_map< art::Ptr<recob::Track>, std::vector< art::Ptr<anab::ParticleID> > >;
    Association trackToPIDs;
    std::cout<<"DEBUG GSAFT Point 2"<<std::endl;

    // Retrieval of a handle to the track collection and the PID association
    const auto handle = event.getValidHandle< std::vector<recob::Track> >(TrackLabel);
    std::cout<<"DEBUG GSAFT Point 3"<<std::endl;
    const art::FindManyP<anab::ParticleID> assoc(handle, event, PIDLabel);
    std::cout<<"DEBUG GSAFT Point 4"<<std::endl;

    // Loop over all the tracks in the track collection and associate them with their corresponding PIDs
    for (unsigned int i = 0; i < handle->size(); ++i) {
        const art::Ptr<recob::Track> objectL(handle, i);
        std::cout<<"DEBUG GSAFT Point 4.1"<<std::endl;
        for (const auto &objectR : assoc.at(objectL.key())) {
            trackToPIDs[objectL].push_back(objectR);
            std::cout<<"DEBUG GSAFT Point 4.2"<<std::endl;
        }
    }
    std::cout<<"DEBUG GSAFT Point 5 - trackToPIDs size: "<<trackToPIDs.size()<<std::endl;
    // Todo remove
    for (const auto& kv : trackToPIDs) {
    // Print out the key (track) and value (vector of ParticleID objects)
    std::cout << "DEBUG trackToPIDs - Track " << kv.first.key() << " - length: "<<kv.first->Length()<< ":\n";
    for (const auto& pid : kv.second) {
        for (const auto& algo : pid->ParticleIDAlgScores()) {
            std::cout<<"DEBUG trackToPIDs algo.fAlgName: " << algo.fAlgName << " - algo.fValue: " << algo.fValue << " - algo.fVariableType: " << algo.fVariableType << " - algo.fTrackDir: " << algo.fTrackDir << " - algo.fAssumedPdg: " << algo.fAssumedPdg << " - algo.fPlaneMask: " << algo.fPlaneMask << std::endl;
            break;
        }
    }
    }
    std::cout<<std::endl;


    // Find the track associated with the input track and return its PID pointer
    const auto iter = trackToPIDs.find(track);
    std::cout<<"DEBUG GSAFT Point 6"<<std::endl;
    if (iter == trackToPIDs.end()) {
        return art::Ptr<anab::ParticleID>();
    }
    std::cout<<"DEBUG GSAFT Point 7"<<std::endl;
    const auto& objects = iter->second;
    if (objects.size() != 1) {
        throw cet::exception("GetSingleAssociatedFromTrack") << " - Found " << objects.size() << " objects associated to the input track. Expected 1." << std::endl;
    }
    std::cout<<"DEBUG GSAFT Point 8"<<std::endl;
    return objects.front();
}

// art::Ptr<anab::ParticleID> TrackAnalysis::GetSingleAssociatedFromTrack(const art::Ptr<recob::Track> &track, const art::Event &event)
// {
//     // template <typename L, typename R>
//     // using Association = std::unordered_map< art::Ptr<L>, Collection<R> >;  ///< An association from L->R
//     // template <typename T>
//     // using Collection = std::vector< art::Ptr<T> >;  ///< A collection of objects of type T

//     const std::string trackLabel = "pandora";
//     const std::string pidLabel = "pandoracalipidSCE";

//     std::unordered_map< art::Ptr<recob::Track>, std::vector< art::Ptr<anab::ParticleID> > > trackToPIDs;
//     const auto handle = event.getValidHandle< std::vector<recob::Track> >(trackLabel);
//     const art::FindManyP<anab::ParticleID> assoc(handle, event, pidLabel);

//     for (unsigned int i = 0; i < handle->size(); ++i)getpid
//     {
//         const art::Ptr<recob::Track> objectL(handle, i);

//         for (const auto &objectR : assoc.at(objectL.key()))
//             trackToPIDs[objectL].push_back(objectR);
//     }

//     const auto iter = trackToPIDs.find(track);
//     const auto objects = iter == trackToPIDs.end() ? std::vector< art::Ptr<anab::ParticleID> >() : iter->second;

//     if (objects.size() != 1)
//         throw cet::exception("TrackAnalysis::GetSingleAssociatedFromTrack") << " - Found " << objects.size() << " objects associated to the input object. Expected 1." << std::endl;

//     return objects.front();
// }

float TrackAnalysis::GetMultiplaneBraggLikelihood(const art::Ptr<anab::ParticleID> pid, const int pdg, const anab::kTrackDir dir, const float yzAngle, const float sin2AngleThreshold, const int nHitsU, const int nHitsV)
{
  // When no PID object is provided:
  if (!pid) return -std::numeric_limits<float>::max();

  const auto piBy3 = std::acos(0.5f);
  const auto likelihoodW = searchingfornues::PID(pid, "BraggPeakLLH", anab::kLikelihood, dir, pdg, 2);
  std::cout<<"DEBUG GetMultiplaneBraggLikelihood Point -1"<<std::endl;
  const auto likelihoodWTest = TrackAnalysis::GetBraggLikelihood(pid, pdg, geo::kW, dir);
  std::cout<<"DEBUG GetMultiplaneBraggLikelihood: likelihoodW = "<<likelihoodW<<", likelihoodWTest = "<<likelihoodWTest<<std::endl;
  const bool isTrackAlongWWire = (std::pow(std::sin(yzAngle), 2) < sin2AngleThreshold);
  const auto hasW = (likelihoodW > -1.f && !isTrackAlongWWire);

  std::cout<<"DEBUG GetMultiplaneBraggLikelihood  Point 0"<<std::endl;
  if (hasW)
      return likelihoodW;
  std::cout<<"DEBUG GetMultiplaneBraggLikelihood  Point 1"<<std::endl;

  // Otherwise the track is along the W wire direction, so just average the other two planes weighted by the number of degrees of freedom
  const auto likelihoodU = searchingfornues::PID(pid, "BraggPeakLLH", anab::kLikelihood, dir, pdg, 0);
  const auto likelihoodUTest = TrackAnalysis::GetBraggLikelihood(pid, pdg, geo::kU, dir);
  std::cout<<"\tDEBUG GetMultiplaneBraggLikelihood: likelihoodU = "<<likelihoodU<<", likelihoodUTest = "<<likelihoodUTest<<std::endl;
  const bool isTrackAlongUWire = (std::pow(std::sin(yzAngle - piBy3), 2) < sin2AngleThreshold);
  const auto hasU = (likelihoodU > -1.f && !isTrackAlongUWire);

  const auto likelihoodV = searchingfornues::PID(pid, "BraggPeakLLH", anab::kLikelihood, dir, pdg, 1);
  const auto likelihoodVTest = TrackAnalysis::GetBraggLikelihood(pid, pdg, geo::kV, dir);
  std::cout<<"\tDEBUG GetMultiplaneBraggLikelihood: likelihoodV = "<<likelihoodV<<", likelihoodVTest = "<<likelihoodVTest<<std::endl;
  const bool isTrackAlongVWire = (std::pow(std::sin(yzAngle + piBy3), 2) < sin2AngleThreshold);
  const auto hasV = (likelihoodV > -1.f && !isTrackAlongVWire);

  const auto dofU = hasU ? static_cast<float>(nHitsU) : 0.f;
  const auto dofV = hasV ? static_cast<float>(nHitsV) : 0.f;
  const auto dofUV = dofU + dofV;

  const auto weightU = hasU ? (likelihoodU * dofU) : 0.f;
  const auto weightV = hasV ? (likelihoodV * dofV) : 0.f;

  if ((!hasU && !hasV) || (nHitsU == 0 && nHitsV == 0) || dofUV <= std::numeric_limits<float>::epsilon())
      return -std::numeric_limits<float>::max();

  const auto likelihood = (weightU + weightV) / dofUV;
  if(likelihood > -1.f) return likelihood;
  else return -std::numeric_limits<float>::max();
}

DEFINE_ART_CLASS_TOOL(TrackAnalysis)
} // namespace analysis

#endif